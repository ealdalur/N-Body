# Thread Pool Implementation for Parallel N-Body Computation

## The Problem: Per-Frame Thread Creation Overhead

Gravitational acceleration calculations are embarrassingly parallel — each body's acceleration is independent of the others (given a fixed tree). The naive approach spawns N threads per frame:

```cpp
// Naive: create and destroy threads every frame
for (int i = 0; i < numThreads; i++) {
    threads[i] = std::thread(CalcAccelRange, iStart, iEnd);
}
for (int i = 0; i < numThreads; i++) {
    threads[i].join();
}
```

Thread creation on Windows involves kernel calls, stack allocation, and scheduler registration — roughly 10-50 microseconds per thread. At 14 threads × 6 dispatch points per frame × 100 FPS, that's 8400 thread creates/destroys per second, wasting ~0.5-1ms per frame on overhead alone.

A thread pool creates threads once at startup and reuses them for the lifetime of the program.

## Architecture Overview

```
Main Thread                    Worker Threads (created once)
    |                              |  |  |  ... |
    |-- submit(task1) ------------>|  |  |      |
    |-- submit(task2) -------------->  |  |      |
    |-- submit(task3) ----------------->  |      |
    |-- submit(task4) -------------------->      |
    |                              |  |  |      |
    |-- waitAll() -----> blocks    |  |  |      |
    |                    until     |  executing  |
    |                    all done  |  tasks...   |
    |<---- notify ---- last task finishes        |
    |                              |  |  |      |
    | (workers go back to sleep,   |  |  |      |
    |  waiting for next submit)    |  |  |      |
```

## Implementation

### Data Members

```cpp
class ThreadPool
{
    std::vector<std::thread> workers;           // persistent OS threads
    std::vector<std::function<void()>> tasks;   // pending work queue
    std::mutex mtx;                             // protects tasks + pending_tasks
    std::condition_variable cv_work;            // "wake up, there's work"
    std::condition_variable cv_done;            // "all work is finished"
    int pending_tasks = 0;                      // count of submitted but not-yet-completed tasks
    bool stop = false;                          // shutdown signal
};
```

### Mutexes and Condition Variables

A **mutex** (mutual exclusion) ensures only one thread can access shared state at a time. Without it, two threads reading/writing `tasks` or `pending_tasks` simultaneously would produce data races (undefined behavior).

A **condition variable** enables efficient waiting. Without it, a thread waiting for work would have to spin in a loop:

```cpp
// Bad: busy-wait burns CPU cycles doing nothing
while (tasks.empty()) { /* spin */ }
```

With a condition variable, the thread sleeps (consuming zero CPU) and is woken only when the condition might have changed:

```cpp
// Good: thread sleeps until notified
cv_work.wait(lock, [&]() { return !tasks.empty(); });
```

### The Worker Loop

Each worker thread runs this loop for the lifetime of the pool:

```cpp
while (true) {
    std::function<void()> task;
    {
        std::unique_lock<std::mutex> lock(mtx);
        cv_work.wait(lock, [this]() { return stop || !tasks.empty(); });
        if (stop && tasks.empty()) return;
        task = std::move(tasks.back());
        tasks.pop_back();
    }
    // Mutex is released here — task executes without holding the lock
    task();
    {
        std::lock_guard<std::mutex> lock(mtx);
        if (--pending_tasks == 0) {
            cv_done.notify_one();
        }
    }
}
```

**Step by step:**

1. **Acquire lock.** The worker takes ownership of the mutex, guaranteeing exclusive access to `tasks`.

2. **Wait for work.** `cv_work.wait(lock, predicate)` does the following atomically:
   - Evaluates the predicate (`stop || !tasks.empty()`)
   - If true, returns immediately (still holding the lock)
   - If false, releases the mutex and puts the thread to sleep
   - When another thread calls `cv_work.notify_one()`, the sleeping thread wakes, reacquires the mutex, and re-checks the predicate
   - If predicate is now true, returns; otherwise sleeps again (spurious wakeup protection)

3. **Check shutdown.** If `stop` is set and no tasks remain, the thread exits (used during destructor).

4. **Take a task.** Move the task out of the queue (O(1) from back of vector) and release the lock.

5. **Execute.** The task runs without the mutex held — this is critical for parallelism. If we held the lock during execution, only one task could run at a time.

6. **Signal completion.** Reacquire the lock, decrement `pending_tasks`. If it reaches zero, notify the main thread via `cv_done`.

### submit()

```cpp
void submit(F&& f)
{
    {
        std::lock_guard<std::mutex> lock(mtx);
        tasks.emplace_back(std::forward<F>(f));
        pending_tasks++;
    }
    cv_work.notify_one();
}
```

1. **Lock.** Protect the task queue and counter.
2. **Enqueue.** Add the callable to the queue and increment the pending counter.
3. **Unlock.** The `lock_guard` destructor releases the mutex.
4. **Notify.** Wake one sleeping worker to pick up the task.

Note: `notify_one()` is called **after** releasing the lock. This is fine here — it's an optimization to avoid "hurry up and wait" where a woken thread immediately blocks on the mutex that the notifier still holds.

### waitAll()

```cpp
void waitAll()
{
    std::unique_lock<std::mutex> lock(mtx);
    cv_done.wait(lock, [this]() { return pending_tasks == 0; });
}
```

The main thread blocks until all submitted tasks have completed. The predicate `pending_tasks == 0` is re-checked after every notification, protecting against spurious wakeups.

## The Critical Design Decision: Why pending_tasks Must Be Mutex-Protected

An earlier version used `std::atomic<int>` for the task counter:

```cpp
// BROKEN on hybrid CPUs:
std::atomic<int> active_tasks{0};

task();
if (--active_tasks == 0) {    // atomic decrement, no lock
    cv_done.notify_one();      // notification outside mutex
}
```

This has a race condition. The issue is that `condition_variable::wait()` internally does two things that are only atomic **with respect to the mutex**:

1. Check the predicate
2. Release the mutex and sleep

If the decrement and notify happen without the mutex, there's a window:

```
Main thread:   checks predicate → pending_tasks is 1 → decides to sleep
Worker thread: decrements pending_tasks to 0
Worker thread: calls notify_one() → nobody is sleeping yet!
Main thread:   releases mutex and sleeps → never woken
```

With the mutex protecting the decrement:

```
Main thread:   holds mutex, checks predicate → pending_tasks is 1 → releases mutex and sleeps
Worker thread: acquires mutex, decrements to 0, calls notify_one()
Main thread:   wakes up, reacquires mutex, checks predicate → 0 → returns
```

The mutex creates a **happens-before** relationship: the decrement is guaranteed to be visible to the waiter because they both go through the same mutex. The `condition_variable` specification only guarantees notification delivery with respect to its associated mutex — any state changes that bypass the mutex are invisible to the wait/notify protocol.

This race manifested specifically on Intel Core Ultra hybrid CPUs (P-core + E-core) where cross-core visibility latency widened the timing window enough to trigger the deadlock reliably.

## Usage Pattern in the Simulation

The simulation dispatches work in a fork-join pattern:

```cpp
// Fork: submit work chunks to all threads
int chunk_size = N_BODIES / numThreads;
for (int i = 0; i < numThreads; i++) {
    int iStart = i * chunk_size;
    int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
    pool->submit([this, iStart, iEnd]() {
        CalcAccelRangeOct(iStart, iEnd);
    });
}

// Join: wait for all chunks to complete
pool->waitAll();
```

Each frame has this sequence of fork-join points:
1. **LeapFrog positions** — update all body positions (trivial O(N), but parallelized)
2. **Prepare derivatives** — set up pointer aliases for the state arrays
3. **Compute accelerations** — the expensive tree traversal (dominant cost)
4. **LeapFrog velocities + outputs** — update velocities and compute diagnostics

Steps 1 and 4 are fused where possible to reduce the number of synchronization barriers.

## Thread Count Selection

```cpp
numThreads = hardware_concurrency - 2;
if (numThreads < 4) numThreads = 4;
```

We reserve 2 cores: one for the main thread (which does serial tree building and rendering) and one for the OS/driver. On hybrid architectures, this also avoids over-subscribing when efficiency cores have lower throughput than performance cores — equal-sized work chunks assigned to unequal cores means the fastest core waits idle for the slowest.

## Shutdown

```cpp
~ThreadPool()
{
    {
        std::lock_guard<std::mutex> lock(mtx);
        stop = true;
    }
    cv_work.notify_all();
    for (auto &w : workers) w.join();
}
```

1. Set `stop = true` under the lock (so workers see it when they next check their predicate).
2. Wake all workers (`notify_all`) so they can observe the stop flag.
3. Wait for all workers to exit their loops via `join()`.

Workers check `if (stop && tasks.empty()) return;` — they finish any remaining tasks before exiting, ensuring no work is lost during shutdown.
