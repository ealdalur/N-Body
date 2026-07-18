#pragma once

#include <vector>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <atomic>

class ThreadPool
{
private:
	std::vector<std::thread> workers;
	std::vector<std::function<void()>> tasks;
	std::mutex mtx;
	std::condition_variable cv_work;
	std::condition_variable cv_done;
	std::atomic<int> active_tasks{0};
	bool stop = false;

public:
	ThreadPool(int numThreads)
	{
		for (int i = 0; i < numThreads; i++) {
			workers.emplace_back([this]() {
				while (true) {
					std::function<void()> task;
					{
						std::unique_lock<std::mutex> lock(mtx);
						cv_work.wait(lock, [this]() { return stop || !tasks.empty(); });
						if (stop && tasks.empty()) return;
						task = std::move(tasks.back());
						tasks.pop_back();
					}
					task();
					if (--active_tasks == 0) {
						cv_done.notify_one();
					}
				}
			});
		}
	}

	~ThreadPool()
	{
		{
			std::unique_lock<std::mutex> lock(mtx);
			stop = true;
		}
		cv_work.notify_all();
		for (auto &w : workers) w.join();
	}

	template<typename F>
	void submit(F&& f)
	{
		{
			std::unique_lock<std::mutex> lock(mtx);
			tasks.emplace_back(std::forward<F>(f));
			active_tasks++;
		}
		cv_work.notify_one();
	}

	void waitAll()
	{
		std::unique_lock<std::mutex> lock(mtx);
		cv_done.wait(lock, [this]() { return active_tasks == 0; });
	}
};
