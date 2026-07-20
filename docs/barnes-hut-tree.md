# Barnes-Hut Octree Implementation

## The N-Body Problem

Given N bodies with masses and positions, compute the gravitational acceleration on each body due to all other bodies. The direct approach computes every pair:

```
a_i = G * Σ(j≠i) m_j * (r_j - r_i) / |r_j - r_i|³
```

This is O(N²) — for 40,000 bodies that's 1.6 billion interactions per frame. At ~10 FLOPS per interaction, that's 16 GFLOPS per frame, far beyond real-time capability on a CPU.

## The Barnes-Hut Approximation

The key insight: a distant cluster of bodies can be approximated as a single point mass at its center of mass. If a group of 1000 bodies is far enough away, their individual contributions blur together — computing one interaction with their combined mass introduces negligible error compared to computing 1000 separate interactions.

**The opening angle criterion:** A cluster represented by a tree node of width `d` at distance `r` can be approximated if:

```
d / r ≤ θ
```

where θ is the "opening angle" parameter (typically 0.5-1.0). Smaller θ means more accurate but slower. Our implementation uses θ = 1.0 (equivalent to `d² ≤ r²`), which is aggressive but fast.

This reduces the complexity from O(N²) to O(N log N) — each body traverses the tree visiting O(log N) nodes instead of O(N) bodies.

## The Octree Structure

An octree recursively subdivides 3D space into 8 equal octants (the 3D analogue of a quadtree's 4 quadrants). Each node represents a cubic region of space:

```
        +--------+--------+
       /  BNW   /  BNE   /|
      /        /        / |
     +--------+--------+  |
    /  FNW   /  FNE   /|  +
   /        /        / | /|
  +--------+--------+  |/ |
  |        |        |  +  |
  |  FNW   |  FNE   | /|  +
  |        |        |/ | /
  +--------+--------+  |/
  |  FSW   |  FSE   |  +
  |        |        | /
  +--------+--------+/
```

The 8 octants are named by their spatial position:
- **F/B** — Front (z ≥ center) / Back (z < center)
- **N/S** — North (y ≥ center) / South (y < center)
- **E/W** — East (x ≥ center) / West (x < center)

```cpp
enum Octant { FNE=0, FNW, FSW, FSE, BNE, BNW, BSW, BSE };
```

## Node Data Structure

```cpp
struct BHNode
{
    float p_min[3], p_max[3];  // axis-aligned bounding box
    float p_c[3];              // geometric center (midpoint of min/max)
    float p_cm[3];             // center of mass (mass-weighted average position)
    float Mass;                // total mass of all bodies in this node's subtree
    int BodyIdx;               // index of the body (leaf nodes with exactly 1 body)
    int Num_Bodies;            // number of bodies contained in subtree
    int Octants[8];            // child node indices (-1 = no child)
};
```

**Size: 96 bytes.** All spatial data is single-precision float to minimize cache footprint during tree traversal (the hot path). The final acceleration accumulator remains double-precision for numerical stability in summation.

A node is one of three types:
- **Empty** (`Num_Bodies == 0`): no bodies in this region of space
- **Leaf** (`Num_Bodies == 1`): contains exactly one body, stored in `BodyIdx`
- **Internal** (`Num_Bodies > 1`): contains multiple bodies; `p_cm` and `Mass` represent their aggregate

## Memory Management: Arena Allocator

Instead of `new`/`delete` for each node (which would fragment the heap and destroy cache locality), the tree uses a flat `std::vector<BHNode>` as an arena:

```cpp
std::vector<BHNode> nodes;  // pre-allocated pool (512,000 nodes)
int nodeCount;              // next free index

int AllocNode(float *p_min, float *p_max) {
    int idx = nodeCount++;
    if (idx >= nodes.size()) nodes.resize(nodes.size() * 2);
    InitNode(idx, p_min, p_max);
    return idx;
}

void Reset(float *p_min, float *p_max) {
    nodeCount = 0;            // "free" all nodes instantly
    AllocNode(p_min, p_max);  // allocate root
}
```

**Why this is fast:**
- `Reset()` is O(1) — no destructors, no deallocation, just reset the counter
- Nodes are contiguous in memory — sequential allocation means sequential memory layout
- No heap fragmentation — the vector is one contiguous block
- The pre-allocation (512k nodes) avoids resize during normal operation

The tree is completely rebuilt every frame. This sounds expensive, but maintaining a persistent tree (with insertions, deletions, and rebalancing) is more complex and offers similar performance because body positions change every frame anyway.

## Tree Construction

### Phase 1: Determine Bounds

Before building the tree, we compute the axis-aligned bounding cube that contains all bodies:

```cpp
for each body i:
    p_min = componentwise_min(p_min, pos[i])
    p_max = componentwise_max(p_max, pos[i])

// Expand to a cube (equal extent in all axes)
halfwidth = max(p_max - p_min) / 2
center = (p_min + p_max) / 2
p_min = center - halfwidth
p_max = center + halfwidth
```

The bounding region must be a cube (not a rectangular box) because the octree subdivides equally in all three dimensions. A non-cubic root would produce non-cubic children, breaking the octant geometry.

### Phase 2: Insert Bodies

Bodies are inserted one at a time. Each insertion traverses from the root until it finds an empty leaf:

```cpp
void InsertBody(float *p, float m, int BodyIdx)
```

The algorithm uses an explicit stack instead of recursion (avoids function call overhead in a tight loop called 40k times):

**Case 1: Empty node** — Simply place the body here.
```cpp
if (Num_Bodies == 0) {
    p_cm = position;
    Mass = m;
    BodyIdx = body_index;
    Num_Bodies = 1;
}
```

**Case 2: Leaf node (1 body already here)** — Two bodies share the same cell. Subdivide: push the existing body into the appropriate child octant, then continue inserting the new body deeper.
```cpp
if (Num_Bodies == 1) {
    // Determine which octant the existing body belongs to
    oct_existing = DetermineOctant(node, existing_position);
    CreateChild(node, oct_existing);
    // Move existing body into child

    if (existing and new are at the same position) {
        // Merge: combine masses in the child
    } else {
        // Push new body to continue descending
        oct_new = DetermineOctant(node, new_position);
        CreateChild(node, oct_new);
        push(child for new body);
        depth++;
    }
}
```

**Case 3: Internal node (>1 bodies)** — Determine which octant the new body belongs to and descend.
```cpp
if (Num_Bodies > 1) {
    oct = DetermineOctant(node, position);
    CreateChild(node, oct);  // no-op if child exists
    Num_Bodies++;
    push(child);
    depth++;
}
```

### Critical Implementation Detail: Dangling Reference Prevention

`CreateChild` calls `AllocNode`, which may trigger `vector::resize()` if the pool is exhausted. A resize moves the entire array to a new memory location, invalidating all existing pointers and references into it.

The code must never hold a `BHNode &reference` across a `CreateChild` call:

```cpp
// WRONG: n is a reference into the vector
BHNode &n = nodes[idx];
CreateChild(idx, oct);     // may resize → n is now dangling
n.Num_Bodies++;            // undefined behavior!

// CORRECT: copy needed values before allocation, use nodes[idx] after
float existing_cm[3];
vcopyf(nodes[idx].p_cm, existing_cm);   // copy BEFORE potential resize
CreateChild(idx, oct);                    // may resize
nodes[idx].Num_Bodies++;                  // re-index into (possibly moved) vector
```

### Depth Limiting

If two bodies have nearly identical positions (within floating-point resolution), the tree would subdivide indefinitely trying to separate them. A depth limit prevents this:

```cpp
if (depth >= 60) {
    // Just lump the mass into the deepest node reached
    nodes[deepest].Mass += m;
    nodes[deepest].Num_Bodies++;
    break;
}
```

This introduces negligible error (bodies at the same position within 2⁻⁶⁰ of the simulation extent) while preventing stack overflow.

### Octant Determination

Given a node's center point `p_c` and a body position `p`, the octant is determined by three comparisons:

```cpp
Octant DetermineOctant(int idx, float *p) {
    if (p[2] >= p_c[2]) {          // Front vs Back
        if (p[1] >= p_c[1]) {      // North vs South
            if (p[0] >= p_c[0])    // East vs West
                return FNE;
            else
                return FNW;
        } else { ... }
    } else { ... }
}
```

This maps naturally to a 3-bit number: `(z_bit << 2) | (y_bit << 1) | x_bit`, though the implementation uses nested ifs for clarity.

### Child Creation

Creating a child node computes its bounding box as one octant of the parent:

```cpp
int CreateChild(int idx, Octant Oct) {
    // For octant FNE (front-north-east):
    //   min = parent's center
    //   max = parent's max
    // For octant BSW (back-south-west):
    //   min = parent's min
    //   max = parent's center
    // Other octants mix min/center/max per axis

    int child = AllocNode(p_min_new, p_max_new);
    nodes[idx].Octants[Oct] = child;
    return child;
}
```

Each child is exactly half the parent's width in each dimension, maintaining the cubic property.

## Phase 3: Compute Center-of-Mass (CalcMasses)

After all bodies are inserted, internal nodes only know their body count and bounding box — not their aggregate center of mass. `CalcMasses` computes this bottom-up.

Because children always have higher indices than their parents (a consequence of the arena allocator — children are allocated after parents), a simple reverse iteration processes all children before their parents:

```cpp
void CalcMasses() {
    for (int idx = nodeCount - 1; idx >= 0; idx--) {
        if (Num_Bodies <= 1) continue;  // leaves already have their p_cm set

        Mass = 0;
        p_cm = {0, 0, 0};
        for each child:
            Mass += child.Mass;
            p_cm += child.p_cm * child.Mass;
        p_cm /= Mass;  // weighted average
    }
}
```

This replaces a recursive approach — it eliminates function call overhead for potentially 100k+ nodes and produces the same result because the index ordering guarantees children are processed first.

## Force Calculation: Tree Traversal

The core algorithm — computing gravitational acceleration on one body by traversing the tree:

```cpp
void CalcAcceleration(float *p, int BodyIdx, float G, float r_soft, double *a)
```

The traversal uses an explicit stack (512 entries) and visits nodes from root downward:

```cpp
while (stack not empty) {
    pop node

    if (node is a leaf with 1 body) {
        if (not self-interaction) {
            compute direct force from this body
        }
    } else {
        compute distance to node's center of mass

        if (node_width² ≤ distance²) {
            // Far enough — use approximation
            treat entire node as point mass at its center of mass
        } else {
            // Too close — must look at individual children
            push all non-empty children onto stack
        }
    }
}
```

### The Opening Angle Test

```cpp
float d = node.p_max[0] - node.p_min[0];  // node width
float dsq = distance_squared(node.p_cm, body_position);

if (d * d <= dsq) {
    // Use multipole approximation
}
```

This is equivalent to `d/r ≤ 1.0` but avoids computing a square root. Squaring both sides of `d ≤ r` gives `d² ≤ r²`, which uses only quantities we already have.

### Gravitational Force Computation

For a body at position `p` and a node with center of mass `p_cm` and total mass `M`:

```
F = G * M * (p_cm - p) / |p_cm - p|³
```

The denominator `|v|³` where `v = p_cm - p` can be rewritten as `(v·v)^(3/2)` or equivalently `dsq * sqrt(dsq)`. We compute:

```cpp
float dsq = v·v + r_soft;            // softened distance squared
float r3_inv = fast_r3_inv(dsq);     // ≈ 1 / (dsq * sqrt(dsq))
a += v * (G * M * r3_inv);           // accumulate acceleration
```

### Gravitational Softening

At very small separations, gravitational force diverges to infinity (`1/r²` as `r→0`). In a discrete simulation, two bodies passing close can receive enormous accelerations that launch them to infinity.

Softening adds a small constant `ε` to the distance squared:

```
dsq_soft = |r|² + ε
```

This caps the maximum force at `~1/ε` instead of infinity. The parameter `r_soft = 0.1` is chosen to be small relative to typical inter-body distances but large enough to prevent numerical blowup during close encounters.

### Fast Inverse Cube Root

The inner loop requires `1 / (dsq * sqrt(dsq))` — a division and a square root, the two slowest floating-point operations (~12-20 cycles each). We replace them with an SSE hardware approximate reciprocal square root plus one Newton-Raphson refinement:

```cpp
float fast_rsqrtf(float x) {
    // _mm_rsqrt_ss: hardware estimate of 1/sqrt(x), ~12-bit accuracy, ~5 cycles
    // One Newton-Raphson iteration: est = est * (3 - x*est²) / 2
    // Result: ~23-bit accuracy (sufficient for float)
}

float fast_r3_inv(float dsq) {
    float rsqrt = fast_rsqrtf(dsq);   // 1/sqrt(dsq)
    return rsqrt * rsqrt * rsqrt;      // (1/sqrt(dsq))³ = 1/(dsq^(3/2))
}
```

This computes `1/(dsq * sqrt(dsq))` using three multiplies and one hardware rsqrt estimate, versus a `sqrtf` call plus a division in the naive approach.

### Accumulator Precision

The acceleration accumulator `a[3]` is `double` even though the tree uses `float`. This prevents catastrophic cancellation: when summing thousands of small force contributions (each a float), rounding errors accumulate. A double accumulator has 53 bits of mantissa versus float's 23 bits, keeping the sum accurate to ~15 decimal places throughout the traversal.

### Stack Overflow Protection

The traversal stack is sized at 512 entries. In the worst case (all 8 children pushed at each level, depth 60), this could theoretically overflow. A guard prevents this:

```cpp
if (n.Octants[i] != -1 && top < 510) {
    stack[top++] = n.Octants[i];
}
```

In practice, the opening angle criterion ensures most branches are pruned well before reaching deep levels, keeping the stack depth under 100 even for pathological distributions.

## Complete Per-Frame Pipeline

```
1. CalcLeapFrogPositions()     — advance positions by velocity + 0.5*acceleration*dt²
2. PrepareDerivatives()        — set up state array pointers
3. CalcAccelRangeOct()         — parallel tree traversal (each thread: chunk of bodies)
4. CalcLeapFrogVelocities()    — advance velocities by average of old and new acceleration
5. CalcOutputs()               — compute pos², vel², acc² for diagnostics and tree bounds
6. BuildOctree()               — serial: bounds → morton sort → insert 40k bodies → CalcMasses
```

Steps 1-5 are parallelized across all worker threads. Step 6 is serial (tree insertion is inherently sequential), making it the dominant cost at ~55% of frame time.

## Performance Characteristics

With 40,000 bodies on a 16-thread machine:

| Phase | Time | Parallel? |
|-------|------|-----------|
| LeapFrog position | 0.12ms | Yes |
| Acceleration (tree traversal) | 3.8ms | Yes (14 threads) |
| LeapFrog velocity + outputs | 0.15ms | Yes |
| Octree build | 5.0ms | No (serial) |
| **Total** | **~9.1ms** | **~110 FPS** |

The tree traversal benefits from Morton-order body processing (bodies sorted by spatial locality access overlapping tree nodes, improving cache hit rate by ~25%).

## Complexity Analysis

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Tree build (insert all bodies) | O(N log N) | Each insertion is O(depth) ≈ O(log N) |
| CalcMasses | O(M) | M = number of nodes ≈ 2-3N |
| CalcAcceleration (one body) | O(log N) | Opening angle prunes most branches |
| CalcAcceleration (all bodies) | O(N log N) | Parallelized across threads |
| Morton sort (radix) | O(N) | 3 passes over 30-bit keys |
| **Total per frame** | **O(N log N)** | vs O(N²) for direct summation |

For N = 40,000: O(N²) = 1.6 billion operations. O(N log N) ≈ 600,000 operations. A reduction factor of ~2,700×.
