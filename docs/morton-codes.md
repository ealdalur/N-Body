# Morton Codes and Space-Filling Curves for N-Body Simulation

## The Problem: Cache Misses in Tree Traversal

In a Barnes-Hut N-Body simulation, each body traverses the octree to compute its gravitational acceleration. The tree has ~100k+ nodes stored in a flat array. When bodies are processed in arbitrary order, consecutive bodies access completely different regions of the tree, causing constant L1/L2/L3 cache evictions.

Consider two bodies on opposite sides of the simulation space. Their tree traversals share almost no nodes — the first body's traversal loads nodes into cache, then the second body's traversal evicts all of them to load different nodes. This "cache thrashing" means nearly every node access is a cache miss (~100 cycles penalty on modern CPUs).

If instead we process bodies that are **spatially close together** consecutively, their tree traversals will visit many of the same nodes. The first body loads those nodes into cache, and the next body finds them already there — a cache hit (~4 cycles).

The question becomes: how do we sort 3D positions into an order that keeps spatially-nearby points adjacent in the sorted sequence?

## Space-Filling Curves

A space-filling curve is a continuous curve that passes through every point in a bounded region of space. The key property we exploit: points that are near each other on the curve tend to be near each other in space (though the converse isn't always true).

Several space-filling curves exist (Hilbert, Peano, Z-order), but the Z-order curve (Morton curve) is the most practical because computing a point's position on the curve requires only bit manipulation — no lookup tables or recursive computation.

## The Z-Order Curve (Morton Curve)

### 1D Intuition

Consider sorting numbers on a number line. The natural sort order (0, 1, 2, 3, ...) preserves spatial locality perfectly — adjacent numbers are adjacent in memory.

### 2D Intuition

In 2D, we want to map (x, y) coordinates to a single integer such that nearby points get nearby integers. The Morton code achieves this by **interleaving the bits** of the x and y coordinates:

```
x = 5 = 101 (binary)
y = 3 = 011 (binary)

Morton code: interleave y and x bits
  y1 x1 y0 x0 y2 x2
   0  1  1  0  1  1 = 100111 (binary) = 39
```

Wait — let me show this more carefully. For 2D with x=5 (101) and y=3 (011):

```
x bits: 1 0 1
y bits: 0 1 1

Interleaved (x in even positions, y in odd positions):
bit 5: y2 = 0
bit 4: x2 = 1
bit 3: y1 = 1
bit 2: x1 = 0
bit 1: y0 = 1
bit 0: x0 = 1

Result: 010111 = 23
```

### Why This Preserves Locality

The interleaving creates a recursive quadrant structure. The highest two bits encode which quadrant (top-left, top-right, bottom-left, bottom-right) the point falls in. The next two bits encode the sub-quadrant within that quadrant. And so on.

Visually, the Z-curve traces through a 4x4 grid like this:

```
 0  1  4  5
 2  3  6  7
 8  9 12 13
10 11 14 15
```

Notice the Z-shaped pattern at every scale — hence the name. Points within the same quadrant always have adjacent Morton codes, which means sorting by Morton code groups spatially-nearby points together.

### Extension to 3D

For 3D (our N-Body simulation), we interleave three coordinates instead of two. Given x, y, z each as 10-bit integers:

```
x bits: x9 x8 x7 x6 x5 x4 x3 x2 x1 x0
y bits: y9 y8 y7 y6 y5 y4 y3 y2 y1 y0
z bits: z9 z8 z7 z6 z5 z4 z3 z2 z1 z0

Morton code (30 bits):
x9 y9 z9 x8 y8 z8 x7 y7 z7 ... x0 y0 z0
```

Each group of 3 bits selects one of 8 octants (the 8 children of an octree node), which directly corresponds to our BHTree's octant structure. The highest 3 bits select the top-level octant, the next 3 bits select the sub-octant, and so on for 10 levels of subdivision.

## Computing Morton Codes: Bit Expansion

The core operation is "expanding" a 10-bit integer so its bits are spaced 3 apart (with zeros between them), allowing the three coordinates to be OR'd together:

```
Input:        ---- ---- --98 7654 3210   (10 bits packed)
Goal:    9--8 --7- -6-- 5--4 --3- -2-- 1--0   (10 bits spread with 2 gaps each)
```

This is done in stages, each moving groups of bits into their final positions:

```c
uint32_t expandBits(uint32_t v)
{
    v = (v | (v << 16)) & 0x030000FF;  // ---- --98 ---- ---- ---- ---- 7654 3210
    v = (v | (v <<  8)) & 0x0300F00F;  // ---- --98 ---- ---- 7654 ---- ---- 3210
    v = (v | (v <<  4)) & 0x030C30C3;  // ---- --98 ---- 76-- --54 ---- 32-- --10
    v = (v | (v <<  2)) & 0x09249249;  // ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return v;
}
```

Let's trace through with v = 0b1111111111 (all 10 bits set, value 1023):

**Step 1:** `v = (v | (v << 16)) & 0x030000FF`
```
v           = 0000 0000 0000 0000 0000 0011 1111 1111
v << 16     = 0000 0011 1111 1111 0000 0000 0000 0000
v | (v<<16) = 0000 0011 1111 1111 0000 0011 1111 1111
& 030000FF  = 0000 0011 0000 0000 0000 0000 1111 1111
```
Bits 9,8 are now in their final upper position; bits 7-0 remain in the lower byte.

**Step 2:** `v = (v | (v << 8)) & 0x0300F00F`
```
v           = 0000 0011 0000 0000 0000 0000 1111 1111
v << 8      = 0000 0000 0000 0000 1111 1111 0000 0000
v | (v<<8)  = 0000 0011 0000 0000 1111 1111 1111 1111
& 0300F00F  = 0000 0011 0000 0000 1111 0000 0000 1111
```
Now 9,8 are in place; 7,6,5,4 are grouped; 3,2,1,0 are grouped.

**Step 3:** `v = (v | (v << 4)) & 0x030C30C3`
```
Splits each 4-bit group into two 2-bit groups, spaced 4 apart.
Result:       0000 0011 0000 1100 0011 0000 1100 0011
              (98)      (76)      (54)      (32)  (10)
```

**Step 4:** `v = (v | (v << 2)) & 0x09249249`
```
Splits each 2-bit group into individual bits, spaced 3 apart.
Result:       0--9 --8- -7-- 6--5 --4- -3-- 2--1 --0-
              (exact: 1001 0010 0100 1001 0010 0100 1001)
```

The final Morton code is assembled by shifting and OR-ing:
```c
uint32_t morton3D(float x, float y, float z)
{
    uint32_t ix = /* x normalized to 0-1023 */;
    uint32_t iy = /* y normalized to 0-1023 */;
    uint32_t iz = /* z normalized to 0-1023 */;
    return (expandBits(ix) << 2) | (expandBits(iy) << 1) | expandBits(iz);
}
```

The `<< 2`, `<< 1`, `<< 0` offsets place x in the highest bit of each triplet, y in the middle, and z in the lowest.

## Normalization: Mapping Float Positions to 10-bit Integers

Before computing Morton codes, we need to map continuous float positions to the discrete 0-1023 range:

```c
float inv = 1.0f / (p_max[0] - p_min[0]);  // bounds are a cube (equal extent in all axes)
uint32_t ix = (uint32_t)((x - p_min[0]) * inv * 1023.0f);
uint32_t iy = (uint32_t)((y - p_min[1]) * inv * 1023.0f);
uint32_t iz = (uint32_t)((z - p_min[2]) * inv * 1023.0f);
```

We use 10 bits per axis because 3 * 10 = 30 bits fits in a `uint32_t`. This gives 1024 discrete positions per axis, or ~1 billion unique Morton codes. For our simulation with 40k bodies, this is far more resolution than needed.

## Radix Sort: O(N) Sorting by Morton Code

Standard comparison sort (quicksort, merge sort) is O(N log N). For 40k bodies, that's ~640k comparisons. But Morton codes are integers with a fixed bit-width (30 bits), so we can use **radix sort** which is O(N * k) where k is the number of digit passes.

We split the 30-bit Morton code into three 10-bit chunks and sort by each chunk from least-significant to most-significant (LSD radix sort):

```c
for (int shift = 0; shift < 30; shift += 10) {
    // Count occurrences of each 10-bit value
    int counts[1024] = {};
    for (int i = 0; i < numActive; i++)
        counts[(mortonCodes[sortedIdx[i]] >> shift) & 0x3FF]++;

    // Convert counts to prefix sums (output positions)
    int total = 0;
    for (int i = 0; i < 1024; i++) {
        int c = counts[i];
        counts[i] = total;
        total += c;
    }

    // Scatter elements to their sorted positions
    for (int i = 0; i < numActive; i++) {
        int bi = sortedIdx[i];
        sortTemp[counts[(mortonCodes[bi] >> shift) & 0x3FF]++] = bi;
    }

    memcpy(sortedIdx, sortTemp, numActive * sizeof(int));
}
```

**Why 3 passes of 10 bits instead of 30 passes of 1 bit?**

Each pass requires iterating over all N elements three times (count, prefix-sum, scatter). With 1-bit passes we'd need 30 passes = 90N memory accesses. With 10-bit passes we need 3 passes = 9N memory accesses, plus 3 * 1024 for the prefix sums (negligible). The 1024-entry count array fits comfortably in L1 cache.

**Why LSD (least-significant-digit first)?**

LSD radix sort has a key property: it's **stable** — elements with equal keys maintain their relative order from the previous pass. This means after sorting by bits 0-9, then by bits 10-19, then by bits 20-29, the final order is correctly sorted by the full 30-bit key. MSD (most-significant first) requires recursion or extra bookkeeping to maintain sub-bucket ordering.

## Why Morton Order Helps the Barnes-Hut Traversal

When `CalcAcceleration` processes body `i`, it traverses the tree from root to some set of nodes (those close enough to require expansion). If body `i+1` is spatially near body `i`, its traversal visits nearly the same set of nodes. Since body `i` just loaded those nodes into the CPU cache, body `i+1` finds them already resident.

Measured impact in this simulation:
- **Without Morton sort:** Accel phase takes ~5.1ms (40k bodies, 14 threads)
- **With Morton sort:** Accel phase takes ~3.8ms (25% faster)

The improvement scales with N — larger simulations have deeper trees with more nodes, making cache misses more expensive and Morton ordering more beneficial.

## Relationship Between Morton Codes and the Octree Structure

The 3D Morton code has a direct correspondence to the octree:

- Bits 29-27 (highest triplet) = which octant of the root node the body falls in
- Bits 26-24 = which sub-octant of that octant
- Bits 23-21 = next level down
- ... and so on for 10 levels

This means Morton-sorted bodies are naturally in a depth-first octree traversal order. Bodies in the same leaf node have identical Morton code prefixes. Bodies in the same branch of the tree have similar Morton codes.

This is also why Morton-order insertion produces better tree construction cache behavior — consecutive insertions navigate the same branch of the tree, keeping those nodes in cache.

## Limitations

1. **Not perfectly locality-preserving.** The Z-curve has "jumps" where spatially-adjacent points get distant Morton codes (at quadrant/octant boundaries). The Hilbert curve is better at preserving locality, but its computation is significantly more complex (lookup tables, state machines) and the practical difference is small.

2. **Resolution limited by bit width.** With 10 bits per axis, two bodies must be more than 1/1024th of the simulation extent apart to get different Morton codes. Bodies closer than this get the same code and their relative ordering is arbitrary (but stable within a radix sort pass).

3. **Sorting cost.** The radix sort itself takes ~0.5ms for 40k bodies — this is "free" only if the traversal savings exceed it. For very small N (<1000), the sort overhead may exceed the cache benefit.
