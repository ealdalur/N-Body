# LeapFrog Integration

## Overview

LeapFrog integration (also known as the Störmer-Verlet or Velocity Verlet method) is a second-order symplectic integrator used to advance Newton's equations of motion in time. It is the standard choice for gravitational N-body simulation because it conserves energy over arbitrarily long timescales — a property that non-symplectic methods (Euler, RK4) fundamentally lack.

## The Equations of Motion

An N-body gravitational system evolves under:

```
dx/dt = v
dv/dt = a(x) = -∇Φ(x)
```

where `Φ(x)` is the gravitational potential. The acceleration `a` depends only on positions, not velocities — this is what makes the system Hamiltonian and allows symplectic integration.

## Derivation

### Starting from Taylor expansion

Expand position forward and backward in time around `t`:

```
x(t + dt) = x(t) + v(t)·dt + ½·a(t)·dt² + ⅙·j(t)·dt³ + ...
x(t - dt) = x(t) - v(t)·dt + ½·a(t)·dt² - ⅙·j(t)·dt³ + ...
```

Adding these two expansions cancels all odd-order terms:

```
x(t + dt) = 2·x(t) - x(t - dt) + a(t)·dt²  +  O(dt⁴)
```

This is the **Störmer form** (position-only, second-order accurate with fourth-order local error). It requires no velocity at all — positions "leapfrog" over each other. However, we need velocities for physical diagnostics (kinetic energy, momentum) and for the initial kick.

### Velocity Verlet form

The Velocity Verlet reformulation makes velocities explicit. Start from:

```
x(t + dt) = x(t) + v(t)·dt + ½·a(t)·dt²
```

This is just the Taylor expansion truncated at the acceleration term. Then compute the new acceleration:

```
a(t + dt) = f(x(t + dt))
```

And update velocity using the trapezoidal rule (average of old and new acceleration):

```
v(t + dt) = v(t) + ½·[a(t) + a(t + dt)]·dt
```

This is mathematically equivalent to the Störmer form but carries velocities synchronized with positions. The algorithm per step is:

```
1. x_new = x + v·dt + ½·a·dt²
2. a_new = f(x_new)
3. v_new = v + ½·(a + a_new)·dt
```

### Why "LeapFrog"?

In the original leapfrog formulation, velocities are stored at half-integer time steps (`t + dt/2`) while positions are at integer steps (`t`). They alternate, each "leaping" over the other:

```
v(t + dt/2) = v(t - dt/2) + a(t)·dt
x(t + dt)   = x(t) + v(t + dt/2)·dt
```

The Velocity Verlet form we use is algebraically identical but keeps positions and velocities synchronized at the same time step, which is more convenient for output and initialization.

## Symplectic Integration

### What "symplectic" means

A Hamiltonian system has a conserved quantity — the Hamiltonian `H = T + V` (kinetic plus potential energy). The exact flow of the system preserves the symplectic 2-form: it conserves phase-space volume and the geometric structure of the Hamiltonian.

A **symplectic integrator** is a discrete map that exactly preserves this structure. It does not conserve the *true* Hamiltonian exactly, but it conserves a *nearby* Hamiltonian (called the shadow Hamiltonian) `H̃ = H + O(dt²)` exactly. This means:

- Energy oscillates around the true value but never drifts
- The bounded energy error is proportional to `dt²` and remains bounded for all time
- Phase-space volume is exactly preserved (no artificial dissipation or growth)

### Why this matters for N-body

In a galaxy simulation running for billions of dynamical times, energy conservation is the critical requirement. A method that drifts even slightly will:

- Cause orbits to spiral inward (losing energy) or outward (gaining energy)
- Destroy the virial equilibrium of the system
- Produce unphysical heating or cooling of the stellar "fluid"

LeapFrog's symplectic property guarantees none of this happens regardless of how long the simulation runs.

## Comparison with Other Methods

### Forward Euler

```
x_new = x + v·dt
v_new = v + a(x)·dt
```

- First-order accurate: error ~ O(dt)
- **Not symplectic**: energy drifts monotonically. Orbits spiral outward (Euler adds energy)
- Requires impractically small dt to control drift
- Useless for long-duration gravitational simulations

### Classical RK4 (4th-order Runge-Kutta)

```
k1 = f(y)
k2 = f(y + dt/2 · k1)
k3 = f(y + dt/2 · k2)
k4 = f(y + dt · k3)
y_new = y + (dt/6)·(k1 + 2·k2 + 2·k3 + k4)
```

- Fourth-order accurate: local error ~ O(dt⁵)
- **Not symplectic**: energy drifts, typically losing energy (orbits spiral inward)
- Requires 4 force evaluations per step (vs 1 for LeapFrog)
- Higher instantaneous accuracy per step, but the accuracy advantage is irrelevant when the simulation must run for thousands of orbital periods

#### RK4's hidden cost

RK4's higher order means it can take larger time steps for the same *per-step* accuracy. But for N-body:

- The force evaluation dominates runtime (O(N log N) for a tree code). Four evaluations per step means RK4 costs 4× per step.
- To match LeapFrog's long-term energy behavior, RK4 would need such small steps that the 4× cost makes it far more expensive overall.
- RK4 also requires storing intermediate state estimates (additional memory proportional to N), while LeapFrog only needs the current acceleration and the previous acceleration.

#### When RK4 is appropriate

RK4 is the right choice for dissipative systems (systems with friction, drag, or damping) where energy is not conserved anyway, or for short-duration high-accuracy problems (spacecraft trajectory optimization over a single orbit). For conservative gravitational systems evolved over many dynamical times, it is strictly inferior to LeapFrog.

### Summary table

| Property | Euler | RK4 | LeapFrog |
|----------|-------|-----|----------|
| Order of accuracy | 1 | 4 | 2 |
| Symplectic | No | No | Yes |
| Energy behavior (long run) | Diverges (grows) | Drifts (shrinks) | Bounded oscillation |
| Force evaluations per step | 1 | 4 | 1 |
| Time-reversible | No | No | Yes |
| Memory per body | 6 doubles | 24+ doubles | 12 doubles |

## Relationship to Velocity Verlet

LeapFrog and Velocity Verlet are the **same algorithm** expressed differently:

- **LeapFrog** (kick-drift-kick): velocity is half-stepped, then position is full-stepped, then velocity is half-stepped again. Velocities naturally live at half-integer times.
- **Velocity Verlet** (drift-kick): position is full-stepped using current velocity and acceleration, then new acceleration is computed, then velocity is full-stepped using the average of old and new acceleration.

They produce identical trajectories to machine precision. The Velocity Verlet form is used in this codebase because it keeps positions and velocities synchronized, making output and initialization straightforward.

## Time Reversibility

LeapFrog is time-reversible: if you negate all velocities and run the same number of steps, you return exactly to the starting configuration (to floating-point precision). This is a consequence of the symplectic property and the symmetric structure of the update equations.

This provides a useful verification test: run forward N steps, negate velocities, run forward N more steps, and check that you recover the initial positions.

## Implementation in This Codebase

### Data layout

Each body `i` has:
- `pos[i]` → pointer to 3 doubles (x, y, z position) in `pos_data`
- `vel[i]` → pointer to 3 doubles (vx, vy, vz velocity) in `vel_data`
- `acc[i]` → pointer to 3 doubles (ax, ay, az current acceleration) in `acc_data`
- `acc_prev[i]` → pointer to 3 doubles (previous step's acceleration) in `acc_prev_data`

### Step function

```cpp
void Simulation::Step()
{
    CalcLeapFrogPositions();           // also drifts the inertial halo centres
    CalcAccelerations();               // includes IntegrateHaloCenters()
    CalcLeapFrogVelocitiesAndOutputs(); // also kicks the inertial halo centres

    BuildOctree();  // rebuild spatial index for next step
    t += dt;
}
```

The halo centres are integrated with the same velocity-Verlet step as the particles: `CalcLeapFrogPositions()` drifts them, `IntegrateHaloCenters()` (called inside `CalcAccelerations()`) sets their acceleration from the disc back-reaction and the other halos, and `CalcLeapFrogVelocitiesAndOutputs()` applies the velocity kick. See `docs/dark-matter-halo.md`.

### Phase 1: Advance positions (`CalcLeapFrogPositionsRange`)

```cpp
for (int i = iStart; i <= iEnd; i++)
{
    vcopy(acc[i], acc_prev[i]);           // save current a for velocity update

    vscaleadd(vel[i], dt, pos[i]);        // pos += vel * dt
    vscaleadd(acc[i], 0.5*dt*dt, pos[i]); // pos += 0.5 * a * dt²
}
```

This implements: **x(t+dt) = x(t) + v(t)·dt + ½·a(t)·dt²**

The acceleration is saved to `acc_prev` before being overwritten in Phase 3.

### Phase 2: Compute new accelerations (`CalcAccelerations`)

Evaluates `a(t+dt) = f(x(t+dt))` using either the Barnes-Hut octree (O(N log N)) or direct particle-to-particle summation (O(N²)), plus the analytic dark matter halo forces.

### Phase 3: Advance velocities (`CalcLeapFrogVelocitiesAndOutputsRange`)

```cpp
for (int i = iStart; i <= iEnd; i++)
{
    vadd(acc[i], acc_prev[i], a);         // a = a_new + a_old
    vscaleadd(a, 0.5*dt, vel[i]);         // vel += 0.5 * (a_new + a_old) * dt

    pos_sq[i] = vdot(pos[i], pos[i]);     // diagnostics (fused for cache efficiency)
    vel_sq[i] = vdot(vel[i], vel[i]);
    acc_sq[i] = vdot(acc[i], acc[i]);
}
```

This implements: **v(t+dt) = v(t) + ½·[a(t) + a(t+dt)]·dt**

The velocity and diagnostic computations are fused into a single loop to avoid a second pass over all body data.

### Parallelization

All three phases are parallelized across the thread pool. Each worker processes a contiguous chunk of bodies. The phases are separated by barriers (`pool->waitAll()`):

1. Position update (all threads) → barrier
2. Pin central bodies + compute accelerations (all threads for tree traversal) → barrier
3. Velocity update + diagnostics (all threads)

The barriers ensure that all positions are updated before forces are computed, and all forces are computed before velocities are updated.

### Initialization

At construction time, positions and velocities are set by the initial conditions (galaxy disc, solar system, etc.). The first call to `CalcAccelerations()` populates `acc[i]` from the initial positions. On the first `Step()`, `acc_prev[i]` is set from this initial acceleration, so the velocity update on the first step is exact (no "cold start" error).

### Cost per step

| Phase | Work | Parallelized |
|-------|------|-------------|
| Position update | O(N) | Yes |
| Force computation | O(N log N) or O(N²) | Yes |
| Velocity update + diagnostics | O(N) | Yes |
| Octree rebuild | O(N log N) | No (serial insertion) |

The force computation dominates. With a single force evaluation per step, LeapFrog extracts maximum science per FLOP.
