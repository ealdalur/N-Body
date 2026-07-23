# Dark Energy (FDE Parameter)

## What Is Dark Energy?

### Observational Discovery

In 1998, two independent teams (the Supernova Cosmology Project and the High-z Supernova Search Team) measured the distances and redshifts of Type Ia supernovae and found that the expansion of the universe is *accelerating*. This was unexpected — gravity is attractive, so a universe filled with matter should be decelerating. Something is driving space apart faster and faster.

This unknown repulsive influence was named "dark energy." It accounts for approximately 68% of the total energy content of the observable universe, with dark matter comprising ~27% and ordinary (baryonic) matter only ~5%.

### The Cosmological Constant

The simplest model of dark energy is Einstein's cosmological constant Lambda, introduced as a constant energy density permeating all of space. In general relativity, Lambda enters the Einstein field equations as:

```
G_uv + Lambda * g_uv = (8*pi*G/c^4) * T_uv
```

The cosmological constant produces a repulsive effect because a uniform positive energy density with negative pressure (p = -rho*c^2) generates gravitational repulsion in GR. The key property is that Lambda's energy density remains constant as space expands — unlike matter, which dilutes as volume increases.

### Hubble's Law and Expansion

The expansion of the universe is described by Hubble's law:

```
v = H_0 * d
```

where v is the recession velocity of a distant object, d is its distance, and H_0 ~ 70 km/s/Mpc is the Hubble constant. Every point in space recedes from every other point at a rate proportional to the separation — there is no center.

In an accelerating universe (dark energy dominated), the Hubble "constant" actually increases over time, causing exponential expansion at late times.

### Effect on Structure

Dark energy competes with gravity:
- On small scales (galaxies, galaxy clusters), gravity dominates and structures remain bound.
- On large scales (hundreds of Mpc), dark energy dominates and drives the accelerating expansion.

The transition scale is roughly where the gravitational attraction between structures equals the expansion rate — objects bound more tightly than this threshold are unaffected by dark energy.

---

## Implementation in This Simulator

### The FDE Parameter

The simulation implements dark energy as a single scalar parameter `FDE` (Force Dark Energy). When non-zero, it applies an additional acceleration to every body:

```
acceleration += FDE * position
```

In code (`Simulation.cpp`, applied in both P2P and Octree acceleration loops):
```cpp
vscaleadd(pos_t[i], FDE, acc_t[i]);
```

where `vscaleadd(v, s, vr)` computes `vr += s * v`. This adds `FDE * pos_t[i]` to the acceleration of body i.

### Physical Interpretation

This creates a force that is:
- **Proportional to distance from the origin** — the farther a body is from (0,0,0), the stronger the force
- **Directed radially outward** (when FDE > 0) — pushing bodies away from the origin
- **Uniform in all directions** — isotropic about the origin

This mimics the Hubble flow: every particle is pushed away from the center at a rate proportional to its distance, exactly as Hubble's law dictates (v ~ H*d implies a ~ H_dot*d + H*v, which at the simplest level is proportional to d).

### Setting the Parameter

In a `.sim` script file:
```
FDE    0.001
```

The default is 0.0 (no dark energy). Positive values produce expansion (repulsion), negative values would produce an additional attractive force toward the origin (contraction).

### Dimensional Analysis

The parameter FDE has dimensions of inverse time squared:

```
acceleration [length/time^2] = FDE [1/time^2] * position [length]
```

For a body at rest at distance r from the origin, the outward acceleration is FDE*r. This matches the form of Hubble expansion where the "effective acceleration" in a de Sitter universe is:

```
a = (Lambda*c^2 / 3) * r = H^2 * r
```

So FDE corresponds physically to H^2 (Hubble parameter squared) or equivalently Lambda*c^2/3 in appropriate units.

---

## How This Simplifies Dark Energy

### What This Model Captures

1. **Hubble-like expansion**: Bodies recede from the origin with velocity proportional to distance (once they reach the expansion-dominated regime).

2. **Competition with gravity**: For a cluster of bodies near the origin, gravity pulls them together while FDE pushes them apart. There exists a critical radius beyond which expansion wins — this mimics the real turnaround radius of galaxy clusters.

3. **Acceleration**: The expansion is genuinely accelerating (not just constant velocity). A body that begins at rest at distance r will accelerate outward, gaining speed.

4. **Qualitative structure formation**: With appropriate initial conditions (expanding cloud of particles with gravity), the simulation can show filamentary structure formation with voids — gravity pulls material into clusters and filaments while dark energy expands the voids.

### What This Model Omits

The implementation makes several drastic simplifications relative to true cosmological dark energy:

**1. Preferred center (broken homogeneity)**

Real dark energy (and the Hubble flow) has no center — every point is equivalent. The simulator's implementation has a privileged origin at (0,0,0). A body at the origin feels zero dark energy force, while a body far away feels a strong force. In the real universe, any observer sees the same expansion pattern centered on themselves.

This means the simulation correctly models expansion *as seen from the origin* but introduces spurious asymmetries for structures located far from (0,0,0). A galaxy cluster sitting at (1000, 0, 0) would experience a net outward tidal force that a cluster at the origin would not.

**2. No metric expansion (Newtonian framework)**

Real dark energy operates through the expansion of spacetime itself — it is a property of the metric, not a force. In general relativity, galaxies are not "pushed apart by a force" but rather the space between them grows. Photon wavelengths stretch (cosmological redshift), and the maximum observable distance is finite (cosmic horizon).

The Newtonian force approximation cannot capture:
- Cosmological redshift
- The cosmic event horizon
- Light propagation effects
- The distinction between "recession velocity" and "peculiar velocity"

**3. No scale factor evolution**

In real cosmology, the expansion rate H(t) evolves over time as the universe transitions from radiation-dominated to matter-dominated to dark-energy-dominated eras. The simulator's FDE is a constant — it doesn't evolve with the state of the system.

A proper cosmological simulation would solve the Friedmann equations to evolve H(t) and apply co-moving coordinates with a time-dependent scale factor a(t).

**4. No relativistic effects at large distances**

When recession velocities become comparable to c (which happens at large distances in an expanding universe), special and general relativistic effects become important. The Newtonian approximation breaks down entirely — recession "velocities" exceeding c are physically meaningful in GR but would be nonsensical in the simulator's framework.

**5. No effect on bound systems**

In reality, dark energy has a negligible effect on gravitationally bound systems (the solar system, individual galaxies) because the local gravitational field is enormously stronger than the cosmological expansion. The simulator does not enforce this separation of scales — if FDE is set too large relative to the gravitational forces, it will incorrectly tear apart bound structures.

**6. No self-consistent cosmological initial conditions**

A proper cosmological N-body simulation (like Gadget or AREPO) starts from a nearly uniform density field with small perturbations derived from the CMB power spectrum, and evolves forward using co-moving coordinates with an expanding background metric. This simulator simply applies a radial force to whatever initial conditions are provided.

---

## When This Approximation Is Valid

The FDE force model is a reasonable approximation when:

1. **The region of interest is small compared to c/H** — so relativistic effects are negligible and the Newtonian approximation holds.

2. **Structures are arranged roughly symmetrically around the origin** — minimizing artifacts from the preferred center.

3. **FDE is small compared to gravitational accelerations within bound structures** — so galaxies and clusters remain self-consistently bound while the inter-cluster space expands.

4. **The qualitative behavior (expansion vs. collapse) matters more than quantitative cosmological accuracy** — the simulation can demonstrate the concept of dark energy competing with gravity without reproducing precise Lambda-CDM cosmological evolution.

For pedagogical or visual demonstrations of expansion, void formation, and the large-scale fate of matter in an accelerating universe, this simple model is effective. For quantitative cosmological predictions, it is not.

---

## Example Usage

A simulation of a uniform sphere of particles with dark energy:

```
G          1.0
FDE        0.0005
dt         0.01
r_soft     0.5

Solver     LeapFrog
Gravity    Octree

N_SystemBodies  10000
Camera     0.0  500.0  0.0

SphericalUniverse  0   0.0 0.0 0.0   0.0 0.0 0.0   100.0 200.0 0.0  0.0 0.0
```

With FDE > 0 and a sufficiently sparse initial distribution, the outer shells expand while the denser inner regions collapse into filaments and clusters — a qualitative analog of cosmic structure formation.

---

## Relationship to the Friedmann Equations

For context, the proper cosmological treatment derives the acceleration of a shell of matter at radius r in a homogeneous universe:

```
d²r/dt² = -(4*pi*G/3)*rho*r + (Lambda*c²/3)*r
```

The first term is gravitational deceleration (from matter enclosed within radius r). The second term is the dark energy acceleration. In a matter-free universe (rho=0), the acceleration is purely (Lambda*c²/3)*r — exactly the form `FDE * r` implemented in the simulator.

The FDE parameter thus corresponds to:
```
FDE = Lambda*c² / 3 = H² * Omega_Lambda
```

in the appropriate unit system, where Omega_Lambda ~ 0.68 is the dark energy density parameter.
