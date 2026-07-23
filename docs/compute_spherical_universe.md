# compute_spherical_universe.py

## 1. Cosmological Background: Universe Types and Expansion

### 1.1 The Expanding Universe

The universe is expanding: every point in space recedes from every other point at a rate proportional to the separation between them. This is described by Hubble's law:

```
v = H * r
```

where v is the recession velocity, r is the distance, and H is the Hubble parameter. This relationship — velocity proportional to distance — is called **Hubble flow**. It is not motion through space but the expansion of space itself, though in a Newtonian simulation we model it as literal particle velocities.

### 1.2 The Critical Density and H_crit

For a uniform sphere of matter expanding with Hubble flow, there is a critical expansion rate at which the total kinetic energy exactly equals the gravitational binding energy. This defines the **critical Hubble parameter**:

```
H_crit = sqrt(2*G*M / R^3)
```

Derivation:
- Total kinetic energy of Hubble flow in a uniform sphere:
  ```
  KE = integral_0^R (1/2)*rho*(H*r)^2 * 4*pi*r^2 dr = (3/10)*M*H^2*R^2
  ```
- Gravitational self-energy (binding energy) of a uniform sphere:
  ```
  PE = -(3/5)*G*M^2 / R
  ```
- Setting KE = |PE|:
  ```
  (3/10)*M*H^2*R^2 = (3/5)*M^2/R
  H_crit^2 = 2*M/R^3        [G=1 in code units]
  ```

### 1.3 Total Energy and Universe Classification

The total energy E = KE + PE determines the fate of the universe:

```
E = (3/10)*M*H^2*R^2 - (3/5)*M^2/R
  = |PE| * [(H/H_crit)^2 - 1]
```

The ratio KE/|PE| = (H/H_crit)^2 determines the classification:

| H vs H_crit | KE/|PE| | Total E | Fate | Classification |
|-------------|---------|---------|------|----------------|
| H > H_crit | > 1 | E > 0 | Expands forever | **Open** |
| H = H_crit | = 1 | E = 0 | Barely expands forever | **Flat** |
| H < H_crit | < 1 | E < 0 | Expands then recollapses | **Closed** |

### 1.4 The Density Parameter Omega

In cosmology, the density parameter Omega is defined as:

```
Omega = rho / rho_crit
```

where rho_crit is the density that gives H = H_crit for the observed expansion rate. Equivalently, Omega = |PE|/KE for the observed Hubble flow:
- Omega > 1: closed universe (too much matter for the expansion rate)
- Omega = 1: flat universe
- Omega < 1: open universe (not enough matter to halt expansion)

In our simulation, with H/H_crit = 0.8, we have KE/|PE| = 0.64, which corresponds to Omega = 1/0.64 = 1.56 — a closed universe with 56% more mass than needed for flatness.

### 1.5 Dark Energy and the Cosmological Constant

In the real universe, a third component — dark energy — acts as a repulsive force proportional to distance. In our simulation, this is the FDE parameter:

```
a_dark_energy = FDE * r    (outward, proportional to distance)
```

Dark energy competes with gravity:
```
a_gravity = G*M_enclosed / r^2 = M/r^2    [G=1, for matter outside the sphere]
```

The **FDE balance** is the value at which dark energy exactly cancels gravity for the mean density at the initial radius:

```
FDE_balance = M/R^3 = (4/3)*pi*G*rho    [G=1]
```

- FDE = 0: pure gravity (closed universe recollapses)
- FDE < FDE_balance: gravity still dominates at the initial scale
- FDE = FDE_balance: static (Einstein static universe analog, unstable)
- FDE > FDE_balance: dark energy dominates, expansion accelerates

### 1.6 The Dominance Ratio

At any radius r (with all mass M enclosed inside), the ratio of dark energy acceleration to gravitational acceleration is:

```
ratio = a_FDE / a_grav = (FDE * r) / (M / r^2) = FDE * r^3 / M
```

This ratio grows as r^3 — dark energy becomes increasingly dominant as the universe expands. The crossover radius where ratio = 1 (dark energy equals gravity) is:

```
r_crossover = (M / FDE)^(1/3)
```

For our simulation (FDE=0.10, M=5e7):
```
r_crossover = (5e7 / 0.10)^(1/3) = 794 code units
```

At the initial radius (r=200):  ratio = 0.10 * 200^3 / 5e7 = 0.016 = **1.6%**
At r=100:  ratio = 0.10 * 100^3 / 5e7 = 0.002 = **0.2%**

This means dark energy is completely negligible during structure formation (which occurs at r ~ 100-300) but becomes significant once the universe has expanded to r ~ 800, delaying the recollapse.

---

## 2. Gravitational Instability and Structure Formation

### 2.1 How Filaments Form

In a uniformly expanding medium, if one region has slightly higher density than average, that region experiences stronger self-gravity. It decelerates faster than the mean expansion, eventually halting its local expansion and collapsing. Meanwhile, underdense regions decelerate less and become voids.

The collapse is anisotropic: matter doesn't collapse equally in all directions (that would form spheres). It first collapses along the shortest dimension, forming sheets (pancakes/Zel'dovich pancakes), then filaments (collapse along two axes), and finally nodes (three-axis collapse). The result is the "cosmic web" — a network of filaments connecting dense nodes, with voids between them.

### 2.2 The Jeans Length

Not all overdensities collapse. Thermal motions (random velocities) provide pressure support that resists collapse below a critical scale. The **Jeans length** is this minimum scale:

```
lambda_J = v_thermal * sqrt(pi / (G * rho))    [G=1]
```

- Perturbations larger than lambda_J: gravity wins, collapse occurs
- Perturbations smaller than lambda_J: thermal pressure prevents collapse

In our Hubble-flow simulation, the effective "thermal" velocity arises from the Gaussian position scatter in the particle generator. A particle displaced by delta_r from its ideal position has a velocity difference of H*delta_r relative to its neighbors. With sigma=10:

```
v_thermal = H * sigma = 2.83 * 10 = 28.3 km/s
```

This gives lambda_J = 41 code units. Structure forms on scales larger than this; smaller-scale perturbations are washed out by thermal motion.

### 2.3 The Jeans Collapse Time

Once an overdensity exceeds the Jeans scale, it collapses on the free-fall timescale:

```
t_Jeans ~ 1 / sqrt(G * rho) = 1 / sqrt(rho)    [G=1]
```

For our parameters: t_Jeans = 1/sqrt(1.49) = 0.82 code units.

This is the characteristic time for filaments to form. Structure becomes visible around 1-3 Jeans times (t ~ 0.8-2.5 code units).

### 2.4 Thermal Pressure and Small-Scale Collapse

**Thermal pressure** in this context refers to the random kinetic energy of particles that acts as an isotropic pressure resisting gravitational collapse — analogous to gas pressure in a real fluid.

In our simulation, the thermal motion arises from the Gaussian position scatter (sigma=10). Displaced particles have different Hubble velocities than their neighbors, creating a velocity dispersion of H*sigma ~ 28 km/s.

If the thermal velocity is too large relative to the virial velocity, pressure support overwhelms gravity at all scales and no structure forms. In an earlier failed attempt, a random isotropic velocity of 200 km/s was used instead of Hubble flow — this was ~28% of the escape velocity and completely destroyed structure by washing out the Gaussian-seeded density perturbations via thermal mixing.

Our choice (v_thermal ~ 28 km/s vs v_virial ~ 500 km/s) ensures that thermal pressure only prevents collapse below ~41 code units, while larger-scale filaments form freely.

### 2.5 Relationship to Escape Velocity

The escape velocity from the sphere is:
```
v_escape = sqrt(2*G*M/R) = sqrt(2*M/R) = 707 km/s    [G=1]
```

Note that H_crit * R = v_escape — this is not a coincidence but a direct consequence of the critical condition KE = |PE|. The edge velocity at H_crit equals the escape velocity.

In our Hubble flow, the edge velocity is H*R = 566 km/s (80% of v_escape). This is coherent radial outflow, not random thermal motion. The actual thermal dispersion is only 28 km/s (4% of v_escape), which is why structure forms — the random component is tiny relative to the gravitational binding energy.

---

## 3. Inter-Particle Spacing and Density Perturbation Seeding

### 3.1 Inter-Particle Spacing

For N particles uniformly distributed in a sphere of radius R, the mean distance between nearest neighbors is approximately:

```
spacing = R * (V_sphere / N)^(1/3) = R * ((4/3)*pi / N)^(1/3)
```

For N=100,000, R=200: spacing = 6.9 code units.

### 3.2 Gaussian Position Scatter

The particle generator adds a random Gaussian displacement (sigma=10 code units in each axis) to each particle's ideal position. This creates local density fluctuations that seed structure formation.

The ratio **scatter/spacing** determines whether the perturbations are coherent:
- scatter/spacing < 1: particles barely overlap their neighbors' territory; perturbations are weak and noisy
- scatter/spacing ~ 1: particles intermingle, creating moderate coherent overdensities
- scatter/spacing > 1: strong mixing, well-defined density clumps

For N=100,000: scatter/spacing = 10/6.9 = 1.44. This is in the good regime — each particle is displaced by more than one inter-particle spacing, creating coherent clumps of ~100 particles that serve as seeds for filamentary collapse.

### 3.3 Why Particle Count Matters

| N | Spacing | Scatter/Spacing | Particles per Jeans Volume | Quality |
|---|---------|-----------------|---------------------------|---------|
| 20,000 | 11.9 | 0.84 | 22 | Poor — noisy, weak perturbations |
| 40,000 | 9.4 | 1.06 | 43 | Marginal — visible but grainy |
| 100,000 | 6.9 | 1.44 | 108 | Good — clean filaments |
| 200,000 | 5.5 | 1.81 | 216 | Excellent — sharp, well-resolved |

More particles provide: better sampling of the density field, more particles per collapsing structure (smoother filaments), and stronger coherent perturbation seeding from the Gaussian scatter.

---

## 4. Python Script Implementation

### 4.1 Structure

The script is organized into:
1. Parameter definitions (R, N, M_total)
2. Critical Hubble parameter derivation
3. Hubble parameter selection (0.8 * H_crit)
4. Dark energy analysis (FDE dominance ratio at various radii)
5. Jeans analysis (thermal velocity, Jeans length, collapse time)
6. Timescale comparison (free-fall, Jeans, turnaround)
7. Particle count resolution analysis
8. Final parameter output

### 4.2 Key Derivations in the Script

**Critical Hubble parameter:**
```python
H_crit = math.sqrt(2.0 * M_total / R**3)
```
From KE = |PE| for a uniform sphere with Hubble flow.

**Energy classification:**
```python
KE_frac = H_frac**2  # = (H/H_crit)^2 = KE/|PE|
```
Since KE = (3/10)*M*H^2*R^2 and |PE| = (3/5)*M^2/R, the ratio is H^2*R^3/(2M) = (H/H_crit)^2.

**Dominance ratio:**
```python
ratio_at_200 = FDE_use * 200**3 / M_total
```
The fraction of gravitational acceleration that dark energy represents at radius 200.

**Jeans length:**
```python
v_thermal = H_use * 10.0  # H * sigma_position
lambda_J = v_thermal * math.sqrt(math.pi / rho)
```
Thermal velocity from position scatter, then standard Jeans formula.

**Turnaround radius (without FDE):**
```python
r_max_no_fde = R / (1.0 - KE_frac)
```
From energy conservation: E = (1/2)*v^2 - M/r, at turnaround v=0 so r_max = M/|E|. For the outer shell: E per unit mass = (1/2)*(H*R)^2 - M/R = (M/R)*(KE_frac - 1), giving r_max = R/(1-KE_frac).

### 4.3 How Parameters Were Determined

The parameter selection process was iterative:

1. **R=200, M=5e7 chosen** to make inter-particle spacing (9.4 for N=40k, 6.9 for N=100k) comparable to or smaller than the Gaussian scatter sigma=10. This ensures coherent overdensity seeding.

2. **H = 0.8 * H_crit chosen** to make the system bound (closed) but with enough expansion energy that filaments form before recollapse. H too close to H_crit (e.g., 0.95) gives very slow expansion where everything collapses into a blob. H too low (e.g., 0.5) gives rapid collapse before filaments can differentiate.

3. **FDE=0.10 chosen empirically** after testing several values:
   - FDE=0.44: too strong, dominated expansion during structure formation (ratio=7% at r=200), prevented filament formation
   - FDE=0.25: still too aggressive, expansion outpaced local collapse
   - FDE=0.10: negligible during formation (1.6% at r=200), only relevant at r>794, successfully delays recollapse while preserving filamentary structure

4. **N=100,000 chosen** to give scatter/spacing=1.44 (well above 1.0, ensuring strong coherent perturbation seeding) and 108 particles per Jeans volume (smooth filaments).

---

## 5. Expected Simulation Behavior

### 5.1 Timeline

| Phase | Sim time | Physical time | Description |
|-------|----------|---------------|-------------|
| Initial expansion | t = 0-0.3 | 0-18 Myr | Sphere expands uniformly. Density perturbations begin growing but are not yet visible. |
| Filament emergence | t = 0.3-1.0 | 18-59 Myr | Linear overdensities become visible as streaks. Voids begin opening between them. |
| Cosmic web peak | t = 1.0-3.0 | 59-176 Myr | Fully developed web with nodes (dense clusters at filament intersections), filaments (elongated collapsed structures), and voids (expanding empty regions). |
| Late expansion | t > 3 | > 176 Myr | Bulk expansion slows. FDE provides gentle outward support. Cosmic web structure persists. Eventually (t >> 10) the closed nature means recollapse, but on a very long timescale. |

Note: The absolute physical timescales derive from the code unit system (1 time unit = 58.7 Myr) but are not cosmologically realistic. The real universe's structure formation spans ~10 Gyr. The simulation captures the correct physical mechanism and qualitative progression, not the quantitative timeline.

### 5.2 What to Watch For

- **Filaments**: Linear chains of particles connecting denser nodes. These are anisotropic collapses (two-axis collapse, one axis still expanding).
- **Nodes**: Dense clusters at the intersections of 3+ filaments. These are three-axis collapses.
- **Voids**: Large empty regions between filaments that continue expanding.
- **Hierarchy**: Small structures form first (bottom-up), merging into larger filaments and clusters over time.
- **Velocity field**: Particles in filaments stream inward along the filament axis toward the nearest node.

### 5.3 Comparison to Previous (Failed) Approach

| Parameter | Old (evaporated) | New (works) |
|-----------|-----------------|-------------|
| Velocity type | Random directions | Hubble flow (v = H*r) |
| FDE | 0.25 (dominant) | 0.10 (negligible early) |
| What drives expansion | FDE | Initial Hubble flow |
| What seeds structure | Random velocity mixing | Gaussian position scatter |
| Why it works | N/A (failed) | Differential deceleration of over/underdense regions |
| Failure mode | Center collapses, exterior evaporates | N/A (correct structure) |

The key insight: structure formation requires **coherent expansion** (Hubble flow) with **differential deceleration** (gravity stronger in overdense regions). Random velocities + FDE cannot achieve this because FDE drives expansion from a preferred center, not uniformly.
