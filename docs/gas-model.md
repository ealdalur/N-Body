# Gas Modelling Methodology

This document describes the dissipative gas model used in the simulator — the
"sticky particle" scheme with a kinetic Monte-Carlo (DSMC-style) collision step,
the spatial grid used for collision detection, and the per-species gravitational
softening that keeps the gas from fragmenting. It records the algorithm precisely
as implemented (`Simulation::ProcessGasCollisions` and `Simulation::SoftSq`) and
cites the literature the method descends from.

All quantities are in the code unit system (see `docs/compute_M51.md`): 1 length
unit = 60 pc, 1 velocity unit = 1 km/s, 1 time unit = 58.7 Myr, and `G = 1`.
Because 1 velocity unit is 1 km/s, dispersions quoted below in "km/s" are also the
raw code values.

---

## 1. Motivation and overview

A collisionless N-body disc reproduces the *stellar* response to a tidal
encounter, but the interstellar gas behaves differently: it is dissipative. Gas
clouds collide, radiate away the kinetic energy of their relative motion, and
settle into thin, cool, high-contrast features — the narrow dust/gas lanes and
sharpened spiral arms and bridges that dominate observed interacting-galaxy
morphology. To capture this without solving the full hydrodynamic equations we use
the **sticky particle** approximation: a subset of the disc particles is tagged as
gas, gravitates exactly like the stars, and *additionally* undergoes inelastic
collisions that remove energy from random motion.

This mirrors Salo & Laurikainen (2000), whose M51 model this project reproduces.
Following that paper:

- Each disc's baryonic mass is split into collisionless **stars** and dissipative
  **gas** at the ratio `Nstar : Ngas = 4 : 1` (20 % gas by count and, at equal
  particle mass, by mass).
- The gas is initialised **warm** — the same Toomre `Q = 1.5` Gaussian velocity
  dispersion as the stars (Toomre 1964) — not cold. The cool `σ_gas ≈ 5–10 km/s`
  equilibrium the paper reports is an *emergent* result of the collisional cooling
  during the isolation phase, reached where collisions are frequent (the dense
  inner disc) while low-density outer gas stays warm.

Gas particles are the **last** `round(gas_fraction · N)` bodies of each system's
index range (see the `GalaxyDisc` command in `docs/script-format.md`); a `char`
array `is_gas[i]` tags them.

---

## 2. Where the method comes from

The sticky-particle idea entered galactic dynamics in the early 1980s (e.g.
Combes & Gerin 1985) as a cheap stand-in for gas dissipation: represent the ISM by
particles with a finite collision cross-section and bleed off the normal component
of their relative velocity on impact.

The specific scheme we follow is Salo's. Salo & Laurikainen (2000) state only the
per-impact rule — "the component of the relative velocity in the direction joining
the particle centres is reversed and reduced by a factor α" — and refer to
**Salo (1991)** and **Salo & Laurikainen (1993, "SL93")** for the full method,
including how collision partners are found and how often collisions are resolved.
Those details (a grid-based, rate-controlled collision search) are exactly the
structure of the **Direct Simulation Monte Carlo (DSMC)** method that Bird
developed for rarefied gas dynamics (Bird 1963, 1994). Our implementation is an
explicit DSMC step — Bird's **No-Time-Counter (NTC)** pair selection — grafted
onto the leapfrog integrator. This is a deliberate, faithful reconstruction of the
paper's collision physics rather than a verbatim reimplementation of its code
(which we do not have).

Why DSMC and not a literal geometric-overlap test is explained in §5.

---

## 3. Operator splitting

Collisions are handled by an **operator-split** step. Within one leapfrog step the
sequence is:

1. gravity drift/kick (positions and velocities advanced under the smoothed
   potential — stars and gas together, in the Barnes–Hut octree);
2. **`ProcessGasCollisions()`** — the DSMC collision step, which modifies gas
   *velocities only*, leaving positions untouched.

Splitting is legitimate here because the collision time-scale and the orbital
time-scale are well separated, and it keeps the collision operator independent of
the gravity solver. The collision step runs on **every** step, including the
pre-`t = 0` isolation ("warmup") phase, so an isolated gas disc cools to its
collisional equilibrium before the encounter begins — as in the paper's setup.
Star–star and star–gas pairs never collide; only gas–gas pairs do.

---

## 4. The grid / collision-detection mechanism

DSMC decouples *collision partnering* from exact particle positions. Instead of
asking "which pairs physically overlap this instant?", it partitions space into
**collision cells** and, within each cell, treats every gas particle as a possible
partner for every other, selecting actual collisions stochastically at the correct
physical rate.

### 4.1 The spatial hash grid

Each step we bin the gas into a uniform cubic grid of edge `L = Gas_Cell_Size`
(code units). The grid is a hash map keyed on the integer cell coordinates, built
over gas particles only, so it costs `O(N_gas)` and holds only occupied cells:

```
cell(p) = ( floor(p.x / L), floor(p.y / L), floor(p.z / L) )
```

The three signed cell indices are packed into a single 64-bit key
(21 bits each, offset by 2^20 to keep them non-negative) and used as the hash-map
key. All gas particles sharing a key are collision candidates for each other.

**Cell-size choice matters, but not for the rate.** `Gas_Cell_Size` must be large
enough to contain several *in-plane* neighbours (so in-plane collisions can occur)
yet small compared with the disc/arm scale and the scale of velocity gradients
(so a cell is roughly uniform). Crucially, the collision **rate is independent of
`L`**: the per-particle rate scales with the local *number density* `n = Nc/Vc`
(§4.3), and both the cell population `Nc` and the estimated volume `Vc` grow
together as `L` grows. `L` only sets the spatial locality of collision partners,
not how often collisions happen.

### 4.2 Per-cell kinematics

For each cell with `Nc ≥ 2` gas particles we compute, in a single pass:

- the mean cell velocity `v̄`;
- an **upper bound on the pairwise relative speed**,
  `v_rel,max = 2 · max_i |v_i − v̄|` (since `|v_i − v_j| ≤ |v_i − v̄| + |v_j − v̄|`);
- the **bounding box** of the cell's particles, `[lo, hi]` per axis.

The collision volume is the bounding box, floored on each axis at the collision
diameter `2·Gas_Radius`:

```
Vc = Π_k  max( hi[k] − lo[k],  2·Gas_Radius )
```

Using the occupied bounding box rather than the full `L³` is deliberate: a galactic
gas disc is a thin pancake, so its particles fill only a slab within a cubic cell.
Taking `Vc = L³` would badly *under*-estimate the local density and suppress the
collision rate; the bounding box tracks the true occupied volume (and adapts
automatically to a tilted disc). This is an approximation — see §8.

### 4.3 NTC candidate selection

The expected number of collisions in a cell over a step `Δt` is
`½ Nc(Nc−1) · ⟨σ v_rel⟩ · Δt / Vc`. Bird's No-Time-Counter method computes this by
selecting a fixed number of **candidate pairs**,

```
M = ½ · Nc · (Nc − 1) · σ · v_rel,max · Δt / Vc ,
```

and accepting each candidate with probability `σ v_rel / (σ v_rel,max)`. Since our
cross-section `σ` is a constant, the acceptance probability is simply
`|v_rel| / v_rel,max`. Averaged over candidates this recovers the correct rate
`∝ ⟨σ v_rel⟩` while only ever evaluating `M ≈ O(Nc)` pairs — never the full
`O(Nc²)` — so the cost is linear in `N_gas`.

The **collision cross-section** is set by `Gas_Radius` (`r`): two gas particles of
radius `r` collide when their centres pass within `2r`, giving

```
σ = π · (2r)² .
```

So `Gas_Radius` is the single knob on the collision *rate* (rate `∝ σ ∝ r²`); it no
longer gates collisions geometrically as it did in the earlier overlap scheme.

Implementation details:
- `M` is generally non-integer; its fractional part is carried **stochastically**
  (round up with probability equal to the fraction), so the mean is exact even when
  `M < 1` — which is the norm at our small `Δt`.
- A safety cap `M ≤ 4·Nc` guards against a pathologically small `Vc`.
- Candidate partners are drawn uniformly at random and distinct within the cell.
- A fixed-seed `std::mt19937` supplies the random numbers, so runs are
  reproducible.

### 4.4 Note on time-step and detection

The paper's grid scheme and ours both decouple collisions from instantaneous
geometry, which is what makes the method robust to time-step. The paper uses
`Δt = 0.01` outer-disc crossing times (≈ 0.0136 code units); this project uses a
much smaller `Δt = 0.0005` for accurate inner-disc orbits. Because DSMC's collision
count scales with `Δt`, the *physical* collision rate is the same at either
time-step — we simply resolve fewer collisions per step over more steps.

---

## 5. Resolving a collision

For each accepted pair `(i, j)` with relative velocity `g = v_i − v_j`
(`|g| = grel`):

### 5.1 Sampling the line of centres

In a real hard-sphere collision the impulse acts along the **line of centres**
(the "apse line") at contact, whose orientation depends on the impact parameter.
DSMC *samples* this direction from the hard-sphere distribution rather than reading
it from particle positions. For hard spheres the angle `θ` between the apse line
`n̂` and the incoming direction `−ĝ` is distributed as `p(θ) ∝ sinθ cosθ`, so

```
μ = cos θ = √U₁ ,     azimuth ε = 2π U₂ ,     (U₁, U₂ uniform on [0,1))
```

and `n̂` is constructed in an orthonormal frame about `ĝ`. By construction
`g · n̂ = −μ · grel < 0` (the pair is approaching along `n̂`).

**This sampling is the whole reason we use DSMC.** A geometric-overlap test in a
razor-thin disc detects almost exclusively *vertical* encounters (particles are
far apart in-plane but crammed together vertically), so it damps `σ_z` but barely
touches the in-plane `σ_r` — the gas ends up vertically cold but radially warm
(~30–40 km/s), and the arms stay diffuse. Sampling the apse line makes the cooling
**isotropic**: `σ_r` relaxes toward the collisional equilibrium alongside `σ_z`,
reproducing the paper's cool gas disc.

### 5.2 The inelastic impulse

The normal (line-of-centres) component of the relative velocity is reversed and
scaled by the restitution coefficient `α = Gas_Restitution`: post-collision it
becomes `−α (g·n̂)`. This is applied as an equal-and-opposite impulse along `n̂`,
using the reduced mass `m_red = m_i m_j / (m_i + m_j)`:

```
J   = −(1 + α) · (g · n̂) · m_red
v_i +=  (J / m_i) · n̂
v_j −=  (J / m_j) · n̂
```

- `α = 0` (the paper's value, and the default) is **fully inelastic**: the entire
  normal relative velocity is removed, maximising dissipation.
- `α = 1` is elastic (no energy removed).
- The tangential component is untouched.

Because the impulse is equal and opposite, **linear momentum is conserved
exactly** for every collision; energy is removed (for `α < 1`) only from the normal
relative motion, which is the physical dissipation channel.

---

## 6. Per-species gravitational softening

Cold, dissipative, self-gravitating gas is gravitationally (Jeans/Toomre)
unstable. On a *grid* gravity solver like the paper's (48 × 36 × 19 polar grid,
softening `ε = 0.05` length units ≈ 0.93 kpc) sub-cell collapse simply cannot be
represented, so the cooled gas forms smooth density enhancements. On our
high-resolution **Barnes–Hut tree** with small softening, the same cold gas
fragments into dense clumps ("balls") — an artefact of resolving a collapse the
paper's method smooths away. No collision parameter can prevent this: gas falling
into a forming clump reaches relative speeds far above any threshold, and a
gravitationally bound clump is not un-bound by collisions.

The fix is in the *gravity*, via a **per-species softening**. The force on a body
uses Plummer-style softening `a ∝ 1/(|Δx|² + ε²)^{3/2}`, and we let the softening
length depend on the **species of the body being accelerated** (`SoftSq(i)`):

```
ε(i) = Gas_Softening   if  is_gas[i]  and  Gas_Softening > 0
     = r_soft          otherwise
```

Setting `Gas_Softening ≫ r_soft` means the gas feels a potential smoothed over
~kpc, so its self-gravity cannot collapse below `Gas_Softening`, while the
collisionless stars keep the small `r_soft` and retain sharp spiral structure.
This **decouples** the two requirements: sharp stellar arms *and* non-fragmenting
gas.

Properties and caveats of the sink-based scheme:

- **Symmetry.** Because the softening is chosen from the *sink* body, gas–gas and
  star–star interactions use a single softening in both directions and remain
  symmetric (Newton's third law holds). Only **star↔gas cross pairs** are
  asymmetric (a star feels a gas particle at `r_soft`; the gas feels the star at
  `Gas_Softening`), a small, bounded momentum non-conservation acting only at close
  range. A fully symmetric pair rule (e.g. `ε_ij = max(ε_i, ε_j)`) would require
  per-node effective softenings or separate trees and is not currently used.
- **Convention.** `r_soft` and `Gas_Softening` are physical **lengths** in code
  units. The force kernels need `ε²`, so the two possible squared values
  (`r_soft²` and `Gas_Softening²`) are precomputed once (`UpdateSofteningSq`) and
  `SoftSq(i)` is a cached branch-and-read on the hot path.

---

## 7. Parameters

| Script command | Symbol | Meaning | Default |
|---|---|---|---|
| `Gas_Restitution` | `α` | Restitution coefficient; normal relative velocity → `−α ×` prior. `0` = fully inelastic (paper). | `0` |
| `Gas_Radius` | `r` | Gas collision radius (code units); sets `σ = π(2r)²` and hence the collision **rate**. `0.155` = `0.0005·Rd` (paper). | `0.155` |
| `Gas_Cell_Size` | `L` | DSMC collision-cell edge (code units). Locality of partners only; does not affect the rate. | `6.0` |
| `Gas_Softening` | `ε_gas` | Gravitational softening **length** applied to gas sinks; `> r_soft` suppresses gas fragmentation. `0` = gas uses `r_soft`. | `0` |
| `r_soft` | `ε_star` | Gravitational softening **length** for all particles (gas falls back to this when `Gas_Softening = 0`). | — |

Gas particles themselves are declared per disc through the `gas_mass` /
`gas_fraction` fields of `GalaxyDisc`. The whole collision machinery is inactive
unless some system declares gas.

The M51 case (`scripts/M51.sim`) uses the paper-faithful collision parameters
`α = 0`, `r = 0.155`, and relies on `Gas_Softening = 15.5` (≈ 0.93 kpc, the paper's
softening) to keep the gas from fragmenting while the stars stay sharp at
`r_soft = 3.873`.

---

## 8. Relationship to the paper, and limitations

**Faithful to the paper:** the split into 4:1 star:gas, the warm `Q = 1.5`
initialisation, the fully-inelastic (`α = 0`) line-of-centres collision rule, the
gas radius `0.0005 Rd`, the grid-based rate-controlled collision search, and the
emergent `σ_gas ≈ 5–10 km/s` cool disc. The DSMC/NTC formulation is the same class
of method Salo (1991)/SL93 describe.

**Where we necessarily differ / known limitations:**

- **Density estimate.** The collision rate uses the per-cell bounding-box volume as
  a local density proxy. This is an approximation (noisy for sparsely populated
  cells, and it can over- or under-estimate the true occupied volume); the rate can
  be calibrated in practice through `Gas_Radius`.
- **Tree vs grid gravity.** Our tree resolves gas self-gravity the paper's coarse
  grid cannot, hence the need for the per-species `Gas_Softening`. Even with it, the
  gas can concentrate into thin filaments / softening-limited (~kpc) features under
  strong tidal compression — partly physical, partly a residue of the
  sticky-particle-on-a-tree approach.
- **No pressure, no shocks, no driving.** The scheme captures *dissipation* only.
  Real ISM turbulence is continuously *driven* (supernovae, stellar winds,
  instabilities) and maintains a pressure-supporting `σ ≈ 5–10 km/s` equilibrium
  dynamically; here there is no such driving term, so the gas can only cool. There
  is no thermal pressure, radiative heating/cooling network, magnetic field, or star
  formation.
- **Sink-based softening** is not perfectly momentum-conserving for star↔gas pairs
  (§6).
- **State save/load** does not currently persist `is_gas` or `Gas_Softening`; the
  gas model is intended for fresh runs, not restarts.

---

## 9. References

- **Bird, G. A. 1963**, *Physics of Fluids*, 6, 1518 — original Direct Simulation
  Monte Carlo method.
- **Bird, G. A. 1994**, *Molecular Gas Dynamics and the Direct Simulation of Gas
  Flows* (Oxford: Clarendon Press) — the standard DSMC reference, including the
  No-Time-Counter (NTC) collision-pair selection used here.
- **Combes, F. & Gerin, M. 1985**, *A&A*, 150, 327 — early use of sticky-particle
  gas dynamics in a galactic (barred-spiral) context.
- **Salo, H. 1991**, *A&A*, 243, 118 — full description of the sticky-particle
  method referenced by Salo & Laurikainen (2000).
- **Salo, H. & Laurikainen, E. 1993**, *ApJ*, 410, 586 (SL93) — sticky-particle
  treatment of gas in interacting-galaxy modelling.
- **Salo, H. & Laurikainen, E. 2000**, *MNRAS*, 319, 377 — the M51 (NGC 5194/5195)
  model this project reproduces; source of the gas prescription (their Section 2.1).
- **Toomre, A. 1964**, *ApJ*, 139, 1217 — the Toomre `Q` stability parameter used
  for the warm initial velocity dispersion.
- **Toomre, A. & Toomre, J. 1972**, *ApJ*, 178, 623 (TT72) — foundational tidal
  interaction models motivating the M51 study.

---

*Implementation: `Simulation::ProcessGasCollisions()` (collision step) and
`Simulation::SoftSq()` / `UpdateSofteningSq()` (per-species softening) in
`Simulation.cpp` / `Simulation.h`. Script parameters are documented in
`docs/script-format.md`; the M51 application in `docs/compute_M51.md`.*
