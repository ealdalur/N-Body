# Dark Matter Halo — Cored Isothermal Sphere

## The Problem: Why Gravity Alone Fails

A disc galaxy composed only of baryonic matter (stars and gas) with a central point mass produces a Keplerian rotation curve:

```
v(r) = sqrt(G * M_enclosed / r)
```

At large radii where most of the mass is interior, this gives `v ~ 1/sqrt(r)` — velocity drops off. But observed galaxies (including the Milky Way) have **flat rotation curves**: the orbital velocity remains roughly constant out to many disc scale lengths.

This discrepancy is one of the strongest pieces of evidence for dark matter. The dark matter halo provides an extended mass distribution that dominates at large radii, producing the flat rotation curve that stabilizes the outer disc.

Without a halo, an N-body disc galaxy will:
- Have Keplerian orbits that are too slow at the outer edge
- Be Toomre-unstable (Q < 1) and fragment into clumps
- Fail to sustain coherent spiral structure
- Produce unrealistic tidal tail morphology during interactions

## The Cored Isothermal Sphere Model

### Density Profile

The singular isothermal sphere has density:

```
ρ(r) = v_c² / (4πG * r²)
```

This produces a perfectly flat rotation curve `v(r) = v_c` at all radii, but diverges at r=0. Adding a core radius `r_c` regularizes the center:

```
ρ(r) = v_c² / (4πG * (r² + r_c²))
```

### Enclosed Mass

Integrating the cored density profile over a sphere of radius r:

```
M_halo(r) = (v_c² / G) * (r - r_c * arctan(r/r_c))
```

For `r >> r_c`, this approaches `M_halo ~ v_c² * r / G` (linear growth), which is what produces the flat rotation curve. For `r << r_c`, `M_halo ~ v_c² * r³ / (3G * r_c²)` (uniform-density core).

### Gravitational Potential

The potential of the cored isothermal sphere is:

```
Φ(r) = (v_c² / 2) * ln(r² + r_c²)
```

### Acceleration (Force per Unit Mass)

Taking the gradient of the potential:

```
a_halo = -∇Φ = -v_c² * r_vec / (r² + r_c²)
```

where `r_vec` is the position vector relative to the halo center and `r² = |r_vec|²`. This is the formula we implement. It has the properties:

- **At large r** (`r >> r_c`): `a ~ v_c²/r` → circular velocity = `v_c` (flat rotation curve)
- **At small r** (`r << r_c`): `a ~ v_c² * r / r_c²` → solid body rotation (harmonic potential)
- **At r = r_c**: transition between the two regimes

### Circular Velocity from the Halo

For a particle on a circular orbit, centripetal acceleration equals gravitational pull:

```
v²/r = v_c² * r / (r² + r_c²)
```

Solving for the orbital speed contributed by the halo:

```
v_halo²(r) = v_c² * r² / (r² + r_c²)
```

This rises from 0 at the center, reaches `v_c/sqrt(2)` at `r = r_c`, and asymptotes to `v_c` at large r.

### Total Circular Velocity

The total circular velocity for a particle at radius r combines the baryonic and halo contributions:

```
v_total²(r) = G * M_baryonic(r) / r + v_c² * r² / (r² + r_c²)
```

where `M_baryonic(r)` is the enclosed baryonic mass (central body + disc particles interior to r).

## Parameter Selection

### v_c — Asymptotic Circular Velocity

This sets the "height" of the flat rotation curve. For a galaxy with total baryonic mass `M` and disc radius `R`, the Keplerian edge velocity is:

```
v_kep = sqrt(G * M / R)
```

Choosing `v_c ≈ v_kep` means the halo doubles the circular velocity at the disc edge relative to a no-halo scenario. This is consistent with observed galaxies where dark matter contributes roughly equal mass to baryonic matter within the optical radius, and dominates beyond it.

For our primary galaxy: `v_kep = sqrt(1.0 * 1e7 / 250) = 200`, so `v_c = 200`.

### r_c — Core Radius

The core radius controls where the rotation curve transitions from rising (inner solid-body regime) to flat:

- **Too small** (r_c << inner disc radius): the halo is nearly singular, forces are very strong near the center, and the innermost particles are dominated by the halo rather than the central mass.
- **Too large** (r_c >> R): the halo barely affects the disc because particles never reach the flat regime; the rotation curve is still rising at the disc edge.
- **Optimal**: `r_c` ≈ 10–25% of the disc radius. The rotation curve becomes flat within the disc, stabilizing the outer regions while letting the central mass dominate the bulge.

For our primary galaxy: `r_c = 50` (20% of R=250).

### Mass Ratio

The effective halo mass within radius r is `M_halo(r) ≈ v_c² * r / G` (for r >> r_c). At the disc edge:

```
M_halo(R) ≈ v_c² * R / G = 200² * 250 / 1.0 = 1.0e7
```

This gives a halo-to-baryonic mass ratio of ~1:1 within the disc radius — the halo approximately doubles the gravitating mass. In reality, the halo extends far beyond the disc, so the total halo mass is much larger, but only the mass within the particle orbits affects dynamics.

## Implementation

### Data Structures

Per-system (per-galaxy) halo parameters stored in `Simulation`:

```cpp
std::vector<double> halo_vc;      // [N_Systems]     v_c for each system's halo
std::vector<double> halo_rc_sq;   // [N_Systems]     r_c² (precomputed, avoids multiply in hot loop)
std::vector<int>    halo_central; // [N_Systems]     body index of each system's central mass
std::vector<double> halo_center;  // [N_Systems * 3] current halo center, recomputed each substep
std::vector<int>    body_system;  // [N_Bodies]      maps each body to its parent system
```

`body_system` is assigned once at load time and is **never updated thereafter** — a particle
stripped from one galaxy and captured by the other still counts as a member of its original
system for the rest of the run. See the Limitations section for why this matters.

### Initialization (LoadGalaxyDiscState)

When creating a galaxy disc, the halo parameters are set and system membership is assigned:

```cpp
halo_vc[system] = haloVc;
halo_rc_sq[system] = haloRc * haloRc;
halo_central[system] = sysIdx;

for (int i = 0; i < N_SYSTEM_BODIES[system]; i++)
    body_system[sysIdx + i] = system;
```

The initial orbital velocities account for both the baryonic and halo contributions:

```cpp
double vc_sq = G * m_orbit / r + haloVc * haloVc * r * r / (r * r + haloRc * haloRc);
vm = sqrt(vc_sq);
```

Without this correction, particles initialized at Keplerian velocity would be too slow for the combined potential and would fall inward, causing the disc to collapse on the first few timesteps.

### Halo Centering (ComputeHaloCenters)

The halo is **not** anchored to the central body. `ComputeHaloCenters()` runs at the top of every
derivative evaluation — so four times per RK4 step, once per LeapFrog step — and places each
halo at the mass-weighted centroid of all particles belonging to that system:

```cpp
halo_center[sys] = Σ mᵢ · pos_t[i] / Σ mᵢ        // i over that system's particle range only
```

Two consequences worth being explicit about:

- The centroid is computed from `pos_t`, the integrator's *current* intermediate state, so the
  halo position is consistent within each RK4 substage rather than lagging a full step.
- The halo has **no state of its own**. It carries no position or velocity variable, no inertia,
  and no momentum. It is re-derived from the baryons every substep. "How the halo moves" is
  entirely a statement about how its member particles move.

Under LeapFrog, `PinCentralBodies()` additionally snaps each central body onto this same centroid
(position and velocity), so the central body and the halo center coincide. Under RK4,
`PinCentralBodies()` is **not** called and the central body is free to drift off the centroid.

### Force Computation (CalcAccelRangeOct / CalcAccelRangeP2P)

After computing the N-body gravitational acceleration from the tree (or direct sum), each body
receives its own system's halo acceleration, and then the halo acceleration of every *other*
system:

```cpp
// Own halo
int sys = body_system[bi];
double *hc = &halo_center[sys * 3];
r_halo = hc - pos_t[bi];                           // r_vec points from body toward halo center
double rsq = vmagsq(r_halo);
if (rsq > 1e-10) {                                 // skip a body sitting exactly at the center
    double halo_scale = halo_vc[sys] * halo_vc[sys] / (rsq + halo_rc_sq[sys]);
    vscaleadd(r_halo, halo_scale, acc_t[bi]);      // acc += halo_scale * r_vec
}

// Cross-halo: every other system's halo
for (int s = 0; s < N_Systems; s++) {
    if (s == sys || halo_vc[s] == 0.0) continue;
    double *hc2 = &halo_center[s * 3];
    r_halo = hc2 - pos_t[bi];
    rsq = vmagsq(r_halo);
    double halo_scale = halo_vc[s] * halo_vc[s] / (rsq + halo_rc_sq[s]);
    vscaleadd(r_halo, halo_scale, acc_t[bi]);
}
```

This computes `a = v_c² * r_vec / (r² + r_c²)` toward each halo center. The sign is attractive
because `r_halo` points from the body toward the center. No softening is needed: the `r_c²` term
in the denominator already regularizes `r → 0`.

The **cross-halo term is what makes multi-galaxy encounters work**. Because a halo's field acts
nearly uniformly across the other galaxy at large separation, and because each halo tracks its own
galaxy's centroid, the relative acceleration of the pair comes out correct:

```
a_rel = G(M_baryA + M_haloA + M_baryB + M_haloB) / d²
```

That is, the pair's two-body orbit is right even though neither halo is a dynamical object. This
is the main thing the barycenter-tracking design buys over a halo anchored to a single particle.

The central body is **no longer specially excluded**. The old `bi != halo_central[sys]` test has
been replaced by the `rsq > 1e-10` guard against the halo center. Under LeapFrog these are
equivalent, since `PinCentralBodies()` puts the central body exactly at the centroid. Under RK4
they are not: a central body displaced by `δ` from the centroid feels a restoring acceleration
`v_c² δ / r_c²`, i.e. a harmonic oscillation about the centroid at angular frequency `ω = v_c/r_c`
(≈ 2.5 per code time unit for M51's primary). This is broadly physical — a nucleus does sit in its
halo's potential — but note the restoring force points at the *particle centroid*, which includes
tidal debris, not at the visible nucleus.

### Computational Cost

The halo acceleration is O(N_Systems) per particle — for each halo, one subtraction, one dot
product, one division, one scale-add. With the handful of systems these scripts use this is
effectively O(1), so the added cost is O(N) per timestep against O(N log N) for the tree
traversal. `ComputeHaloCenters()` adds one more O(N) pass per derivative evaluation. Both are
negligible.

Note that `ComputeHaloCenters()` runs single-threaded between the two parallel phases of
`CalcDerivatives`, since every worker thread needs all halo centers finalized before any body's
force is computed.

### Design: Static vs. Live Halo

This is a **rigid analytic halo whose center is slaved to its own system's baryonic barycenter**.
The shape never changes; only the center moves, and it moves because its member particles moved.

Advantages:
- Zero additional particles (no memory or tree cost)
- Exact force computation (no tree approximation errors)
- Trivially parallelizable (no data dependencies between bodies)
- Produces the correct rotation curve by construction
- Reproduces the correct two-body orbit for an interacting pair (see Force Computation above)

Limitations — these are ordered by how much they actually distort results:

1. **No dynamical friction.** This is the dominant error for interacting systems. In reality an
   intruder raises a trailing density wake in the other galaxy's halo, and that overdensity pulls
   backward — Chandrasekhar friction, which is the primary orbital-decay channel in a merger.
   A rigid, spherical potential comoving with its own source is symmetric by construction and
   exerts **zero net drag**. Some friction survives from the live baryonic particles, but in these
   scripts the baryons are a minority of the mass. Consequence: **orbital decay is far too slow
   and pairs do not sink on realistic timescales.**

2. **Infinite extent, no truncation.** `M_halo(r) = (v_c²/G)(r − r_c·arctan(r/r_c))` grows without
   bound — mass increases linearly with radius forever, with no virial radius at which the halo
   ends. This is not a coding oversight; it is forced by the profile. A flat rotation curve
   *requires* `v² = GM(r)/r = const`, hence `M ∝ r`. The isothermal sphere buys its flat rotation
   curve by paying with unbounded mass. Real halos stop near the virial radius r_200, beyond which
   the field reverts to `1/r²`.

   Magnitude of the error: in `Milky_Way_Andromeda_Collision.sim` the two halos supply ~3e8
   code-mass at the initial 3000-unit separation against ~1.6e7 of baryons, so the encounter is
   ~95% halo-driven. That fraction is roughly correct — halos really do dominate at these radii,
   and it should not by itself be read as an error. The actual overestimate is milder: 3000 code
   units is ~180 kpc, still near the real virial radius, where the SIS overshoots the MW's
   enclosed mass by only ~1.5x. The overshoot grows linearly with separation, so it is modest for
   these scripts and becomes severe for any setup started farther apart. The more damaging
   companion problem is (3) — this mass can never be removed once it is there.

3. **No tidal stripping.** `v_c` is read from the script once and is constant for all time, so a
   satellite's halo survives intact no matter how deeply it plunges. M51b retains `v_c = 130` even
   after passing through M51a's disc, when in reality it would lose most of its outer halo on the
   first passage. Note this error and (2) push in opposite directions — too much binding at large
   radius, too much survival at small radius — so they do not usefully cancel.

4. **Centroid slaving becomes a feedback pathology once tails form.** The centroid is taken over
   *all* of a system's particles, including debris, and `body_system` is never reassigned. In
   `M51.sim` the primary has `Mfrac = 4.90`, so the disc outweighs the central body ~5:1 and the
   disc particles — not the nucleus — dominate the centroid. A long asymmetric tidal tail therefore
   drags the entire halo off the nucleus, and stars stripped by the companion keep pulling their
   original halo toward the companion for the rest of the run.

5. **Conservation laws are violated.** The halo acts on particles and nothing acts back on it, so
   there is no third-law partner and linear momentum is not conserved. Because the center moves,
   `Φ` is explicitly time-dependent and energy is not conserved either. Both are negligible for an
   isolated galaxy (the centroid barely moves) and grow precisely during close encounters. Note
   the code currently computes no energy or momentum diagnostic — `Data_Log` writes only `pos_sq`
   — so this drift is entirely unmonitored. If such a diagnostic is ever added, it must include
   the halo potential `Φ = (v_c²/2)ln(r² + r_c²)` for every system, and should be expected to
   drift at pericenter for reasons that are *not* integrator error.

6. **The post-merger state is not meaningful.** Once the populations mix, both centroids converge
   and you have two coincident rigid wells giving `v_c,eff² = v_cA² + v_cB²`, arrived at
   instantaneously with no violent relaxation and no mass loss.

7. **Other missing physics.** The halo is spherical and non-rotating, so there is no halo spin, no
   disc/halo angular momentum exchange (bars are never slowed by halo friction), and no substructure.

**Regime of validity.** For flybys and prograde tidal encounters — the M51 case this code was tuned
for — the approximation is sound and the morphology is trustworthy. Tidal tails, bridges, and
grand-design arms are the *disc's* response to an external tidal field, and both that field and the
trajectory are approximately right over a single passage. The halo's real job there is deepening
the well ~2× so that κ rises, Toomre Q clears 1, and the disc sustains coherent structure; a rigid
analytic well does that perfectly well.

For **mergers** the model fails in a way that changes the outcome rather than merely blurring it,
because the missing physics (friction, stripping) *is* the physics of a merger. The practical
boundary: trust the simulation while the halos remain distinguishable and the encounter is a
perturbation; stop trusting it once the halos interpenetrate. Full fidelity there requires live DM
particles, though adding a truncation radius and a stripping-driven decay in `v_c` would be far
cheaper and would move merger timescales substantially toward reality.

## Physical Significance

### Toomre Stability (Q Parameter)

The Toomre stability parameter for a stellar disc is:

```
Q = σ_r * κ / (3.36 * G * Σ)
```

where `σ_r` is the radial velocity dispersion, `κ` is the epicyclic frequency, and `Σ` is the surface density. The halo increases `κ` (steeper effective potential) without increasing `Σ`, raising Q above 1 and preventing fragmentation. This is why a halo-stabilized disc can maintain coherent spiral structure rather than breaking into clumps.

### Swing Amplification

Spiral arms in disc galaxies form through swing amplification — leading density waves are sheared into trailing arms by differential rotation, amplified by self-gravity. The amplification factor depends on the parameter X:

```
X = κ² * r / (2π * G * Σ * m)
```

where m is the number of spiral arms. The halo modifies κ, shifting X into the range (1 < X < 3) where swing amplification is most efficient for m=2 (two-armed spirals). Without the halo, the disc either fragments (X too small) or is too hot for amplification (X too large after artificial heating).

### Tidal Interaction Morphology

During a galaxy flyby, the tidal force from the companion stretches the primary disc. The depth of the halo potential well determines:
- How far tidal tails extend (deeper well → shorter tails)
- The bridge morphology between galaxies
- Whether the encounter is bound (merger) or unbound (flyby)

The halo makes the potential well ~2× deeper than baryonic mass alone, producing more realistic tidal tail lengths consistent with observations of M51 and similar interacting systems.
