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

### Truncation

By default the profile above extends to all radii, so its enclosed mass grows without bound (the price of a flat rotation curve — `M ∝ r`). A halo can be **truncated** at a radius `Rh` via the `halo_truncation_radius` argument on `GalaxyDisc` / `SphericalUniverse` (Salo & Laurikainen 2000 set `Rh` equal to each galaxy's disc truncation `Rd`). Beyond `Rh` the enclosed mass is frozen at its value there and the field reverts to a point-mass `1/r²`:

```
r <= Rh :  a_halo = v_c² * r_vec / (r² + r_c²)          cored isothermal
r >  Rh :  a_halo = M_halo(Rh) * r_vec / r³             enclosed mass frozen
           M_halo(Rh) = v_c² * Rh³ / (Rh² + r_c²)
```

The two branches agree at `r = Rh`, so the force is continuous. This is what `HaloScale()` implements (returning the scale factor `s` with `a = s * r_vec` for either branch). Set `Rh = 0` to leave the halo untruncated. The M51 scripts truncate at `Rd`; the Milky Way / Andromeda scripts currently leave `Rh = 0`.

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

### Halo Centering (inertial halo centres)

The halo is **not** anchored to the central body, and it is **not** re-derived from the particles
each step. Following Salo & Laurikainen (2000) — *"the coordinate grids are centred on the halo
centres and the disc back-action is taken into account in the halo motion"* — each halo centre is a
**dynamical body integrated under gravity**, carrying its own position, velocity and inertial mass
(the total truncated halo mass). `IntegrateHaloCenters()` computes its acceleration every derivative
evaluation, and it is advanced with the same velocity-Verlet step as the particles.

Its acceleration has two contributions:

- **Disc back-reaction.** Every particle the halo pulls, pulls back on it (Newton's third law). The
  net force the halo exerts on all particles is `Σ mᵢ · HaloScale_s · (hc_s − posᵢ)`, and the
  reaction `−(that) / M_halo` accelerates the centre. This is the "disc back-action" the paper
  refers to, and it makes the interaction momentum-conserving.
- **Halo–halo.** Each centre falls in the field of every other halo exactly as a particle there
  would; beyond the truncation radii the pair acts as equal-and-opposite point masses.

The centre therefore moves *because gravity moves it*, not because a centroid was recomputed — so
tidal debris cannot drag the halo off the galaxy core. During warmup the systems are isolated and
held in place, so the centre stays on the relaxing disc; it begins evolving freely at t = 0.

The central body is a **free particle — it is not pinned** to any centroid;
pinning it to the barycentre would drag the nucleus toward tidal debris (the very
artifact the inertial centre is designed to avoid). Its mass is script-dependent:
some setups (e.g. the Milky Way / Andromeda) use it as a real bulge, while the M51
scripts fold all baryon into the disc and reduce it to a token ~1-particle anchor.
Either way it moves under gravity like any other particle.

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
    double halo_scale = HaloScale(sys, rsq);       // cored isothermal, or frozen 1/r^2 beyond Rh
    vscaleadd(r_halo, halo_scale, acc_t[bi]);      // acc += halo_scale * r_vec
}

// Cross-halo: every other system's halo
for (int s = 0; s < N_Systems; s++) {
    if (s == sys || halo_vc[s] == 0.0) continue;
    double *hc2 = &halo_center[s * 3];
    r_halo = hc2 - pos_t[bi];
    rsq = vmagsq(r_halo);
    double halo_scale = HaloScale(s, rsq);
    vscaleadd(r_halo, halo_scale, acc_t[bi]);
}
```

This computes `a = v_c² * r_vec / (r² + r_c²)` toward each halo center (or the truncated `1/r²`
form beyond `Rh`, see *Truncation*). The sign is attractive because `r_halo` points from the body
toward the center. No softening is needed: the `r_c²` term in the denominator already regularizes
`r → 0`.

The **cross-halo term is what makes multi-galaxy encounters work**. Because a halo's field acts
nearly uniformly across the other galaxy at large separation, and because each halo centre moves
with its own galaxy under gravity, the relative acceleration of the pair comes out correct:

```
a_rel = G(M_baryA + M_haloA + M_baryB + M_haloB) / d²
```

That is, the pair's two-body orbit is right, and — because each halo centre is a dynamical object
that recoils under the disc back-reaction (see *Halo Centering* above) — momentum is conserved as
well, so the orbit no longer decays spuriously.

The central body is not specially excluded from the halo force; the `rsq > 1e-10` guard just skips
a body sitting exactly at the centre. Since the free central body sits near, but not exactly on, the halo
centre (both move independently under gravity), it feels a negligible force from its own halo.

### Computational Cost

The halo acceleration is O(N_Systems) per particle — for each halo, one subtraction, one dot
product, one division, one scale-add. With the handful of systems these scripts use this is
effectively O(1), so the added cost is O(N) per timestep against O(N log N) for the tree
traversal. `IntegrateHaloCenters()` adds one more O(N × N_Systems) pass per derivative evaluation to
accumulate the disc back-reaction that drives each halo centre. Both are negligible.

`IntegrateHaloCenters()` runs single-threaded after the parallel force phase: the halo centres are
fixed while particle forces are computed, then updated once from the accumulated reaction.

### Design: Static vs. Live Halo

This is a **rigid analytic halo carried by an inertial centre**. The shape never changes; the
centre moves as a dynamical body under gravity (the mutual halo field plus the disc back-reaction),
not by being recomputed from the particles.

Advantages:
- Zero additional particles (no memory or tree cost)
- Exact force computation (no tree approximation errors)
- Trivially parallelizable (no data dependencies between bodies)
- Produces the correct rotation curve by construction
- Reproduces the correct two-body orbit for an interacting pair (see Force Computation above)

Limitations — these are ordered by how much they actually distort results:

1. **No halo dynamical friction.** This is the dominant error for interacting systems. In reality an
   intruder raises a trailing density wake in the other galaxy's halo, and that overdensity pulls
   backward — Chandrasekhar friction, the primary orbital-decay channel in a merger. A **rigid**
   spherical potential cannot deform or hold that wake, so it exerts **zero halo drag** regardless of
   how its centre moves. The friction that does survive is the tidal braking of the live baryonic
   discs, a minority of the mass. Consequence: **orbital decay is underestimated and pairs sink more
   slowly than in reality** — but the decay that *is* present is now physical, not the spurious
   barycentre-dragging drag the earlier centring scheme produced.

2. **Unbounded mass when left untruncated.** The cored isothermal profile has `M_halo(r) =
   (v_c²/G)(r − r_c·arctan(r/r_c))`, which grows linearly with radius forever — the price a flat
   rotation curve pays (`v² = GM/r = const` requires `M ∝ r`). This is now **optional**: setting
   `halo_truncation_radius = Rh` freezes the enclosed mass beyond `Rh` and reverts the field to
   `1/r²`, mimicking the virial cutoff of a real halo. The M51 scripts use this (`Rh = Rd` for each
   galaxy, following Salo & Laurikainen); the Milky Way / Andromeda scripts currently leave `Rh = 0`
   (untruncated).

   Where it still bites — an *untruncated* run overestimates long-range attraction. In
   `Milky_Way_Andromeda_Collision.sim` the two halos supply ~3e8 code-mass at the initial 3000-unit
   (~180 kpc) separation against ~1.6e7 of baryons, so the encounter is ~95% halo-driven — roughly
   right, since halos really do dominate at these radii. The SIS overshoots the MW's enclosed mass
   by only ~1.5× at 180 kpc, but the overshoot grows linearly with separation, so a setup started
   farther apart attracts too strongly. Setting `Rh ~ r_200` on those scripts would remove it.

3. **No tidal stripping.** `v_c` and the truncation radius `Rh` are read from the script once and
   held constant, so a satellite's halo keeps its full mass and extent no matter how deeply it
   plunges. M51b retains `v_c = 186.6` and its `Rh = Rd` truncation even after passing through
   M51a's disc, when in reality it would lose much of its outer halo on the first passage — so a
   deeply-plunging companion stays more bound than it should. Truncation (2) caps the *initial*
   extent but does not shrink dynamically as mass is stripped.

4. **Centroid slaving (resolved).** An earlier scheme placed each halo at the mass-weighted
   centroid of *all* its system's particles, recomputed every step. Because the disc outweighs the
   nucleus several-to-one, a long asymmetric tidal tail dragged the entire halo off the core, and
   stripped stars kept pulling their original halo toward the companion — a spurious drag that
   spiralled a close pair into an artificial merger. This is now fixed: the halo centre is an
   inertial body integrated under gravity (see *Halo Centering*), so tidal debris no longer moves it.

5. **Momentum conservation (now enforced).** Because the halo centre is inertial, the force it
   exerts on the particles has a third-law partner — the reaction is applied back to the centre — so
   linear momentum is conserved (exactly for particle–halo pairs, and for halo–halo pairs beyond the
   truncation radii). Energy is not strictly conserved, because the rigid shape cannot absorb tidal
   work, but there is no longer the systematic bulk drift the old centring produced.

6. **The post-merger state is not meaningful.** Once the populations mix, the two halo centres
   converge and you have two coincident rigid wells giving `v_c,eff² = v_cA² + v_cB²`, arrived at
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
because the missing physics (halo friction, stripping) *is* the physics of a merger. The practical
boundary: trust the simulation while the halos remain distinguishable and the encounter is a
perturbation; stop trusting it once the halos interpenetrate. Full fidelity there requires live DM
particles; the cheaper partial fixes are already partly in place — halo truncation is available
(see *Truncation*) and the inertial halo centres remove the old spurious drag — leaving a
stripping-driven decay in `v_c`/`Rh` and an explicit Chandrasekhar friction term as the remaining
inexpensive improvements toward realistic merger timescales.

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
