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
double halo_vc[N_SYSTEMS];      // v_c for each system's halo
double halo_rc_sq[N_SYSTEMS];   // r_c² (precomputed, avoids multiply in hot loop)
int halo_central[N_SYSTEMS];    // body index of each system's central mass
int body_system[N_BODIES];      // maps each body to its parent system
```

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

### Force Computation (CalcAccelRangeOct / CalcAccelRangeP2P)

After computing the N-body gravitational acceleration from the tree (or direct sum), the halo acceleration is added:

```cpp
int sys = body_system[bi];
int ci = halo_central[sys];
if (bi != ci) {
    vsub(pos_t[ci], pos_t[bi], r_halo);       // r_vec = pos_central - pos_body
    double rsq = vmagsq(r_halo);               // |r_vec|²
    double halo_scale = halo_vc[sys] * halo_vc[sys] / (rsq + halo_rc_sq[sys]);
    vscaleadd(r_halo, halo_scale, acc_t[bi]);  // acc += halo_scale * r_vec
}
```

This computes `a = v_c² * r_vec / (r² + r_c²)`, directed toward the central body. The sign is attractive because `r_halo` points from the body toward the center.

The central body itself is excluded (`bi != ci`) — it does not feel its own halo.

### Computational Cost

The halo acceleration is O(1) per particle — one subtraction, one dot product, one division, one scale-add. For N particles this adds O(N) work per timestep, compared to O(N log N) for the tree traversal. The overhead is negligible.

### Design: Static vs. Live Halo

This is a **static analytic halo anchored to the central body**. It moves rigidly with the central mass, which itself is a live N-body particle influenced by all other particles.

Advantages:
- Zero additional particles (no memory or tree cost)
- Exact force computation (no tree approximation errors)
- Trivially parallelizable (no data dependencies between bodies)
- Produces the correct rotation curve by construction

Limitations:
- The halo cannot deform tidally during interactions
- No dynamical friction from halo-halo particle scattering
- The halo has infinite extent (no truncation radius)

For M51-type flyby encounters, these limitations are minor — the defining morphological features (tidal tails, bridges, grand-design spiral arms) are governed by the disc response within the halo potential, not by halo deformation itself. A live halo with particles would be needed for full mergers where the halos physically interpenetrate and merge.

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
