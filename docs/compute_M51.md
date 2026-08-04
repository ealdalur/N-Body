# compute_M51.py

## 1. The M51 System: Overview

The Whirlpool Galaxy (M51, NGC 5194/5195) is the archetypal example of a grand-design spiral galaxy shaped by a tidal interaction. It consists of:

- **NGC 5194 (M51a)**: A large Sbc spiral galaxy with prominent two-arm spiral structure extending ~11 kpc from the nucleus.
- **NGC 5195 (M51b)**: A compact SB0-pec lenticular galaxy that has passed through or near M51a's disc at least once, exciting the grand-design spiral arms.

The M51 system is the clearest known case where spiral structure is demonstrably caused by a companion interaction rather than internal instabilities. The two-arm pattern is coherent across the entire disc, extending from ~2 kpc to the outer edge — a signature of external tidal forcing rather than swing amplification.

### 1.1 Why M51 Matters

Toomre & Toomre (1972) first proposed that galactic spiral arms could be tidally excited by companion passages. M51 is the strongest observational confirmation of this theory: the arm morphology, bridge connecting the two galaxies, and tidal tails all match predictions of a prograde encounter.

### 1.2 Observational Constraints on the Interaction

From N-body modeling (Salo & Laurikainen 2000; Dobbs et al. 2010) and observations:

| Parameter | Value | Source |
|-----------|-------|--------|
| Encounter type | Prograde | Salo & Laurikainen 2000 |
| Pericenter distance (3D) | 15-25 kpc | Salo & Laurikainen 2000 |
| Projected separation (current) | 10-12 kpc | Observed on sky |
| Time since first passage | 300-400 Myr | Salo & Laurikainen 2000 |
| M51b position | Behind M51a disc (farther from Earth) | Dobbs et al. 2010 |
| Orbital inclination to disc | 10-20 degrees | Salo & Laurikainen 2000 |
| Mass ratio (dynamical) | 1:3 to 1:5 | Salo & Laurikainen 2000 |
| Distance to system | 8.0-8.6 Mpc | McQuinn et al. 2016 |

Note: The Salo & Laurikainen range of 15-25 kpc is the **3D pericenter distance** (closest approach in the orbit). The 10-12 kpc **projected separation** observed today is the sky-plane distance at the *current* orbital phase (post-pericenter, M51b has moved to a different position). These are not contradictory — the current projected separation can be smaller than the 3D pericenter because (a) projection removes one spatial dimension and (b) M51b's post-pericenter position differs from the pericenter point.

Our simulation uses 12 kpc pericenter, at the low end of the Salo & Laurikainen range, which produces the best visual match to the observed bridge morphology.

### 1.3 Why Prograde?

The prograde nature of the encounter (companion orbiting in the same sense as disc rotation) is critical for producing strong spiral arms:

- **Prograde encounter**: The companion's tidal field rotates at a rate comparable to disc material's orbital frequency. This produces a near-resonant perturbation that lasts many orbital periods, creating strong, coherent two-arm spirals.
- **Retrograde encounter**: The tidal field sweeps past disc material rapidly (in the counter-rotating sense). The perturbation averages out over less than one orbit, producing only weak, transient disturbances.

This is analogous to pushing a swing: pushing in phase with the motion (prograde) builds amplitude; pushing against the motion (retrograde) does little.

---

## 2. Mass Decomposition Methodology

### 2.1 Rotation Curve Decomposition

Following Salo & Laurikainen (2000) and standard practice for N-body tidal interaction studies, the total observed rotation velocity is decomposed into baryonic and halo contributions:

```
V_total^2 = V_baryon^2 + V_halo^2
```

The key constraint from Salo & Laurikainen (2000) and Dobbs et al. (2010) is that the dark matter halo encloses approximately **twice** the baryonic mass within the disc radius (M_halo/M_baryon ~ 2). This means:

```
V_halo^2 / V_total^2 ~ 2/3    (halo provides 2/3 of centripetal support)
V_baryon^2 / V_total^2 ~ 1/3   (baryons provide 1/3)

Therefore:
V_baryon = V_total / sqrt(3)
V_halo   = V_total * sqrt(2/3)
```

### 2.2 Baryonic Mass vs Photometric Stellar Mass

The baryonic masses used here are **rotation-curve-decomposed** values — the dynamically cold disc mass that participates in spiral arm formation. These are deliberately lower than photometric stellar masses (e.g., Querejeta et al. 2015 gives M_stellar ~ 5x10^10 Msun for M51a):

- Photometric mass includes the dynamically hot thick disc, which acts more like a bulge/halo (doesn't respond to tidal perturbation as a cold thin disc)
- The "disc mass" for tidal interaction modeling is the cold, thin, responsive component
- This is standard practice: Salo & Laurikainen (2000) and Dobbs et al. (2010) both use reduced disc masses relative to photometric estimates

### 2.3 Halo Velocity Correction

The simulator uses a **cored isothermal** halo:
```
V_halo(r) = haloVc * r / sqrt(r^2 + Rc^2)
```

At the disc edge (r = R), the circular velocity of the halo is:
```
V_halo(R) = haloVc * R / sqrt(R^2 + Rc^2)
```

We want V_halo(R) to equal V_total * sqrt(2/3), so:
```
haloVc = V_total * sqrt(2/3) / (R / sqrt(R^2 + Rc^2))
       = V_total * sqrt(2/3) * sqrt(1 + (Rc/R)^2)
```

For small Rc/R (which is the case here — Rc << R), haloVc ~ V_total * sqrt(2/3).

---

## 3. Galaxy Parameters

### 3.1 NGC 5194 (M51a)

| Physical Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | Sbc | de Vaucouleurs et al. | — |
| Stellar disc radius (R25) | 11.2 kpc | NED | 186.7 |
| Observed flat V_rotation | 210 km/s | Sofue et al. 1999; Meidt et al. 2013 | — |
| V_baryon (at disc edge) | 121 km/s | = 210/sqrt(3) | — |
| V_halo (at disc edge) | 172 km/s | = 210*sqrt(2/3) | — |
| Total baryonic mass | 2.74 x 10^10 Msun | = V_b^2 * R / G | 2,744,000 |
| Bulge mass (M) | 4.7 x 10^9 Msun | 17% of baryonic (Sbc) | 470,000 |
| Disc mass | 2.27 x 10^10 Msun | baryonic - bulge | 2,274,000 |
| Mfrac (disc/bulge) | 4.84 | = disc/bulge | — |
| DM halo Vc (haloVc) | 172.1 km/s | corrected for core | — |
| DM halo core radius (haloRc) | 1.0 kpc | Salo: ~0.5-1 kpc | 16.7 |
| M_halo within disc | 5.49 x 10^10 Msun | cored isothermal | 5,488,000 |
| M_halo / M_baryon | 2.00 | Target: ~2 | — |
| Inner hole (bulge region) | 0.3 kpc | — | 5.0 |
| Toomre Q | 1.5 | Warm disc (tidal response) | — |

### 3.2 NGC 5195 (M51b)

| Physical Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | SB0-pec | de Vaucouleurs et al. | — |
| Optical radius | 2.1 kpc | NED | 35.0 |
| Observed V_rotation | 130 km/s | mass-Vc scaling | — |
| V_baryon (at disc edge) | 75.1 km/s | = 130/sqrt(3) | — |
| V_halo (at disc edge) | 106.1 km/s | = 130*sqrt(2/3) | — |
| Total baryonic mass | 1.97 x 10^9 Msun | = V_b^2 * R / G | 197,167 |
| Bulge mass (M) | 8.0 x 10^8 Msun | 40% of baryonic (SB0) | 80,000 |
| Disc mass | 1.17 x 10^9 Msun | baryonic - bulge | 117,167 |
| Mfrac (disc/bulge) | 1.46 | = disc/bulge | — |
| DM halo Vc (haloVc) | 110.4 km/s | corrected for core | — |
| DM halo core radius (haloRc) | 0.6 kpc | compact galaxy | 10.0 |
| M_halo within disc | 3.94 x 10^9 Msun | cored isothermal | 394,333 |
| M_halo / M_baryon | 2.00 | Target: ~2 | — |
| Inner hole | 0.2 kpc | — | 3.3 |
| Toomre Q | 1.5 | Warm disc (tidal response) | — |

### 3.3 Mass Ratio

The dynamical mass ratio at the encounter distance determines the tidal interaction strength. At the pericenter distance (12 kpc = 200 code units), using the cored isothermal formula M_halo(r) = Vc^2 * r^3 / (r^2 + Rc^2):

```
M51a at 12 kpc:
  baryonic:     2,744,000
  halo:         172.1^2 * 200^3 / (200^2 + 16.7^2) = 5,886,000
  total:        8,630,000

M51b at 12 kpc:
  baryonic:     197,000
  halo:         110.4^2 * 200^3 / (200^2 + 10.0^2) = 2,431,000
  total:        2,628,000

Dynamical ratio: 2,628,000 / 8,630,000 = 1 : 3.3
```

This is within the observational range (1:3 to 1:5 from Salo & Laurikainen 2000).

---

## 4. Orbital Computation

### 4.1 Strategy

We place M51b at **apocenter** (the farthest point of its bound orbit) at t=0, with all velocity tangential (v_r = 0). The simulation then runs forward as M51b falls inward toward pericenter, exciting spiral arms in M51a.

This approach has two advantages:
1. Starting at apocenter means the initial velocity is purely tangential — no radial component to specify.
2. The spiral arm excitation happens gradually during infall, producing a realistic build-up of arm strength.

### 4.2 Gravitational Potential

The simulation applies cored isothermal halo accelerations:
```
a_halo = Vc^2 * r / (r^2 + Rc^2)    (attractive, directed inward)
```

The corresponding gravitational potential (with a_r = -dPhi/dr):
```
Phi_halo(r) = +0.5 * Vc^2 * ln(r^2 + Rc^2)
```

This INCREASES outward (a particle at larger r has higher potential energy and falls inward — gravity is attractive). For baryonic point masses:
```
Phi_baryon(r) = -G*M/r
```

The total potential for the mutual orbit includes both halos plus the baryonic masses:
```
Phi(r) = +0.5*Vc_a^2*ln(r^2 + Rc_a^2) + 0.5*Vc_b^2*ln(r^2 + Rc_b^2) - M_baryon/r
```

where M_baryon = M_baryon_a + M_baryon_b = 2,744,000 + 197,167 = 2,941,167.

**Important**: Using a simpler -M_eff/r (Keplerian) potential is incorrect for logarithmic halos. The potential is NOT -M/r; it is a sum of logarithmic terms. The Keplerian approximation significantly underestimates the orbital velocity (by a factor of ~2 for these parameters) because it doesn't account for the potential's slower radial variation.

### 4.3 Energy and Angular Momentum Conservation

At apocenter (r = r_start, v_r = 0):
```
E = (1/2) * v_t^2 + Phi(r_start)
L = r_start * v_t
```

At pericenter (r = r_peri, v_r = 0):
```
E = (1/2) * v_p^2 + Phi(r_peri)
L = r_peri * v_p
```

From angular momentum conservation: `v_p = r_start * v_t / r_peri`

Substituting and solving for v_t:
```
v_t^2 = 2 * (Phi_peri - Phi_start) / (1 - (r_start/r_peri)^2)
```

Both numerator and denominator are negative:
- Phi_peri < Phi_start (pericenter is deeper in the potential well)
- 1 - (r_start/r_peri)^2 < 0 (since r_start > r_peri)

So v_t^2 is positive (negative / negative = positive).

### 4.4 Computed Orbital Parameters

| Parameter | Value |
|-----------|-------|
| Initial separation (apocenter) | 667 code units (40 kpc) |
| Pericenter distance | 200 code units (12 kpc) |
| Phi(r_start) | 267,530 |
| Phi(r_peri) | 206,992 |
| Delta Phi | -60,538 |
| Tangential velocity at apocenter | 109.4 km/s |
| Velocity at pericenter | 364.8 km/s |
| Specific orbital energy | 273,517 |
| Angular momentum L | 72,952 |
| v_circular at apocenter | 215.0 km/s |
| v_t / v_circ | 0.509 |
| Half-orbit time (to pericenter) | ~6 code time units (~363 Myr) |
| Full orbital period | ~12 code time units (~725 Myr) |

The ~363 Myr half-orbit time is consistent with Salo & Laurikainen's (2000) 300-400 Myr between passages.

### 4.5 Coordinate System and Orbital Inclination

- M51a disc lies in the **x-z plane** with disc normal along **+y**
- M51a disc rotates **clockwise** when viewed from +y (because `LoadGalaxyDiscState` uses `v_tan = -vm`, and the tangential direction `t_hat = cross(r_hat, n_hat)` is +z at the +x position, so particles at +x move in -z = CW from +y)
- M51b starts along the **+x axis** at 644 code units (40 kpc * cos(15deg))
- For a **prograde** encounter, M51b must orbit in the same sense as M51a's disc (CW from +y), so its tangential velocity at +x is in the **-z direction**
- Orbital inclination: 15 deg tilt from the disc plane, giving a small **-y velocity** component

Final initial conditions:
```
M51b position: (644.0, 0.0, 0.0)
M51b velocity: (0.0, -28.3, -105.7) km/s
  vel_y = -v_t * sin(15deg) = -109.4 * 0.2588 = -28.3
  vel_z = -v_t * cos(15deg) = -109.4 * 0.9659 = -105.7
|velocity| = 109.4 km/s
```

The negative signs ensure the orbit is prograde (same sense as M51a's disc rotation).

---

## 5. Dark Matter Halos

### 5.1 Role in the Interaction

The dark matter halos play one critical role here, and conspicuously fail to play a second one.

1. **Tidal mass** (modeled): The effective gravitating mass at the interaction distance is dominated by the halos (Vc^2 * r^3/(r^2+Rc^2) >> M_baryonic for r > 100 code units). This sets the orbital velocity and the encounter timescale, and the cross-halo term in the code carries it correctly.

2. **Dynamical friction** (NOT modeled): In reality, as M51b moves through M51a's extended halo it raises a trailing density wake and loses orbital energy to gravitational drag (Chandrasekhar friction). **The simulation does not reproduce this.** The analytic halos are rigid spherical potentials comoving with their own galaxy's barycenter, and such a potential is symmetric by construction — it exerts no net drag. The cross-halo term is a conservative central force and does no secular work on the orbit. The only friction present comes from the live baryonic particles, which are a minority of the mass.

   Consequence for this script: orbital decay is badly underestimated, so M51b will not sink on a realistic timescale. The *first passage* — which is what this simulation is set up to reproduce — is unaffected, since friction has had no time to act. Do not trust the post-encounter orbit, the return time, or any eventual merger.

### 5.2 Halo Parameters

The simulator uses a **cored isothermal** halo profile:
```
a_halo = Vc^2 * r / (r^2 + Rc^2)
```

This gives:
- At r >> Rc: flat rotation curve with V = Vc
- At r << Rc: linear rotation (solid body), V = Vc * r / Rc
- Core radius Rc prevents singular density at the center

| Galaxy | haloVc (km/s) | haloRc (code units) | haloRc (kpc) | M_halo/M_baryon |
|--------|---------------|---------------------|--------------|-----------------|
| M51a | 172.1 | 16.7 | 1.0 | 2.00 |
| M51b | 110.4 | 10.0 | 0.6 | 2.00 |

Note: haloVc is NOT the total observed rotation velocity. It is the halo-only contribution, corrected for the core radius so that V_halo(R_disc) = V_total * sqrt(2/3). The total rotation curve results from the quadrature sum of baryonic and halo contributions.

---

## 6. Velocity Dispersion and Disc Stability

### 6.1 Toomre Q Parameter

Both galaxies use a target Toomre Q = 1.5. This controls the radial and tangential velocity dispersion of disc particles via the Toomre stability criterion:

```
sigma_r(r) = Q * 3.36 * G * Sigma(r) / kappa(r)
sigma_phi(r) = sigma_r * kappa / (2 * Omega)
```

where Sigma(r) is the local exponential surface density and kappa(r) is the epicyclic frequency computed from the full rotation curve (baryons + halo).

Q = 1.5 produces a "warm" disc that:
- Is stable against spontaneous fragmentation and particle-noise-driven multi-arm spirals
- Remains responsive to the strong m=2 tidal perturbation from the companion

See `docs/toomre-q-velocity-dispersion.md` for the full derivation and implementation details.

### 6.2 Why Not Lower Q?

At Q < 1.2, particle noise drives incoherent modes (m=2, 3, 4 all comparable) that create spurious multi-arm spiral structure before the encounter even begins. At Q = 1.5, these noise-driven modes are suppressed while the coherent tidal forcing (which is much stronger than noise) still produces clear two-arm spirals.

---

## 7. Particle Distribution

### 7.1 Counts and Resolution

| Galaxy | Particles | Purpose |
|--------|-----------|---------|
| M51a | 640,000 (639,999 disc + 1 central) | High resolution for spiral arm structure |
| M51b | 160,000 (159,999 disc + 1 central) | Adequate for compact companion |
| Total | 800,000 | |

M51a gets 4x more particles because:
- It is the galaxy where spiral arms must be resolved
- Its disc is ~5x larger in radius (much larger area to sample)
- Spiral arm contrast requires adequate particle density to be visible

### 7.2 Softening and Particle Count

The gravitational softening length should scale as 1/sqrt(N) to maintain consistent two-body relaxation:

| N_M51a | r_soft | Motivation |
|--------|--------|------------|
| 40,000 | 0.3 | Low resolution |
| 160,000 | 0.15 | Medium resolution |
| 640,000 | 0.075 | High resolution |

The M51.sim script uses r_soft = 0.3. When tuning for quality, reducing r_soft in proportion to 1/sqrt(N) improves force resolution.

### 7.3 Disc Generation

The `GalaxyDisc` generator creates particles in a disc with:
- Surface density proportional to 1/r (from radius Ri to R)
- Circular orbital velocity from enclosed mass + halo contribution:
  `v_c^2 = G * M_enclosed/r + haloVc^2 * r^2/(r^2 + haloRc^2)`
- Velocity dispersion set by Toomre Q (radial and tangential, varying with radius)
- Disc orientation set by the normal vector

For M51b, the disc normal is tilted 15 deg from +y to align roughly with its orbital plane:
```
normal = (sin(15deg), cos(15deg), 0) = (0.2588, 0.9659, 0.0)
```

---

## 8. Expected Simulation Behavior

### 8.1 Timeline

| Phase | Sim time | Physical time | What happens |
|-------|----------|---------------|--------------|
| Approach | t = 0-3 | 0-175 Myr | M51b falls inward. Tidal perturbation of M51a's outer disc begins. |
| Pre-pericenter | t = 3-6 | 175-350 Myr | Spiral arms grow stronger as M51b accelerates inward. |
| Pericenter | t ~ 6 | ~360 Myr | Closest approach (12 kpc). Strong tidal torque on disc material. Two-arm spiral pattern is strongly excited. |
| Recession | t = 6-10 | 360-590 Myr | M51b recedes. Grand-design spiral arms are fully developed. Tidal bridge/tail connects the galaxies. |
| Second passage | t ~ 12 | ~725 Myr | M51b returns. Arms may be reinforced or disrupted depending on phase alignment. |

### 8.2 What to Look For

1. **Two-arm spiral arms**: Coherent, symmetric spiral arms extending from M51a's inner disc to its outer edge. These should be visible from t ~ 5 onward.

2. **Tidal bridge**: A stream of particles connecting M51a and M51b, pulled out of M51a's disc by M51b's gravity during closest approach.

3. **Tidal tail**: An arm extending in the opposite direction from the bridge (anti-companion side), formed by angular momentum transfer.

4. **M51b compression**: The compact companion may develop a slightly distorted morphology as it passes through M51a's tidal field.

5. **Arm winding**: After pericenter, the inner parts of the spiral arms rotate faster than the outer parts (differential rotation), causing progressive winding. The pattern should remain open for ~2-3 code time units after pericenter before winding significantly.

### 8.3 Why the Arms Form

The physics of tidal spiral arm formation:

1. M51b's gravitational field creates a **tidal quadrupole** on M51a's disc — stretching it along the line connecting the two galaxies and compressing it perpendicular to that line.

2. Disc material on the **near side** of M51a (facing M51b) experiences stronger attraction toward M51b than the disc center does. It is pulled outward.

3. Disc material on the **far side** experiences weaker attraction than the center. Relative to the center, it appears to be pushed outward (tidal effect).

4. Because the encounter is **prograde**, this quadrupole rotates at roughly the same rate as disc material orbits. The perturbation is coherent over multiple orbital periods, building up a strong two-arm response.

5. The **spiral** shape arises because the inner disc orbits faster than the outer disc (differential rotation). A linear perturbation at different radii gets sheared into a trailing spiral by differential rotation.

---

## 9. Comparison to Published Models

### 9.1 Salo & Laurikainen 2000

Their best-fit model parameters:
- Pericenter: 20-25 kpc
- Mass ratio: 1:3
- M_halo/M_disc: ~2
- Orbital inclination: 10-20 deg
- Time since first passage: 350-400 Myr
- Prograde encounter

Our parameters: pericenter 12 kpc, dynamical mass ratio 1:3.3 (at pericenter), M_halo/M_baryon = 2.0, inclination 15 deg, time-to-pericenter ~363 Myr. Good agreement on mass decomposition, encounter geometry, and timescale. Our closer pericenter (below their published range) was chosen empirically to reproduce the observed bridge morphology at 800k particles — the tidal bridge connecting M51b to M51a's arm tip requires M51b to pass close enough to pull material from M51a's outer disc.

### 9.2 Dobbs et al. 2010

Their SPH + N-body model adds gas physics (ISM, star formation) on top of the gravitational interaction. They find:
- M_halo/M_disc: ~2.46 (lowered Evans model)
- Spiral arms visible ~200 Myr before pericenter
- Strongest arm-interarm contrast at ~100-200 Myr after pericenter
- Pericenter 25 kpc, mass ratio 1:3

Our collisionless (gravity-only) simulation should reproduce the stellar spiral arm morphology but not gas features (HII regions, dust lanes). Our M_halo/M_baryon = 2.0 is consistent with both papers' range of 2-2.5.

### 9.3 Limitations of This Simulation

1. **No gas physics**: Real M51 has ~25% of its disc mass in gas, which responds more strongly to tidal compression (gas shocks, star formation in arms). Our disc is collisionless.

2. **Simplified dark matter halos**: We use analytic cored-isothermal halos rather than live DM particle halos. The halos are rigid spherical backgrounds whose centers track their own galaxy's particle barycenter. This means no self-consistent halo response, no tidal stripping (M51b keeps its halo Vc even after plunging through M51a's disc), and no dynamical friction (section 5.1). The halos also have infinite extent — `M_halo(r) ~ Vc^2*r` never stops growing — so long-range attraction is overestimated relative to a truncated halo. Additionally, because each halo center is the mass-weighted centroid of *all* its member particles and system membership is never reassigned, a long tidal tail drags the halo center off the nucleus, and stars stripped by the companion keep pulling their original halo toward it. This is a first-passage simulation, which is the regime where these limitations are least damaging. See `docs/dark-matter-halo.md` for the full accounting.

3. **No second passage history**: The real M51 may have undergone multiple passages. We simulate from the pre-first-encounter state.

4. **Disc thickness**: The particle generator creates thin discs. Real galaxy discs have finite scale heights (especially after tidal heating). The 2D projection will look correct but the 3D structure is idealized.

---

## 10. Simulation Parameters

### 10.1 Time Step and Softening

- **dt = 0.0005**: Required for the fast orbital velocities near M51a's center (v_circular ~ 210 km/s at ~83 code units). At r=83, orbital period ~ 2*pi*83/210 ~ 2.5 code time units, so dt/T ~ 0.0002 (well-resolved).

- **r_soft = 0.3**: Gravitational softening (18 pc). Prevents close-encounter singularities. The minimum physical scale resolved is ~2*r_soft = 36 pc.

- **BH_Opening_Theta = 0.5**: Barnes-Hut opening angle. A value of 0.5 gives good accuracy for the tidal interaction (force errors < 1%).

### 10.2 PinCentralBodies

The simulation uses `PinCentralBodies`, which pins each galaxy's central body to its system's center of mass (both position and velocity). It is applied every step to every system with `halo_vc > 0`. It prevents:
- The central body from wandering away from the galaxy due to N-body noise
- Unphysical recoil from discrete particle encounters

The central body still feels and transmits the tidal field from the other galaxy.

### 10.3 Cross-Halo Gravity

The code includes cross-system halo gravity (the cross-halo loop in `CalcAccelRangeOct` and `CalcAccelRangeP2P`). Every body feels its own system's halo plus the halo of every other system. This means:
- M51a particles feel M51b's dark matter halo potential
- M51b particles feel M51a's dark matter halo potential
- This creates the correct tidal field for spiral arm excitation

Each halo is centered on its own system's mass-weighted particle barycenter, recomputed every derivative evaluation by `ComputeHaloCenters()` — not anchored to the central body. Because a halo's field is nearly uniform across the other galaxy at these separations, and because each halo tracks its own galaxy, the pair's relative acceleration comes out correct at `G(M_totA + M_totB)/d^2` including halo mass. This is why the analytic orbit derived in section 4 is reproduced by the simulation.

What the cross-halo term does **not** do is produce dynamical friction — see section 5.1. It is a conservative central force between two rigid symmetric potentials.

### 10.4 Halo Monopole Removal

A rigid analytic halo has no inertia, so it cannot obey Newton's third law: particles are pulled toward the halo center and nothing pulls back. For an axisymmetric disc the per-particle forces cancel in the sum, but any asymmetry — above all an m=1 lopsided mode, or a tidal tail dragging the halo centroid off the nucleus — leaves a net force that accelerates the whole system off the origin.

The simulation removes this by subtracting the net halo force (divided by total mass) uniformly from every particle each step (`RemoveHaloMonopole`, enabled by default). Because the subtraction is uniform, every difference `a_i - a_j` is unchanged, so the M51a/M51b relative orbit derived in section 4 is preserved exactly — only the spurious bulk drift is cancelled. The correction is global, not per-system; a per-system correction would cancel the real mutual attraction between the two galaxies.

Note this is a fix for the *drift* symptom, not for the underlying halo-center definition. The halo center remains the mass-weighted centroid of all of that system's particles, so a long tidal tail still pulls it off the nucleus and distorts the halo field the disc feels (section 9.3). Momentum is now conserved regardless.

---

## 11. Data Sources

| Source | Used For |
|--------|----------|
| Salo & Laurikainen 2000, MNRAS 319, 377 | Orbital parameters, pericenter distance, mass ratio, M_halo/M_disc ~ 2 |
| Dobbs et al. 2010, MNRAS 403, 625 | Spiral arm evolution timeline, M_halo/M_disc ~ 2.5, gas response |
| Querejeta et al. 2015, ApJS 219, 5 | M51a photometric mass decomposition (Spitzer 3.6um) |
| Mentuch Cooper et al. 2012, ApJ 755, 165 | Stellar masses for both galaxies |
| McQuinn et al. 2016, ApJ 826, 21 | TRGB distance (8.58 +/- 0.10 Mpc) |
| Sofue et al. 1999, ApJ 523, 136 | M51a rotation curve |
| Meidt et al. 2013, ApJ 779, 45 | M51a mass distribution from CO kinematics |
| Schuster et al. 2007, A&A 461, 143 | M51 gas content (HI + H2) |
| Toomre & Toomre 1972, ApJ 178, 623 | Original tidal tail theory |
| Toomre 1964, ApJ 139, 1217 | Disc stability criterion (Q parameter) |
| de Vaucouleurs et al. 1991 (RC3) | Galaxy classifications and sizes |

---

## 12. Running the Simulation

Load `M51.sim` in the simulator. The camera is positioned along +y (face-on to M51a's disc). Key viewing notes:

- At t=0: Two galaxies visible. M51a (large, center) with M51b (compact, to the right at +x).
- Run forward to t ~ 5-6: Watch for spiral arm development in M51a.
- At t ~ 6 (pericenter): Maximum tidal interaction. M51b closest to M51a.
- At t ~ 7-9: Best M51-like morphology — two strong arms with tidal bridge.

For the best visual match to real M51, look at the system around t = 7-9 (shortly after first pericenter passage). The arm pattern should show:
- Two trailing arms opening clockwise (viewed from +y)
- A bridge of particles connecting toward M51b's position
- An opposing tidal tail on the far side of M51a
