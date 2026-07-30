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
| Pericenter distance | 15-25 kpc | Salo & Laurikainen 2000 |
| Time since first passage | 300-400 Myr | Salo & Laurikainen 2000 |
| M51b position | Behind M51a disc (farther from Earth) | Dobbs et al. 2010 |
| Orbital inclination to disc | 10-20 degrees | Salo & Laurikainen 2000 |
| Mass ratio (dynamical) | 1:3 to 1:5 | Salo & Laurikainen 2000 |
| Distance to system | 8.0-8.6 Mpc | McQuinn et al. 2016 |

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
V_total² = V_baryon² + V_halo²
```

The key constraint from Salo & Laurikainen (2000) and Dobbs et al. (2010) is that the dark matter halo encloses approximately **twice** the baryonic mass within the disc radius (M_halo/M_baryon ≈ 2). This means:

```
V_halo² / V_total² ≈ 2/3    (halo provides 2/3 of centripetal support)
V_baryon² / V_total² ≈ 1/3   (baryons provide 1/3)

Therefore:
V_baryon = V_total / sqrt(3)
V_halo   = V_total × sqrt(2/3)
```

### 2.2 Baryonic Mass vs Photometric Stellar Mass

The baryonic masses used here are **rotation-curve-decomposed** values — the dynamically cold disc mass that participates in spiral arm formation. These are deliberately lower than photometric stellar masses (e.g., Querejeta et al. 2015 gives M_stellar ~ 5×10¹⁰ M☉ for M51a):

- Photometric mass includes the dynamically hot thick disc, which acts more like a bulge/halo (doesn't respond to tidal perturbation as a cold thin disc)
- The "disc mass" for tidal interaction modeling is the cold, thin, responsive component
- This is standard practice: Salo & Laurikainen (2000) and Dobbs et al. (2010) both use reduced disc masses relative to photometric estimates

### 2.3 Halo Velocity Correction

The simulator uses a **cored isothermal** halo:
```
V_halo(r) = haloVc × r / sqrt(r² + Rc²)
```

At the disc edge (r = R), the circular velocity of the halo is:
```
V_halo(R) = haloVc × R / sqrt(R² + Rc²)
```

We want V_halo(R) to equal V_total × sqrt(2/3), so:
```
haloVc = V_total × sqrt(2/3) / (R / sqrt(R² + Rc²))
       = V_total × sqrt(2/3) × sqrt(1 + (Rc/R)²)
```

For small Rc/R (which is the case here — Rc << R), haloVc ≈ V_total × sqrt(2/3).

---

## 3. Galaxy Parameters

### 3.1 NGC 5194 (M51a)

| Physical Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | Sbc | de Vaucouleurs et al. | — |
| Stellar disc radius (R25) | 11.2 kpc | NED | 186.7 |
| Observed flat V_rotation | 210 km/s | Sofue et al. 1999; Meidt et al. 2013 | — |
| V_baryon (at disc edge) | 121 km/s | = 210/sqrt(3) | — |
| V_halo (at disc edge) | 172 km/s | = 210×sqrt(2/3) | — |
| Total baryonic mass | 2.74 × 10¹⁰ M☉ | = V_b² × R / G | 2,744,000 |
| Bulge mass (M) | 4.7 × 10⁹ M☉ | 17% of baryonic (Sbc) | 470,000 |
| Disc mass | 2.27 × 10¹⁰ M☉ | baryonic − bulge | 2,274,000 |
| Mfrac (disc/bulge) | 4.84 | = disc/bulge | — |
| DM halo Vc (haloVc) | 172.1 km/s | corrected for core | — |
| DM halo core radius (haloRc) | 1.0 kpc | Salo: ~0.5-1 kpc | 16.7 |
| M_halo within disc | 5.49 × 10¹⁰ M☉ | cored isothermal | 5,488,000 |
| M_halo / M_baryon | 2.00 | Target: ~2 | — |
| Inner hole (bulge region) | 0.3 kpc | — | 5.0 |

### 3.2 NGC 5195 (M51b)

| Physical Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | SB0-pec | de Vaucouleurs et al. | — |
| Optical radius | 2.1 kpc | NED | 35.0 |
| Observed V_rotation | 130 km/s | mass-Vc scaling | — |
| V_baryon (at disc edge) | 75.1 km/s | = 130/sqrt(3) | — |
| V_halo (at disc edge) | 106.1 km/s | = 130×sqrt(2/3) | — |
| Total baryonic mass | 1.97 × 10⁹ M☉ | = V_b² × R / G | 197,167 |
| Bulge mass (M) | 8.0 × 10⁸ M☉ | 40% of baryonic (SB0) | 80,000 |
| Disc mass | 1.17 × 10⁹ M☉ | baryonic − bulge | 117,167 |
| Mfrac (disc/bulge) | 1.46 | = disc/bulge | — |
| DM halo Vc (haloVc) | 110.4 km/s | corrected for core | — |
| DM halo core radius (haloRc) | 0.6 kpc | compact galaxy | 10.0 |
| M_halo within disc | 3.94 × 10⁹ M☉ | cored isothermal | 394,333 |
| M_halo / M_baryon | 2.00 | Target: ~2 | — |
| Inner hole | 0.2 kpc | — | 3.3 |

### 3.3 Mass Ratio

The dynamical mass ratio at the encounter distance determines the tidal interaction strength. At the pericenter distance (15 kpc = 250 code units), using the cored isothermal formula M_halo(r) = Vc² × r³ / (r² + Rc²):

```
M51a at 15 kpc:
  baryonic:     2,744,000
  halo:         172.1² × 250³ / (250² + 16.7²) = 7,376,000
  total:        10,120,000

M51b at 15 kpc:
  baryonic:     197,000
  halo:         110.4² × 250³ / (250² + 10.0²) = 3,042,000
  total:        3,239,000

Dynamical ratio: 3,239,000 / 10,120,000 = 1 : 3.1
```

This is within the observational range (1:3 to 1:5 from Salo & Laurikainen 2000).

---

## 4. Orbital Computation

### 4.1 Strategy

We place M51b at **apocenter** (the farthest point of its bound orbit) at t=0, with all velocity tangential (v_r = 0). The simulation then runs forward as M51b falls inward toward pericenter, exciting spiral arms in M51a.

This approach has two advantages:
1. Starting at apocenter means the initial velocity is purely tangential — no radial component to specify.
2. The spiral arm excitation happens gradually during infall, producing a realistic build-up of arm strength.

### 4.2 Effective Mass

In a system with dark matter halos, the gravitating mass depends on the separation. For the **cored isothermal** halo, the enclosed mass within radius r is:

```
M_halo(r) = Vc² × r³ / (r² + Rc²)     [G=1]
```

The effective mass determining the orbit at separation r is the sum of both halos' enclosed masses plus all baryonic mass:

```
M_eff(r) = Vc_a² × r³/(r² + Rc_a²) + 0.5 × Vc_b² × r³/(r² + Rc_b²) + M_baryon_a + M_baryon_b
```

(The factor 0.5 on M51b's halo accounts for the fact that at large separations, only the inner portion of M51b's less extended halo contributes.)

At the initial separation (667 code units = 40 kpc):
```
M_eff = 172.1² × 667³/(667² + 16.7²) + 0.5 × 110.4² × 667³/(667² + 10.0²) + 2,744,000 + 197,000
      ≈ 19,745,000 + 4,065,000 + 2,941,000
      ≈ 26,700,000 code mass units
```

At pericenter (250 code units = 15 kpc):
```
M_eff = 172.1² × 250³/(250² + 16.7²) + 0.5 × 110.4² × 250³/(250² + 10.0²) + 2,941,000
      ≈ 7,376,000 + 1,521,000 + 2,941,000
      ≈ 11,838,000 code mass units
```

### 4.3 Energy and Angular Momentum Conservation

At apocenter (r = r_start, v_r = 0):
```
E = (1/2) × v_t² - M_eff(r_start) / r_start
L = r_start × v_t
```

At pericenter (r = r_peri, all velocity is tangential):
```
E = (1/2) × v_p² - M_eff(r_peri) / r_peri
L = r_peri × v_p
```

From angular momentum conservation:
```
v_p = r_start × v_t / r_peri
```

Substituting into the energy equation and solving for v_t:
```
v_t² = 2 × (M_eff(r_start)/r_start - M_eff(r_peri)/r_peri) / (1 - (r_start/r_peri)²)
```

Both numerator and denominator are negative (since r_start > r_peri and M_eff/r increases toward smaller r for halo-dominated potentials), so v_t² is positive.

### 4.4 Computed Orbital Parameters

| Parameter | Value |
|-----------|-------|
| Initial separation (apocenter) | 667 code units (40 kpc) |
| Pericenter distance | 250 code units (15 kpc) |
| Tangential velocity at apocenter | 48.6 km/s |
| Velocity at pericenter | 130 km/s |
| Specific orbital energy | −38,936 (bound) |
| Half-orbit time (to pericenter) | ~8 code time units (~490 Myr) |
| Full orbital period | ~17 code time units (~980 Myr) |

The ~490 Myr half-orbit time is slightly longer than Salo & Laurikainen's (2000) 300-400 Myr, consistent with the fact that our analytic halos are not truncated and extend to infinity — the effective mass at large radius is slightly overestimated, producing a wider orbit. The first-passage morphology is insensitive to this (it depends on pericenter distance and encounter velocity, not the approach timescale).

### 4.5 Coordinate System and Orbital Inclination

- M51a disc lies in the **x-z plane** with disc normal along **+y**
- M51a rotates counter-clockwise when viewed from +y
- M51b starts along the **+x axis** at 644 code units (40 kpc × cos 15°)
- Prograde orbit: M51b's tangential velocity is in the **+z direction** (matching disc rotation sense)
- Orbital inclination: 15° tilt from the disc plane, giving a small **+y velocity** component (12.6 km/s)

The inclination ensures M51b passes slightly above/below the disc plane rather than exactly through it — consistent with observations showing M51b is currently behind M51a's disc.

---

## 5. Dark Matter Halos

### 5.1 Role in the Interaction

The dark matter halos play one critical role here, and conspicuously fail to play a second one.

1. **Tidal mass** (modeled): The effective gravitating mass at the interaction distance is dominated by the halos (Vc² × r³/(r²+Rc²) >> M_baryonic for r > 100 code units). This sets the orbital velocity and the encounter timescale, and the cross-halo term in the code carries it correctly.

2. **Dynamical friction** (NOT modeled): In reality, as M51b moves through M51a's extended halo it raises a trailing density wake and loses orbital energy to gravitational drag (Chandrasekhar friction). **The simulation does not reproduce this.** The analytic halos are rigid spherical potentials comoving with their own galaxy's barycenter, and such a potential is symmetric by construction — it exerts no net drag. The cross-halo term is a conservative central force and does no secular work on the orbit. The only friction present comes from the live baryonic particles, which are a minority of the mass.

   Consequence for this script: orbital decay is badly underestimated, so M51b will not sink on a realistic timescale. The *first passage* — which is what this simulation is set up to reproduce — is unaffected, since friction has had no time to act. Do not trust the post-encounter orbit, the return time, or any eventual merger.

### 5.2 Halo Parameters

The simulator uses a **cored isothermal** halo profile:
```
a_halo = Vc² × r / (r² + Rc²)
```

This gives:
- At r >> Rc: flat rotation curve with V = Vc
- At r << Rc: linear rotation (solid body), V = Vc × r / Rc
- Core radius Rc prevents singular density at the center

| Galaxy | haloVc (km/s) | haloRc (code units) | haloRc (kpc) | M_halo/M_baryon |
|--------|---------------|---------------------|--------------|-----------------|
| M51a | 172.1 | 16.7 | 1.0 | 2.00 |
| M51b | 110.4 | 10.0 | 0.6 | 2.00 |

Note: haloVc is NOT the total observed rotation velocity. It is the halo-only contribution, corrected for the core radius so that V_halo(R_disc) = V_total × sqrt(2/3). The total rotation curve results from the quadrature sum of baryonic and halo contributions.

---

## 6. Particle Distribution

### 6.1 Counts and Resolution

| Galaxy | Particles | Purpose |
|--------|-----------|---------|
| M51a | 160,000 (159,999 disc + 1 central) | High resolution for spiral arm structure |
| M51b | 40,000 (39,999 disc + 1 central) | Adequate for compact companion |
| Total | 200,000 | |

M51a gets 4× more particles because:
- It is the galaxy where spiral arms must be resolved
- Its disc is ~5× larger in radius (much larger area to sample)
- Spiral arm contrast requires adequate particle density to be visible

### 6.2 Softening and Particle Count

The gravitational softening length should scale as 1/sqrt(N) to maintain consistent two-body relaxation:

| N_M51a | r_soft | Motivation |
|--------|--------|------------|
| 40,000 | 0.3 | Low resolution |
| 100,000 | 0.2 | Medium resolution |
| 160,000 | 0.15 | High resolution |

The M51.sim script uses r_soft = 0.3, appropriate for N=40,000. When increasing particle count, reduce r_soft accordingly.

### 6.3 Disc Generation

The `GalaxyDisc` generator creates particles in a disc with:
- Surface density ∝ 1/r (from radius Ri to R)
- Circular orbital velocity from enclosed mass + halo contribution:
  `v_c² = G × M_enclosed/r + haloVc² × r²/(r² + haloRc²)`
- Small velocity tolerance (Vtol = 0.1 = 10% scatter around circular velocity)
- Disc orientation set by the normal vector

For M51b, the disc normal is tilted 15° from +y to align roughly with its orbital plane:
```
normal = (sin(15°), cos(15°), 0) = (0.2588, 0.9659, 0.0)
```

---

## 7. Expected Simulation Behavior

### 7.1 Timeline

| Phase | Sim time | Physical time | What happens |
|-------|----------|---------------|--------------|
| Approach | t = 0-4 | 0-230 Myr | M51b falls inward. Tidal perturbation of M51a's outer disc begins. Weak arm-like features start to develop. |
| Pre-pericenter | t = 4-8 | 230-470 Myr | Spiral arms grow stronger as M51b accelerates inward. |
| Pericenter | t ~ 8 | ~490 Myr | Closest approach (15 kpc). Strong tidal torque on disc material. Two-arm spiral pattern is strongly excited. |
| Recession | t = 8-13 | 490-760 Myr | M51b recedes. Grand-design spiral arms are fully developed. Tidal bridge/tail connects the galaxies. |
| Second passage | t ~ 17 | ~980 Myr | M51b returns. Arms may be reinforced or disrupted depending on phase alignment. |

### 7.2 What to Look For

1. **Two-arm spiral arms**: Coherent, symmetric spiral arms extending from M51a's inner disc to its outer edge. These should be visible from t ~ 6 onward.

2. **Tidal bridge**: A stream of particles connecting M51a and M51b, pulled out of M51a's disc by M51b's gravity during closest approach.

3. **Tidal tail**: An arm extending in the opposite direction from the bridge (anti-companion side), formed by angular momentum transfer.

4. **M51b compression**: The compact companion may develop a slightly distorted morphology as it passes through M51a's tidal field.

5. **Arm winding**: After pericenter, the inner parts of the spiral arms rotate faster than the outer parts (differential rotation), causing progressive winding. The pattern should remain open for ~2-3 code time units after pericenter before winding significantly.

### 7.3 Why the Arms Form

The physics of tidal spiral arm formation:

1. M51b's gravitational field creates a **tidal quadrupole** on M51a's disc — stretching it along the line connecting the two galaxies and compressing it perpendicular to that line.

2. Disc material on the **near side** of M51a (facing M51b) experiences stronger attraction toward M51b than the disc center does. It is pulled outward.

3. Disc material on the **far side** experiences weaker attraction than the center. Relative to the center, it appears to be pushed outward (tidal effect).

4. Because the encounter is **prograde**, this quadrupole rotates at roughly the same rate as disc material orbits. The perturbation is coherent over multiple orbital periods, building up a strong two-arm response.

5. The **spiral** shape arises because the inner disc orbits faster than the outer disc (differential rotation). A linear perturbation at different radii gets sheared into a trailing spiral by differential rotation.

---

## 8. Comparison to Published Models

### 8.1 Salo & Laurikainen 2000

Their best-fit model parameters:
- Pericenter: 20-25 kpc
- Mass ratio: 1:3
- M_halo/M_disc: ~2
- Orbital inclination: 10-20°
- Time since first passage: 350-400 Myr
- Prograde encounter

Our parameters: pericenter 15 kpc, dynamical mass ratio 1:3.1 (at pericenter), M_halo/M_baryon = 2.0, inclination 15°, time-to-pericenter ~490 Myr. Good agreement on mass decomposition and encounter geometry; our slightly closer pericenter and longer approach time reflect the choice to use a bound orbit starting from 40 kpc apocenter rather than a parabolic flyby.

### 8.2 Dobbs et al. 2010

Their SPH + N-body model adds gas physics (ISM, star formation) on top of the gravitational interaction. They find:
- M_halo/M_disc: ~2.46 (lowered Evans model)
- Spiral arms visible ~200 Myr before pericenter
- Strongest arm-interarm contrast at ~100-200 Myr after pericenter
- Pericenter 25 kpc, mass ratio 1:3

Our collisionless (gravity-only) simulation should reproduce the stellar spiral arm morphology but not gas features (HII regions, dust lanes). Our M_halo/M_baryon = 2.0 is consistent with both papers' range of 2–2.5.

### 8.3 Limitations of This Simulation

1. **No gas physics**: Real M51 has ~25% of its disc mass in gas, which responds more strongly to tidal compression (gas shocks, star formation in arms). Our disc is collisionless.

2. **Simplified dark matter halos**: We use analytic cored-isothermal halos rather than live DM particle halos. The halos are rigid spherical backgrounds whose centers track their own galaxy's particle barycenter. This means no self-consistent halo response, no tidal stripping (M51b keeps its halo Vc even after plunging through M51a's disc), and no dynamical friction (§5.1). The halos also have infinite extent — `M_halo(r) ~ Vc²r` never stops growing — so long-range attraction is overestimated relative to a truncated halo. Additionally, because each halo center is the mass-weighted centroid of *all* its member particles and system membership is never reassigned, a long tidal tail drags the halo center off the nucleus, and stars stripped by the companion keep pulling their original halo toward it. This is a first-passage simulation, which is the regime where these limitations are least damaging. See `docs/dark-matter-halo.md` for the full accounting.

3. **No second passage history**: The real M51 may have undergone multiple passages. We simulate from the pre-first-encounter state.

4. **Disc thickness**: The particle generator creates thin discs. Real galaxy discs have finite scale heights (especially after tidal heating). The 2D projection will look correct but the 3D structure is idealized.

---

## 9. Simulation Parameters

### 9.1 Time Step and Softening

- **dt = 0.0005**: Required for the fast orbital velocities near M51a's center (v_circular ~ 210 km/s at ~83 code units). At r=83, orbital period ~ 2π×83/210 ~ 2.5 code time units, so dt/T ~ 0.0002 (well-resolved).

- **r_soft = 0.3**: Gravitational softening (18 pc). Prevents close-encounter singularities. The minimum physical scale resolved is ~2×r_soft = 36 pc.

- **BH_Opening_Theta = 0.5**: Barnes-Hut opening angle. A value of 0.5 gives good accuracy for the tidal interaction (force errors < 1%).

### 9.2 PinCentralBodies

The simulation uses `PinCentralBodies`, which pins each galaxy's central body to its system's center of mass (both position and velocity). It is applied every step to every system with `halo_vc > 0`. It prevents:
- The central body from wandering away from the galaxy due to N-body noise
- Unphysical recoil from discrete particle encounters

The central body still feels and transmits the tidal field from the other galaxy.

### 9.3 Cross-Halo Gravity

The code includes cross-system halo gravity (the cross-halo loop in `CalcAccelRangeOct` and `CalcAccelRangeP2P`). Every body feels its own system's halo plus the halo of every other system. This means:
- M51a particles feel M51b's dark matter halo potential
- M51b particles feel M51a's dark matter halo potential
- This creates the correct tidal field for spiral arm excitation

Each halo is centered on its own system's mass-weighted particle barycenter, recomputed every derivative evaluation by `ComputeHaloCenters()` — not anchored to the central body. Because a halo's field is nearly uniform across the other galaxy at these separations, and because each halo tracks its own galaxy, the pair's relative acceleration comes out correct at `G(M_totA + M_totB)/d²` including halo mass. This is why the analytic orbit derived in §4 is reproduced by the simulation.

What the cross-halo term does **not** do is produce dynamical friction — see §5.1. It is a conservative central force between two rigid symmetric potentials.

---

## 10. Data Sources

| Source | Used For |
|--------|----------|
| Salo & Laurikainen 2000, MNRAS 319, 377 | Orbital parameters, pericenter distance, mass ratio, M_halo/M_disc ≈ 2 |
| Dobbs et al. 2010, MNRAS 403, 625 | Spiral arm evolution timeline, M_halo/M_disc ≈ 2.5, gas response |
| Querejeta et al. 2015, ApJS 219, 5 | M51a photometric mass decomposition (Spitzer 3.6μm) |
| Mentuch Cooper et al. 2012, ApJ 755, 165 | Stellar masses for both galaxies |
| McQuinn et al. 2016, ApJ 826, 21 | TRGB distance (8.58 ± 0.10 Mpc) |
| Sofue et al. 1999, ApJ 523, 136 | M51a rotation curve |
| Meidt et al. 2013, ApJ 779, 45 | M51a mass distribution from CO kinematics |
| Schuster et al. 2007, A&A 461, 143 | M51 gas content (HI + H₂) |
| Toomre & Toomre 1972, ApJ 178, 623 | Original tidal tail theory |
| de Vaucouleurs et al. 1991 (RC3) | Galaxy classifications and sizes |

---

## 11. Running the Simulation

Load `M51.sim` in the simulator. The camera is positioned along +y (face-on to M51a's disc). Key viewing notes:

- At t=0: Two galaxies visible. M51a (large, center) with M51b (compact, to the right at +x).
- Run forward to t ~ 6-8: Watch for spiral arm development in M51a.
- At t ~ 8 (pericenter): Maximum tidal interaction. M51b closest to M51a.
- At t ~ 9-10: Best M51-like morphology — two strong arms with tidal bridge.

For the best visual match to real M51, look at the system around t = 9-10 (shortly after first pericenter passage). The arm pattern should show:
- Two trailing arms opening counter-clockwise (viewed from +y)
- A bridge of particles connecting toward M51b's position
- An opposing tidal tail on the far side of M51a
