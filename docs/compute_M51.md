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
| Mass ratio (dynamical) | 1:3 to 1:5 | Querejeta et al. 2015 |
| Distance to system | 8.0-8.6 Mpc | McQuinn et al. 2016 |

### 1.3 Why Prograde?

The prograde nature of the encounter (companion orbiting in the same sense as disc rotation) is critical for producing strong spiral arms:

- **Prograde encounter**: The companion's tidal field rotates at a rate comparable to disc material's orbital frequency. This produces a near-resonant perturbation that lasts many orbital periods, creating strong, coherent two-arm spirals.
- **Retrograde encounter**: The tidal field sweeps past disc material rapidly (in the counter-rotating sense). The perturbation averages out over less than one orbit, producing only weak, transient disturbances.

This is analogous to pushing a swing: pushing in phase with the motion (prograde) builds amplitude; pushing against the motion (retrograde) does little.

---

## 2. Galaxy Parameters

### 2.1 NGC 5194 (M51a)

| Physical Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | Sbc | de Vaucouleurs et al. | — |
| Stellar disc radius (R25) | 11.2 kpc | NED | 186.7 |
| Bulge mass | 1.0 × 10¹⁰ M☉ | Querejeta et al. 2015 | 1,000,000 |
| Disc stellar mass | 4.0 × 10¹⁰ M☉ | Mentuch Cooper et al. 2012 | 4,000,000 |
| Gas mass (HI + H₂) | 0.9 × 10¹⁰ M☉ | Schuster et al. 2007 | 900,000 |
| Flat rotation velocity | 210 km/s | Sofue et al. 1999; Meidt et al. 2013 | 210.0 |
| DM halo core radius | ~5 kpc | Adopted from rotation curve fitting | 83.3 |
| Inner hole (bulge region) | 0.3 kpc | — | 5.0 |

The disc mass fraction (Mfrac) represents the ratio of disc+gas mass to the central concentration:

```
Mfrac = (disc_mass + gas_mass) / bulge_mass = (4.0e10 + 0.9e10) / 1.0e10 = 4.90
```

This high Mfrac (compared to ~0.5 for a typical compact galaxy) means M51a is strongly disc-dominated, consistent with its Sbc classification and active star formation.

### 2.2 NGC 5195 (M51b)

| Physical Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | SB0-pec | de Vaucouleurs et al. | — |
| Optical radius | 2.1 kpc | NED | 35.0 |
| Baryonic mass (total) | 1.0 × 10¹⁰ M☉ | Mentuch Cooper et al. 2012 | 1,000,000 |
| Disc component | 5.0 × 10⁹ M☉ | Estimated | 500,000 |
| Rotation velocity | ~130 km/s | From mass-Vc scaling | 130.0 |
| DM halo core radius | ~2.5 kpc | Adopted | 41.7 |
| Inner hole | 0.2 kpc | — | 3.3 |

M51b is compact and bulge-dominated (Mfrac = 0.50). Its disturbed morphology ("pec" in its classification) is a direct result of the interaction — tidal forces have stripped and heated its outer disc material.

### 2.3 Mass Ratio

The baryonic mass ratio is:

```
M51b_baryonic / M51a_baryonic = 1.5e10 / 5.9e10 = 1 : 3.9
```

This is within the observational range (1:3 to 1:5). The effective dynamical mass ratio (including dark matter halos at the interaction distance) is slightly higher since M51a's more massive halo dominates.

---

## 3. Orbital Computation

### 3.1 Strategy

We place M51b at **apocenter** (the farthest point of its bound orbit) at t=0, with all velocity tangential (v_r = 0). The simulation then runs forward as M51b falls inward toward pericenter, exciting spiral arms in M51a.

This approach has two advantages:
1. Starting at apocenter means the initial velocity is purely tangential — no radial component to specify.
2. The spiral arm excitation happens gradually during infall, producing a realistic build-up of arm strength.

### 3.2 Effective Mass

In a system with dark matter halos, the gravitating mass is not constant but depends on the separation. For an isothermal halo, the enclosed mass within radius r is:

```
M_halo(r) = Vc² × r     [G=1]
```

The effective mass determining the orbit at separation r is the sum of both halos' enclosed masses plus all baryonic mass:

```
M_eff(r) = Vc_a² × r + 0.5 × Vc_b² × r + M_baryonic_a + M_baryonic_b
```

(The factor 0.5 on M51b's halo accounts for the fact that at large separations, only the inner portion of M51b's less extended halo contributes.)

At the initial separation (667 code units = 40 kpc):
```
M_eff = 210² × 667 + 0.5 × 130² × 667 + 5,900,000 + 1,500,000
      = 29,400,000 + 5,630,000 + 7,400,000
      = 42,400,000 code mass units
```

At pericenter (333 code units = 20 kpc):
```
M_eff = 210² × 333 + 0.5 × 130² × 333 + 7,400,000
      = 14,700,000 + 2,820,000 + 7,400,000
      = 24,900,000 code mass units
```

### 3.3 Energy and Angular Momentum Conservation

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

### 3.4 Computed Orbital Parameters

| Parameter | Value |
|-----------|-------|
| Initial separation (apocenter) | 667 code units (40 kpc) |
| Pericenter distance | 333 code units (20 kpc) |
| Tangential velocity at apocenter | 86 km/s |
| Velocity at pericenter | 172 km/s |
| Specific orbital energy | -59,950 (bound) |
| Half-orbit time (to pericenter) | ~7.5 code time units (~440 Myr) |
| Full orbital period | ~15 code time units (~880 Myr) |

The ~440 Myr half-orbit time is consistent with Salo & Laurikainen (2000), who find the first passage occurred 300-400 Myr ago with M51b on a moderately eccentric orbit.

### 3.5 Coordinate System and Orbital Inclination

- M51a disc lies in the **x-z plane** with disc normal along **+y**
- M51a rotates counter-clockwise when viewed from +y
- M51b starts along the **+x axis** at 644 code units (40 kpc × cos 15°)
- Prograde orbit: M51b's tangential velocity is in the **+z direction** (matching disc rotation sense)
- Orbital inclination: 15° tilt from the disc plane, giving a small **+y velocity** component (22.3 km/s)

The inclination ensures M51b passes slightly above/below the disc plane rather than exactly through it — consistent with observations showing M51b is currently behind M51a's disc.

---

## 4. Dark Matter Halos

### 4.1 Role in the Interaction

The dark matter halos play two critical roles:

1. **Dynamical friction**: As M51b moves through M51a's extended halo, it loses orbital energy via gravitational drag (Chandrasekhar dynamical friction). This is modeled implicitly through the particle interactions and the analytic halo potential (cross-halo gravitational term in the code).

2. **Tidal mass**: The effective gravitating mass at the interaction distance is dominated by the halos (Vc² × r >> M_baryonic for r > 100 code units). This determines the orbital velocity and timescale.

### 4.2 Halo Parameters

The simulator uses a **cored isothermal** halo profile:
```
a_halo = Vc² × r / (r² + Rc²)
```

This gives:
- At r >> Rc: flat rotation curve with V = Vc
- At r << Rc: linear rotation (solid body), V = Vc × r / Rc
- Core radius Rc prevents singular density at the center

| Galaxy | Vc (km/s) | Rc (code units) | Rc (kpc) |
|--------|-----------|-----------------|----------|
| M51a | 210 | 83.3 | 5.0 |
| M51b | 130 | 41.7 | 2.5 |

---

## 5. Particle Distribution

### 5.1 Counts and Resolution

| Galaxy | Particles | Purpose |
|--------|-----------|---------|
| M51a | 40,000 (39,999 disc + 1 central) | High resolution for spiral arm structure |
| M51b | 10,000 (9,999 disc + 1 central) | Adequate for compact companion |
| Total | 50,000 | |

M51a gets 4× more particles because:
- It is the galaxy where spiral arms must be resolved
- Its disc is ~5× larger in radius (much larger area to sample)
- Spiral arm contrast requires adequate particle density to be visible

### 5.2 Disc Generation

The `GalaxyDisc` generator creates particles in a disc with:
- Surface density ∝ 1/r (from radius Ri to R)
- Circular orbital velocity from enclosed mass + halo contribution
- Small velocity tolerance (Vtol = 0.1 = 10% scatter around circular velocity)
- Disc orientation set by the normal vector

For M51b, the disc normal is tilted 15° from +y to align roughly with its orbital plane:
```
normal = (sin(15°), cos(15°), 0) = (0.2588, 0.9659, 0.0)
```

---

## 6. Expected Simulation Behavior

### 6.1 Timeline

| Phase | Sim time | Physical time | What happens |
|-------|----------|---------------|--------------|
| Approach | t = 0-4 | 0-230 Myr | M51b falls inward. Tidal perturbation of M51a's outer disc begins. Weak arm-like features start to develop. |
| Pericenter | t ~ 7 | ~440 Myr | Closest approach (20 kpc). Strong tidal torque on disc material. Two-arm spiral pattern is strongly excited. |
| Recession | t = 7-11 | 440-650 Myr | M51b recedes. Grand-design spiral arms are fully developed. Tidal bridge/tail connects the galaxies. |
| Second passage | t ~ 15 | ~880 Myr | M51b returns. Arms may be reinforced or disrupted depending on phase alignment. |

### 6.2 What to Look For

1. **Two-arm spiral arms**: Coherent, symmetric spiral arms extending from M51a's inner disc to its outer edge. These should be visible from t ~ 5 onward.

2. **Tidal bridge**: A stream of particles connecting M51a and M51b, pulled out of M51a's disc by M51b's gravity during closest approach.

3. **Tidal tail**: An arm extending in the opposite direction from the bridge (anti-companion side), formed by angular momentum transfer.

4. **M51b compression**: The compact companion may develop a slightly distorted morphology as it passes through M51a's tidal field.

5. **Arm winding**: After pericenter, the inner parts of the spiral arms rotate faster than the outer parts (differential rotation), causing progressive winding. The pattern should remain open for ~2-3 code time units after pericenter before winding significantly.

### 6.3 Why the Arms Form

The physics of tidal spiral arm formation:

1. M51b's gravitational field creates a **tidal quadrupole** on M51a's disc — stretching it along the line connecting the two galaxies and compressing it perpendicular to that line.

2. Disc material on the **near side** of M51a (facing M51b) experiences stronger attraction toward M51b than the disc center does. It is pulled outward.

3. Disc material on the **far side** experiences weaker attraction than the center. Relative to the center, it appears to be pushed outward (tidal effect).

4. Because the encounter is **prograde**, this quadrupole rotates at roughly the same rate as disc material orbits. The perturbation is coherent over multiple orbital periods, building up a strong two-arm response.

5. The **spiral** shape arises because the inner disc orbits faster than the outer disc (differential rotation). A linear perturbation at different radii gets sheared into a trailing spiral by differential rotation.

---

## 7. Comparison to Published Models

### 7.1 Salo & Laurikainen 2000

Their best-fit model parameters:
- Pericenter: 20-25 kpc
- Mass ratio: 1:3
- Orbital inclination: 10-20°
- Time since first passage: 350-400 Myr
- Prograde encounter

Our parameters: pericenter 20 kpc, mass ratio 1:3.9, inclination 15°, time-to-pericenter 440 Myr. Excellent agreement.

### 7.2 Dobbs et al. 2010

Their SPH + N-body model adds gas physics (ISM, star formation) on top of the gravitational interaction. They find:
- Spiral arms visible ~200 Myr before pericenter
- Strongest arm-interarm contrast at ~100-200 Myr after pericenter
- Pericenter 25 kpc, mass ratio 1:3

Our collisionless (gravity-only) simulation should reproduce the stellar spiral arm morphology but not gas features (HII regions, dust lanes).

### 7.3 Limitations of This Simulation

1. **No gas physics**: Real M51 has ~25% of its disc mass in gas, which responds more strongly to tidal compression (gas shocks, star formation in arms). Our disc is collisionless.

2. **Simplified dark matter halos**: We use analytic cored-isothermal halos rather than live DM particle halos. This means no self-consistent halo response, tidal stripping, or dynamical friction from halo particles. The halos are rigid backgrounds.

3. **No second passage history**: The real M51 may have undergone multiple passages. We simulate from the pre-first-encounter state.

4. **Disc thickness**: The particle generator creates thin discs. Real galaxy discs have finite scale heights (especially after tidal heating). The 2D projection will look correct but the 3D structure is idealized.

---

## 8. Simulation Parameters

### 8.1 Time Step and Softening

- **dt = 0.0005**: Required for the fast orbital velocities near M51a's center (v_circular ~ 210 km/s at ~83 code units). At r=83, orbital period ~ 2π×83/210 ~ 2.5 code time units, so dt/T ~ 0.0002 (well-resolved).

- **r_soft = 0.1**: Gravitational softening (6 pc). Prevents close-encounter singularities. The minimum physical scale resolved is ~2×r_soft = 12 pc.

- **BH_Opening_Theta = 0.5**: Barnes-Hut opening angle. A value of 0.5 gives good accuracy for the tidal interaction (force errors < 1%).

### 8.2 PinCentralBodies

The simulation uses `PinCentralBodies` (enabled by default when `N_Systems > 1`), which pins each galaxy's central body to its system's center of mass. This prevents:
- The central body from wandering away from the galaxy due to N-body noise
- Unphysical recoil from discrete particle encounters

The central body still feels and transmits the tidal field from the other galaxy.

### 8.3 Cross-Halo Gravity

The code includes cross-system halo-halo gravitational attraction (the cross-halo loop in `CalcAccelRangeOct`). This means:
- M51a particles feel M51b's dark matter halo potential
- M51b particles feel M51a's dark matter halo potential
- This creates the correct tidal field for spiral arm excitation

---

## 9. Data Sources

| Source | Used For |
|--------|----------|
| Salo & Laurikainen 2000, A&A 359, 895 | Orbital parameters, pericenter distance, mass ratio, encounter geometry |
| Dobbs et al. 2010, MNRAS 403, 625 | Spiral arm evolution timeline, gas response comparison |
| Querejeta et al. 2015, ApJS 219, 5 | M51a mass decomposition (bulge/disc/halo from Spitzer 3.6μm) |
| Mentuch Cooper et al. 2012, ApJ 755, 165 | Stellar masses for both galaxies |
| McQuinn et al. 2016, ApJ 826, 21 | TRGB distance (8.58 ± 0.10 Mpc) |
| Sofue et al. 1999, ApJ 523, 136 | M51a rotation curve |
| Meidt et al. 2013, ApJ 779, 45 | M51a mass distribution from CO kinematics |
| Schuster et al. 2007, A&A 461, 143 | M51 gas content (HI + H₂) |
| Toomre & Toomre 1972, ApJ 178, 623 | Original tidal tail theory |
| de Vaucouleurs et al. 1991 (RC3) | Galaxy classifications and sizes |

---

## 10. Running the Simulation

Load `M51.sim` in the simulator. The camera is positioned along +y (face-on to M51a's disc). Key viewing notes:

- At t=0: Two galaxies visible. M51a (large, center) with M51b (compact, to the right at +x).
- Run forward to t ~ 5-7: Watch for spiral arm development in M51a.
- At t ~ 7 (pericenter): Maximum tidal interaction. M51b closest to M51a.
- At t ~ 8-10: Best M51-like morphology — two strong arms with tidal bridge.

For the best visual match to real M51, look at the system around t = 8-9 (shortly after first pericenter passage). The arm pattern should show:
- Two trailing arms opening counter-clockwise (viewed from +y)
- A bridge of particles connecting toward M51b's position
- An opposing tidal tail on the far side of M51a
