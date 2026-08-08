# compute_M51.py

## 1. The M51 System

The Whirlpool Galaxy (M51, NGC 5194/5195) is the archetypal grand-design spiral shaped by a tidal interaction:

- **NGC 5194 (M51a)** — a large Sbc spiral with prominent two-arm structure.
- **NGC 5195 (M51b)** — a compact SB0-pec lenticular companion that has crossed M51a's disc plane at least once, exciting those arms.

M51 is the clearest known case where spiral structure is demonstrably caused by a companion rather than internal instability: the two-arm pattern is coherent across the whole disc, the signature of external forcing rather than swing amplification.

This project reproduces the **bound multiple-encounter model** of Salo & Laurikainen (2000), MNRAS 319, 377 (PDF in `docs/`). Every parameter is taken from that paper.

### 1.1 Distance and Unit Conversion

Angular sizes are converted at the paper's **own assumed distance of 9.6 Mpc** (Scoville & Young 1983):

```
1 arcsec = 0.046542 kpc     at 9.6 Mpc
```

The modern TRGB distance is 8.58 Mpc (McQuinn et al. 2016) and is the better measurement of the real galaxy, but it is the wrong choice here. The paper's encounter distances are quoted as ratios of its own disc truncation radius `Rd`, and its `Rcross` = 1.2–1.4 range refers to `Rd` = 18.6 kpc. Mixing a modern distance with those dimensionless ratios would change the encounter geometry.

### 1.2 Observational Constraints

| Parameter | Value | Source |
|-----------|-------|--------|
| Encounter type | Bound, multiple-passage | S&L 2000 |
| Orbital inclination `iorb` | 75–85 deg (near-polar) | Table 2 |
| Eccentricity | e ~ 0.2 | Fig. 1 caption |
| Principal crossing `Rcross` | 1.2–1.4 primary Rd | Table 2 |
| Most recent crossing `Rdown` | 1.2–1.3 primary Rd | Table 2 |
| Time since principal crossing `Tobs` | 5.5–6.5 (440–520 Myr) | Table 2 |
| Time since most recent crossing | 0.5–1.0 (40–80 Myr) | Table 2 |
| Mass ratio `Mp` (total, within 4Re) | 0.5–0.7, nominal 0.55 | Tables 1–2 |
| Toomre `Q_T` | 1.5 | Section 2.2 |
| `M_disc/M_tot` within 4Re | 1/3 | Section 2.2 |

The paper parameterizes by **disc-plane crossing distance** rather than pericenter, "motivated by the high inclination of the relative orbit." Its length unit is the primary `Rd` = 400 arcsec, and its time unit is 80 Myr.

### 1.3 Why the Bound Model

The paper tests two orbit classes and finds both reproduce the gross morphology:

1. **Near-parabolic single passage** (e = 0.67–0.83) — the classical Toomre & Toomre (1972) picture.
2. **Bound multiple encounter** (e ~ 0.2, `iorb` ~ 85 deg) — the companion crosses the disc plane several times.

The bound model is preferred because it also explains kinematics the single-passage model cannot: the apparent counter-rotation of the HI tail (tilted 40–50 deg to the inner disc), gas velocities up to 700 km/s north of the companion, the S-shaped major-axis rotation curve, and the "kink" in the northern spiral arm.

The paper also reports that "any pre-existing spiral arms are washed out by the tidally triggered spiral arms," so the initial disc need not be seeded with structure.

### 1.4 Why High Inclination Matters

At `iorb` = 80 deg the orbit is near-**polar**: the companion plunges through the disc plane roughly perpendicular to it. "Prograde" versus "retrograde" is only weakly meaningful at that inclination, since there is little co-planar component to co-rotate with. The mechanism is a series of impulsive disc-plane crossings, not sustained resonant forcing.

The paper treats the high inclination as essential — it produces the out-of-plane velocities behind the S-shaped rotation curve, the 40–50 deg tail tilt, and the high peculiar velocities north of the companion. A coplanar model reproduces the arms but none of this.

---

## 2. Mass Decomposition

### 2.1 Rotation Curve Decomposition

The observed rotation velocity is split into baryonic and halo contributions:

```
V_total^2 = V_baryon^2 + V_halo^2
```

The paper fixes `M_disc/M_tot` = 1/3 within 4Re, so the halo supplies 2/3 of the centripetal support:

```
V_baryon = V_total / sqrt(3)
V_halo   = V_total * sqrt(2/3)
```

This gives `M_halo/M_baryon` = 2 within the disc radius.

### 2.2 Baryonic Mass vs Photometric Mass

These baryonic masses are **rotation-curve-decomposed** — the dynamically cold disc mass that participates in spiral arm formation. They are deliberately lower than photometric stellar masses (Querejeta et al. 2015 gives M_stellar ~ 5×10^10 M☉ for M51a), because photometric mass includes the dynamically hot thick disc, which responds to a tidal perturbation more like a bulge than a cold thin disc.

### 2.3 Halo Velocity Correction

The simulator uses a cored isothermal halo:

```
V_halo(r) = haloVc * r / sqrt(r^2 + Rc^2)
```

Requiring `V_halo(Rd)` = `V_total * sqrt(2/3)`:

```
haloVc = V_total * sqrt(2/3) * sqrt(1 + (Rc/Rd)^2)
```

Since `Rc << Rd` here, `haloVc` ≈ `V_total * sqrt(2/3)`.

---

## 3. Galaxy Parameters

### 3.1 Disc Geometry: Rd and Re

Two radii describe each disc, and they are independent:

| Symbol | Meaning |
|---|---|
| `Re` (`h_r` in code) | **Scale length** — the e-folding distance of surface density, `Sigma(r) = Sigma_0 * exp(-r/Re)`. Sets how centrally concentrated the disc is. |
| `Rd` | **Truncation radius** — where particle placement stops. Sets the disc's outer edge. |

The ratio `Rd/Re` is **not universal**. Salo & Laurikainen use:

| Galaxy | Re | Rd | Rd/Re |
|---|---|---|---|
| M51a | 100 arcsec | 400 arcsec | **4.00** |
| M51b | 33 arcsec | 240 arcsec | **7.27** |

The paper is explicit for the companion: *"take Rd = 240 arcsec ≈ 7Re for the initial companion disc, assuming Re = 33 arcsec."* Because the ratio differs between the two galaxies, `GalaxyDisc` takes the scale length as a required input independent of the truncation radius. See the `disc_scale_length` column in `docs/script-format.md`.

### 3.2 NGC 5194 (M51a)

| Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | Sbc | de Vaucouleurs et al. | — |
| Scale length Re | 100 arcsec = 4.65 kpc | S&L sect 2.2 | 77.6 |
| Truncation Rd = 4Re | 400 arcsec = 18.62 kpc | S&L sect 2.2 | 310.3 |
| Observed flat V_rotation | 210 km/s | Sofue et al. 1999; Meidt et al. 2013 | — |
| V_baryon at Rd | 121.2 km/s | = 210/sqrt(3) | — |
| V_halo at Rd | 171.5 km/s | = 210*sqrt(2/3) | — |
| Total baryonic mass | 4.56 × 10^10 M☉ | = V_b^2 * Rd / G | 4,561,123 |
| Bulge mass (M) | 7.8 × 10^9 M☉ | 17% of baryonic (Sbc) | 780,000 |
| Disc mass | 3.78 × 10^10 M☉ | baryonic − bulge | 3,781,123 |
| Mfrac (disc/bulge) | 4.85 | — | — |
| DM halo Vc | 171.5 km/s | corrected for core | — |
| DM halo core Rc | 8 arcsec = 0.372 kpc | S&L sect 2.2 | 6.2 |
| M_halo within Rd | 9.12 × 10^10 M☉ | cored isothermal | 9,122,246 |
| M_halo / M_baryon | 2.00 | target ~2 | — |
| Inner hole | 0.3 kpc | — | 5.0 |
| Toomre Q | 1.5 | S&L standard value | — |

`Rd` is the dynamical disc truncation, larger than the R25 optical radius of 11.2 kpc.

### 3.3 NGC 5195 (M51b)

The companion's mass is set by the paper's mass ratio `Mp`, not by an assumed rotation curve — `Mp` is the quantity the paper actually constrains and tests.

| Parameter | Value | Source | Code Units |
|---|---|---|---|
| Morphological type | SB0-pec | de Vaucouleurs et al. | — |
| Scale length Re | 33 arcsec = 1.54 kpc | S&L sect 2.2 | 25.6 |
| Truncation Rd = 7.3Re | 240 arcsec = 11.17 kpc | S&L sect 2.2 (70% of primary) | 186.2 |
| Mass ratio Mp | 0.55 | S&L Table 1 caption | — |
| Total mass within Rd | 7.53 × 10^10 M☉ | = Mp × M_prim | 7,525,853 |
| Total baryonic mass | 2.51 × 10^10 M☉ | = 1/3 of total | 2,508,618 |
| Bulge mass (M) | 1.00 × 10^10 M☉ | 40% of baryonic (SB0) | 1,000,000 |
| Disc mass | 1.51 × 10^10 M☉ | baryonic − bulge | 1,508,618 |
| Mfrac (disc/bulge) | 1.51 | — | — |
| DM halo Vc | 166.4 km/s | inverted from M_halo target | — |
| DM halo core Rc | 40 arcsec = 1.862 kpc | S&L sect 2.2 | 31.0 |
| M_halo within Rd | 5.02 × 10^10 M☉ | = 2/3 of total | 5,017,235 |
| M_halo / M_baryon | 2.00 | target ~2 | — |
| Implied V_rotation | 201 km/s | *derived*, not assumed | — |
| Disc tilt to M51a disc | 32.5 deg | S&L sect 2.2 (observed) | — |
| Inner hole | 0.2 kpc | — | 3.3 |
| Toomre Q | 1.5 | S&L standard value | — |

The paper notes `Mp` = 0.55 deliberately exceeds the *observed disc* mass ratio of ~0.4, because "small galaxies are likely to have more dominant halo components" (citing Persic & Salucci 1990). The 201 km/s implied rotation curve is therefore higher than M51b's observed stellar kinematics suggest — intentional in the source model.

The larger halo core (40 arcsec vs the primary's 8) reflects "the lesser degree of halo concentration expected for smaller galaxies."

**Disc tilt vs orbital inclination** are independent quantities: the 32.5 deg disc tilt comes from the companion's observed position angle and inclination, while the 80 deg orbital inclination is a property of the trajectory.

### 3.4 Mass Ratio Verification

`Mp` = M_tot(companion) / M_tot(primary) within 4Re = Rd:

```
M51a total within Rd:  4,561,123 (baryon) + 9,122,246 (halo) = 13,683,368
M51b total within Rd:  2,508,618 (baryon) + 5,017,235 (halo) =  7,525,853

Mp = 7,525,853 / 13,683,368 = 0.550      <- paper nominal, range 0.5-0.7
```

A different ratio, sometimes quoted, is the enclosed mass within the encounter separation. Because both halos are evaluated at the same radius there and the companion's is less concentrated, that ratio is closer to unity. `Mp` is the one to check against the paper.

---

## 4. Orbital Setup

### 4.1 Strategy

M51b starts at **apocenter** with purely tangential velocity (v_r = 0 by definition). The simulation then runs forward through successive disc-plane crossings.

`Rcross` is the **input**; pericenter is solved to match it. Pericenter is not independently meaningful in the paper's framework.

### 4.2 Gravitational Potential

The simulation applies cored isothermal halo accelerations:

```
a_halo = Vc^2 * r / (r^2 + Rc^2)        attractive, directed inward
```

The corresponding potential (with `a_r = -dPhi/dr`) increases outward:

```
Phi_halo(r) = +0.5 * Vc^2 * ln(r^2 + Rc^2)
```

For baryonic point masses, `Phi_baryon(r) = -G*M/r`. The mutual orbit feels both halos plus both baryonic masses:

```
Phi(r) = +0.5*Vc_a^2*ln(r^2 + Rc_a^2) + 0.5*Vc_b^2*ln(r^2 + Rc_b^2) - M_baryon/r
```

where `M_baryon` = 4,561,123 + 2,508,618 = 7,069,741.

A Keplerian `-M_eff/r` approximation is **not** valid here — the potential is a sum of logarithmic terms, and the Keplerian form underestimates the orbital velocity substantially.

### 4.3 Energy and Angular Momentum

At apocenter and pericenter both have v_r = 0, so with `L = r_apo * v_t = r_peri * v_p`:

```
v_t^2 = 2 * (Phi_peri - Phi_apo) / (1 - (r_apo/r_peri)^2)
```

Both numerator and denominator are negative (pericenter is deeper in the well; `r_apo > r_peri`), so `v_t^2` is positive.

### 4.4 Node Geometry

A **node** is a point where M51b's orbit crosses M51a's disc plane (y = 0). There are two per orbit, and the line joining them through the origin is the **line of nodes**. Gravity here is a central force, so angular momentum is conserved, the orbit plane is fixed, and the line of nodes is fixed for the whole run.

**Node geometry** is where pericenter sits relative to that line — the *argument of pericenter*. It is an independent degree of freedom: pericenter, apocenter, eccentricity and inclination can all be held fixed while the apsides slide around the orbit, changing where in the orbit the crossings occur and at what radius.

The paper's bound model states "the apocentre is between the two disc crossings," which requires:

- apocenter at **maximum height off** the disc plane
- line of nodes ~90 deg away in orbital phase → **argument of pericenter = 90 deg**
- crossings occur on the way in and out, **not** at the apsides

Construction — an orthonormal basis for the orbit plane:

```
e1 = (1, 0, 0)                       along the line of nodes, in the disc plane
e2 = (0, sin(iorb), cos(iorb))       perpendicular to e1, tilted out of the plane

position = r_apo * e2 = (0, 506.5, 89.3)    apocenter, max height off plane
velocity = -v_t * e1  = (-219.0, 0, 0)      purely tangential, in the disc plane
```

The script asserts `|pos|` equals the apocenter the velocity was solved for, and that `pos · vel` ≈ 0.

### 4.5 Computed Parameters

```
target Rcross = 1.20 Rd                <- input
  -> pericenter = 20.57 kpc            <- solved by bisection
  -> apocenter  = r_peri * (1+e)/(1-e) = 30.86 kpc
```

| Parameter | Value |
|-----------|-------|
| Target Rcross | 1.20 Rd (22.34 kpc) |
| Pericenter (solved) | 342.9 code units (20.57 kpc) |
| Apocenter (derived) | 514.3 code units (30.86 kpc) |
| Eccentricity | 0.200 |
| Orbital inclination | 80 deg |
| Argument of pericenter | 90 deg |
| Tidal strength S = Mp*(Rd/r)^3 | 0.318 |
| Tangential velocity at apocenter | 219.0 km/s |
| Velocity at pericenter | 328.4 km/s |
| v_circular at apocenter | 266.0 km/s |
| v_t / v_circ | 0.823 |
| Orbital period | 18.44 code units |

`Rcross` = 1.20 is the low end of the paper's range and therefore the strongest perturbation it allows: *"passages with Mp < 0.55 or Rcross > 1.4 seem to be too weak to account for the observations."* For reference, `Rcross` = 1.30 gives S = 0.250 and 1.40 gives S = 0.200.

### 4.6 Numerical Verification

The crossing radius is not a closed-form function of pericenter — it depends on node geometry, and the apsides precess in a logarithmic potential. The script bisects on pericenter, integrating the relative orbit each iteration (velocity Verlet, dt = 5e-5).

| Quantity | Target | Achieved |
|---|---|---|
| First crossing | 1.200 Rd | 1.200 Rd |
| Eccentricity | 0.200 | 0.200 |
| Apocenter between crossings | required | confirmed |

Event sequence:

| Event | t | Myr | r (kpc) | r / Rd |
|---|---|---|---|---|
| **close crossing #1** (`Rcross`) | 2.87 | 169 | 22.34 | **1.200** |
| pericenter | 3.73 | 219 | 20.57 | 1.105 |
| apocenter | 7.47 | 438 | 30.86 | 1.658 |
| grazing crossing | 7.49 | 440 | 30.86 | 1.657 |
| pericenter | 11.20 | 658 | 20.57 | 1.105 |
| **close crossing #2** (`Rdown`) | 12.09 | 710 | 22.44 | **1.205** |
| apocenter | 14.94 | 877 | 30.86 | 1.658 |
| **close crossing #3** | 17.84 | 1047 | 22.25 | **1.195** |

Close crossings land at 1.195–1.205 Rd, inside the paper's `Rcross` = 1.2–1.4 and `Rdown` = 1.2–1.3. This is an independent check: pericenter was the targeted quantity and the crossing radii fell in range on their own.

**Why crossings alternate close and grazing.** In a logarithmic potential the apsidal angle is not 180 deg, so the apsides **precess** relative to the fixed line of nodes. Successive crossings sample different orbital phases, alternating between ~1.20 Rd (genuine close passages) and ~1.66 Rd (near apocenter, weak). The script classifies a crossing as close when `r <= a`, the semi-major axis.

### 4.7 Best-Morphology Epoch

The paper anchors its observation epoch `Tobs` to the **principal disc-plane crossing**: "the time Tobs elapsed since the disc plane crossing at T = 0." Section 3.4.3 brackets it from both sides:

- `Tobs` ≤ 4–5 — spiral structure and the long HI tail too weak.
- `Tobs` ~ 5–6 — produces the observed **"kink" on the northern arm** toward the companion, from the interplay between the most recent crossing and arms excited by the first passage.
- `Tobs` ~ 7 — that arm "would be opened into a new bridge during the second crossing," overshooting.
- `Tobs` > 6 — the spiral becomes too tightly wound.

Applying `Tobs - Tdown` = 0.5–1.0 (40–80 Myr, i.e. 0.68–1.36 of this project's 58.7 Myr unit) after close crossing #2 at t = 12.09:

**Observe at t = 12.8 to 13.5, centre ~13.1.**

Cross-check on the longer interval: the paper wants 440–520 Myr between successive close crossings; here it is 541 Myr, about 4% above. The untruncated halo (section 9.3) overestimates enclosed mass at large radii, which affects the period for a given orbit shape.

### 4.8 Coordinate System

- M51a disc lies in the **x-z plane**, disc normal **+y**
- M51a disc rotates **clockwise** viewed from +y (`LoadGalaxyDiscState` uses `v_tan = -vm`, and `t_hat = cross(r_hat, n_hat)` is +z at the +x position)
- M51b's line of nodes is along **±x**; its apocenter is on the **+y side**

---

## 5. Dark Matter Halos

### 5.1 Role in the Interaction

The halos play one critical role and conspicuously fail to play a second.

1. **Tidal mass** (modeled): the effective gravitating mass at the interaction distance is halo-dominated, setting the orbital velocity and encounter timescale. The cross-halo term carries this correctly.

2. **Dynamical friction** (NOT modeled): in reality M51b would raise a trailing wake in M51a's halo and lose orbital energy (Chandrasekhar friction). The analytic halos here are rigid spherical potentials comoving with their own galaxy's barycenter, symmetric by construction, so they exert no net drag. The cross-halo term is a conservative central force and does no secular work. The only friction present comes from the live baryonic particles, a minority of the mass.

   Consequence: orbital decay is badly underestimated. The first passage is unaffected, since friction has had no time to act, but do not trust the long-term orbit or any eventual merger.

### 5.2 Halo Parameters

```
a_halo = Vc^2 * r / (r^2 + Rc^2)
```

- At r >> Rc: flat rotation curve with V = Vc
- At r << Rc: linear rotation (solid body), V = Vc * r / Rc
- The core radius Rc prevents a singular central density

| Galaxy | haloVc (km/s) | haloRc (arcsec) | haloRc (kpc) | haloRc (code) | M_halo/M_baryon |
|--------|---------------|-----------------|--------------|---------------|-----------------|
| M51a | 171.5 | 8 | 0.372 | 6.2 | 2.00 |
| M51b | 166.4 | 40 | 1.862 | 31.0 | 2.00 |

`haloVc` is the halo-only contribution, not the total observed rotation velocity. The total rotation curve is the quadrature sum of baryonic and halo terms.

---

## 6. Velocity Dispersion and Disc Stability

Both galaxies use Toomre Q = 1.5, the paper's standard value. This sets the radial and tangential velocity dispersion of disc particles:

```
sigma_r(r)   = Q * 3.36 * G * Sigma(r) / kappa(r)
sigma_phi(r) = sigma_r * kappa / (2 * Omega)
```

where `Sigma(r)` is the local exponential surface density and `kappa(r)` the epicyclic frequency from the full rotation curve.

Q = 1.5 gives a "warm" disc: stable against spontaneous fragmentation and particle-noise-driven multi-arm structure, while remaining responsive to the strong m=2 tidal perturbation. At Q < 1.2 particle noise drives incoherent modes (m = 2, 3, 4 comparable) that create spurious spiral structure before the encounter begins.

See `docs/toomre-q-velocity-dispersion.md` for the full derivation.

---

## 7. Particle Distribution

### 7.1 Counts and Resolution

`M51.sim` carries four commented resolution tiers; uncomment one.

| Tier | M51a | M51b | Total | Use |
|---|---|---|---|---|
| Minimal | 16,000 | 4,000 | 20,000 | Quick parameter checks |
| Low | 64,000 | 16,000 | 80,000 | Fast iteration |
| Medium | 256,000 | 64,000 | 320,000 | Reasonable morphology |
| High | 640,000 | 160,000 | 800,000 | Final renders |

M51a gets 4× the companion's count: it is the galaxy whose arms must be resolved, and its disc covers a larger area.

For reference, Salo & Laurikainen used 200,000 + 60,000 star particles plus 50,000 + 15,000 gas particles, noting this "rather coarse resolution ... [is] sufficient for the study of the gross features excited by tidal perturbation."

### 7.2 Softening

Softening should scale roughly as 1/sqrt(N) to keep two-body relaxation consistent:

| N (M51a) | Suggested r_soft |
|---|---|
| 64,000 | 0.3 |
| 256,000 | 0.15 |
| 640,000 | 0.095 |

`M51.sim` sets r_soft = 0.3. Reduce it when moving to a higher tier.

### 7.3 Disc Generation

`GalaxyDisc` places particles with:

- **Radial distribution** sampled from the exponential disc profile. The radial *number* density is the surface density times ring area, `p(r) ~ r*exp(-r/Re)`, which is a Gamma(2, Re) distribution — sampled exactly as `r = -Re * ln(u1*u2)`, with draws outside [Ri, Rd] rejected.
- **Circular velocity** from enclosed mass plus halo: `v_c^2 = G*M_enc(r)/r + haloVc^2 * r^2/(r^2 + haloRc^2)`
- **Velocity dispersion** from Toomre Q, varying with radius
- **Orientation** from the disc normal vector

M51b's disc normal is tilted 32.5 deg from +y:

```
normal = (sin(32.5deg), cos(32.5deg), 0) = (0.5373, 0.8434, 0.0)
```

---

## 8. Expected Behavior

### 8.1 Timeline

From the orbit integration (section 4.6). `Rd` = 18.62 kpc, `Re` = 4.65 kpc.

| Phase | Sim time | Physical time | What happens |
|-------|----------|---------------|--------------|
| Approach | t = 0–2.5 | 0–145 Myr | M51b falls in from apocenter, descending toward the disc plane from above (y = +506.5). |
| **Close crossing #1** | t ~ 2.87 | ~169 Myr | Crosses at 22.34 kpc = **1.200 Rd**. Principal spiral-arm excitation on an unperturbed disc. |
| First pericenter | t ~ 3.73 | ~219 Myr | Closest approach, 20.57 kpc, below the plane. |
| Arms develop | t = 4–7 | 235–410 Myr | Tidal arms wind up; bridge and tail form. |
| Apocenter | t ~ 7.47 | ~438 Myr | 30.86 kpc, maximum height off the plane. |
| Grazing crossing | t ~ 7.49 | ~440 Myr | Crosses at 1.66 Rd — weak, the companion is at its most distant. |
| Second pericenter | t ~ 11.20 | ~658 Myr | 20.57 kpc. |
| **Close crossing #2** | t ~ 12.09 | ~710 Myr | 22.44 kpc = **1.205 Rd** — the paper's `Rdown`. |
| **Best morphology** | **t = 12.8–13.5** | **750–790 Myr** | 40–80 Myr after `Rdown`. |
| Apocenter | t ~ 14.94 | ~877 Myr | |
| Close crossing #3 | t ~ 17.84 | ~1047 Myr | 22.25 kpc = 1.195 Rd. |

`End_Time` in `M51.sim` is 14.0, covering the best-morphology window with margin. Raise it to ~18 for the third close crossing.

Also worth inspecting: frames shortly after **close crossing #1** (t ~ 3–4), which acts on a clean disc and is the closest analogue to a single-passage encounter.

### 8.2 What to Look For

1. **Two-arm spiral arms** — coherent, symmetric, extending from the inner disc to the outer edge.
2. **Tidal bridge** — a stream of particles connecting M51a and M51b, pulled from M51a's disc during closest approach.
3. **Tidal tail** — an arm extending opposite the bridge, from angular momentum transfer.
4. **M51b distortion** — the companion develops asymmetry passing through M51a's tidal field.
5. **Arm winding** — after a crossing, inner arms rotate faster than outer ones, progressively winding the pattern.

### 8.3 Why the Arms Form

1. M51b's field creates a **tidal quadrupole** on M51a's disc, stretching it along the connecting line and compressing it perpendicular.
2. Disc material on the **near side** feels stronger attraction toward M51b than the disc center does, and is pulled outward.
3. Material on the **far side** feels weaker attraction than the center, so relative to the center it also moves outward.
4. Because the orbit is near-polar, the forcing is a series of impulsive crossings rather than sustained co-rotating perturbation. The paper's argument is that repeated crossings, not resonance, build the observed structure.
5. The **spiral** shape arises from differential rotation shearing the perturbation into a trailing pattern.

---

## 9. Comparison to Published Models

### 9.1 Salo & Laurikainen 2000

| Parameter | Paper (bound model) | This project | Match |
|---|---|---|---|
| Orbit class | Bound, multiple passage | Bound, multiple passage | yes |
| Eccentricity | ~0.2 | 0.200 | yes |
| Orbital inclination | 75–85 deg | 80 deg | yes |
| Crossing distance `Rcross` | 1.2–1.4 Rd | 1.195–1.205 Rd | yes |
| Apocenter between crossings | required | confirmed | yes |
| Mass ratio `Mp` | 0.5–0.7 (nom. 0.55) | 0.550 | yes |
| Toomre `Q_T` | 1.5 | 1.5 | yes |
| `M_disc/M_tot` within 4Re | 1/3 | 1/3 | yes |
| Primary Re / Rd | 100" / 400" | 100" / 400" | yes |
| Companion Re / Rd | 33" / 240" | 33" / 240" | yes |
| Halo Rc (primary / companion) | 8" / 40" | 8" / 40" | yes |
| Companion disc tilt | 32.5 deg | 32.5 deg | yes |
| Distance assumed | 9.6 Mpc | 9.6 Mpc | yes |
| Time between close crossings | 440–520 Myr | 541 Myr | 4% high |
| Halo truncation Rh | 400" (primary) | none (infinite) | **no** |
| Gas component | 65,000 sticky particles | none | **no** |
| Dynamical friction | modelled (sect 4) | none | **no** |
| Asymmetric drift correction | yes | no | **no** |

All the mass, geometry and orbital rows agree by construction. The remaining differences are missing code features, covered in section 9.3.

### 9.2 Dobbs et al. 2010

Their SPH + N-body model adds gas physics. They find `M_halo/M_disc` ~ 2.46, spiral arms visible ~200 Myr before pericenter, strongest arm-interarm contrast ~100–200 Myr after, with pericenter 25 kpc and mass ratio 1:3.

This collisionless simulation should reproduce the stellar arm morphology but not gas features (HII regions, dust lanes).

### 9.3 Limitations

Ordered by how much the paper suggests each matters.

1. **No halo truncation.** The paper truncates the primary halo at Rh = 400 arcsec (= Rd). Here `M_halo(r) ~ Vc^2 * r` grows without bound, so long-range attraction is overestimated. This is why the interval between close crossings is ~4% above the paper's `Tobs` range, and why the enclosed-mass ratio at the encounter separation approaches 1:1 when it should not.

2. **No gas.** The paper runs 50,000 + 15,000 dissipative "sticky" gas particles with fully inelastic collisions (α = 0), reaching 5–10 km/s radial dispersion. Gas produces sharp arm contrast, and the paper's morphological comparisons are largely made on the *gas* particles. Arms here will look smoother and more diffuse than the published figures.

3. **No dynamical friction.** The paper devotes section 4 to this, and it is load-bearing for their argument: friction means earlier passages were *more distant and weaker*, allowing a multiple-encounter history without destroying the disc. They find the observed well-defined far tail requires previous passages "at least 30 per cent more distant than the latest two crossings," achievable only with an extended halo (total `M_halo/M_disc` ~ 4–8). Here the orbit does not decay, so successive crossings are all at similar distance: the disc accumulates more heating than it should and the far tail is more dispersed.

4. **No asymmetric drift correction.** The paper corrects circular velocities for pressure support from random motions. Without it the disc is slightly out of equilibrium at t = 0 and relaxes over the first ~0.5 dynamical times. See `docs/toomre-q-velocity-dispersion.md` section 6.3.

5. **Rigid analytic halos.** Both this code and the paper's main runs use inert analytic halos, so this is not a divergence — but the paper *validates* against live self-consistent halo runs in section 4, which is how they confirm the friction result. Additionally, each halo center here tracks the mass-weighted centroid of all its member particles, so a long tidal tail drags the halo center off the nucleus. Momentum is conserved (section 10.4) but the field shape is affected.

6. **No vertical structure detail.** The paper initializes an isothermal sheet (sech^2) with σ_z/σ_r = 0.7. This code uses a simpler linear height scatter.

7. **No separate bulge in the paper's sense.** The bulge/disc split here (17% primary, 40% companion) is this project's own choice; the paper works with disc + halo.

---

## 10. Simulation Parameters

### 10.1 Time Step and Softening

- **dt = 0.0005** — required for the fast orbital velocities near M51a's center. At r = 78 code units the orbital period is ~2.3 code time units, so dt/T ~ 0.0002.
- **r_soft = 0.3** — gravitational softening (18 pc), preventing close-encounter singularities. Minimum resolved scale ~2 r_soft = 36 pc.
- **BH_Opening_Theta = 0.5** — Barnes-Hut opening angle, giving force errors under ~1%.

### 10.2 Warmup

`InitializationTime 2.0` starts the run at t = −2.0 with the systems dynamically isolated, letting each disc relax out of its initial particle-noise transients before the interaction begins. Bulk velocities are applied at t = 0, so every orbital parameter above still refers to t = 0. See `docs/script-format.md`.

Note that an isolated disc heats through two-body relaxation, which raises effective Q and makes it slightly less responsive to the tidal perturbation. Prefer shorter warmups if chasing maximum arm contrast.

### 10.3 PinCentralBodies

Each galaxy's central body is pinned to its system's center of mass every step (for systems with `halo_vc > 0`), preventing it from wandering due to N-body noise or recoiling from discrete particle encounters. The correction is momentum-conserving: the equal and opposite shift is spread over the system's other particles. The central body still feels and transmits the tidal field from the other galaxy.

### 10.4 Cross-Halo Gravity

Every body feels its own system's halo plus the halo of every other system, so M51a's particles feel M51b's halo and vice versa. This creates the correct tidal field for arm excitation.

Each halo is centered on its own system's mass-weighted particle barycenter, recomputed every derivative evaluation by `ComputeHaloCenters()`. Because a halo's field is nearly uniform across the other galaxy at these separations, the pair's relative acceleration comes out correct at `G(M_totA + M_totB)/d^2` including halo mass — which is why the analytic orbit of section 4 is reproduced by the simulation.

What the cross-halo term does **not** do is produce dynamical friction (section 5.1). It is a conservative central force between two rigid symmetric potentials.

### 10.5 Halo Monopole Removal

A rigid analytic halo has no inertia, so it cannot obey Newton's third law: particles are pulled toward the halo center and nothing pulls back. For an axisymmetric disc the per-particle forces cancel in the sum, but any asymmetry — above all an m=1 lopsided mode, or a tidal tail dragging the halo centroid off the nucleus — leaves a net force that accelerates the whole system off the origin.

The simulation subtracts the net halo force (divided by total mass) uniformly from every particle each step. Because the subtraction is uniform, every difference `a_i - a_j` is unchanged, so the relative orbit is preserved exactly and only the spurious bulk drift is cancelled. The correction is global, not per-system; a per-system correction would cancel the real mutual attraction between the galaxies.

This fixes the drift symptom, not the underlying halo-center definition: the halo center remains the mass-weighted centroid of all of that system's particles, so a long tidal tail still pulls it off the nucleus and distorts the field the disc feels (section 9.3).

---

## 11. Data Sources

| Source | Used For |
|--------|----------|
| Salo & Laurikainen 2000, MNRAS 319, 377 | All orbital and structural parameters (PDF in `docs/`) |
| Dobbs et al. 2010, MNRAS 403, 625 | Spiral arm evolution timeline, gas response |
| Querejeta et al. 2015, ApJS 219, 5 | M51a photometric mass decomposition (Spitzer 3.6μm) |
| Mentuch Cooper et al. 2012, ApJ 755, 165 | Stellar masses |
| Scoville & Young 1983 | Distance 9.6 Mpc (the paper's assumption) |
| McQuinn et al. 2016, ApJ 826, 21 | Modern TRGB distance 8.58 ± 0.10 Mpc (noted, not used) |
| Sofue et al. 1999, ApJ 523, 136 | M51a rotation curve |
| Meidt et al. 2013, ApJ 779, 45 | M51a mass distribution from CO kinematics |
| Persic & Salucci 1990 | Halo dominance in small galaxies |
| Spillar et al. 1992 | Companion scale length |
| Tully 1974 | M51a disc position angle and inclination |
| Toomre & Toomre 1972, ApJ 178, 623 | Original tidal tail theory |
| Toomre 1964, ApJ 139, 1217 | Disc stability criterion (Q parameter) |
| de Vaucouleurs et al. 1991 (RC3) | Galaxy classifications and sizes |

---

## 12. Running the Simulation

Load `M51.sim`. The camera is positioned along +y, face-on to M51a's disc.

- **t = −2 to 0** — warmup. Systems isolated and settling; not recorded to video.
- **t ~ 2.9** — first close crossing. Arm excitation begins.
- **t ~ 3–4** — good morphology on a clean disc.
- **t ~ 12.1** — second close crossing.
- **t = 12.8–13.5** — the paper's nominal observation epoch, the closest analogue to M51 as seen today.

For comparison against observations, note that M51a is seen from Earth at an inclination of ~20 deg from face-on with position angle ~170 deg (Tully 1974), which the paper uses for its projected figures. The default camera here is exactly face-on.
