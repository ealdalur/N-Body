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
| `Tobs` — principal crossing → observation | 5.5–6.5 (440–520 Myr) | Table 2 |
| `Tobs − Tdown` — latest crossing → observation | 0.5–1.0 (40–80 Myr) | Table 2 |
| **`Tdown` — crossing → crossing** (derived) | **4.5–6.0 (360–480 Myr)** | Table 2 |
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

This gives `M_halo/M_baryon` = 2 within the disc radius **for the primary**.

The **companion is not decomposed this way**. As the fainter galaxy it is made more halo-dominated (Persic & Salucci 1990, `M_halo/M_disc ∝ L_disc^-0.5`). The paper sets its masses directly (sect 2.2): `M_disc = 0.13` and `M_halo = 0.42` in units of the primary total within Rd, so `M_halo/M_baryon = 3.23`, not 2. See §3.3.

### 2.2 Baryonic Mass vs Photometric Mass

These baryonic masses are **rotation-curve-decomposed** — the dynamically cold disc mass that participates in spiral arm formation. For the primary the decomposition (1/3 of a `V` = 220 km/s rotation curve) gives `M_baryon` ≈ 7.0×10^10 M☉, matching the paper's disc-mass fraction of its ~2×10^11 M☉ total within Rd. This is comparable to or somewhat above photometric stellar estimates (Querejeta et al. 2015 gives M_stellar ~ 5×10^10 M☉ for M51a); the two are not the same quantity — the dynamical decomposition is fixed by the assumed rotation curve and disc/halo split, not by photometry.

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
| Observed flat V_rotation | 220 km/s | S&L sect 2.2/2.3 (their velocity unit) | — |
| V_baryon at Rd | 127.0 km/s | = 220/sqrt(3) | — |
| V_halo at Rd | 179.6 km/s | = 220*sqrt(2/3) | — |
| Total baryonic mass | 6.98 × 10^10 M☉ | = V_b^2 * Rd / G | 5,005,858 |
| Bulge mass (M) | 1.19 × 10^10 M☉ | 17% of baryonic (Sbc) | 850,000 |
| Disc mass | 5.80 × 10^10 M☉ | baryonic − bulge | 4,155,858 |
| Mfrac (disc/bulge) | 4.89 | — | — |
| DM halo Vc | 179.7 km/s | corrected for core | — |
| DM halo core Rc | 8 arcsec = 0.372 kpc | S&L sect 2.2 | 6.2 |
| M_halo within Rd | 1.40 × 10^11 M☉ | cored isothermal | 10,011,716 |
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
| Total mass within Rd | 1.15 × 10^11 M☉ | = Mp × M_prim | 8,259,666 |
| Total baryonic mass | 2.72 × 10^10 M☉ | Mdisc = 0.13 (0.236 of Mp) | 1,952,285 |
| Bulge mass (M) | 1.09 × 10^10 M☉ | 40% of baryonic (SB0) | 780,000 |
| Disc mass | 1.64 × 10^10 M☉ | baryonic − bulge | 1,172,285 |
| Mfrac (disc/bulge) | 1.50 | — | — |
| DM halo Vc | 186.6 km/s | inverted from M_halo target | — |
| DM halo core Rc | 40 arcsec = 1.862 kpc | S&L sect 2.2 | 31.0 |
| M_halo within Rd | 8.80 × 10^10 M☉ | Mhalo = 0.42 (0.764 of Mp) | 6,307,381 |
| M_halo / M_baryon | 3.23 | paper Mdisc=0.13/Mhalo=0.42 | — |
| Implied V_rotation | 211 km/s | *derived*, not assumed | — |
| Disc tilt to M51a disc | 32.5 deg | S&L sect 2.2 (observed) | — |
| Inner hole | 0.2 kpc | — | 3.3 |
| Toomre Q | 1.5 | S&L standard value | — |

The paper notes `Mp` = 0.55 deliberately exceeds the *observed disc* mass ratio of ~0.4, because "small galaxies are likely to have more dominant halo components" (citing Persic & Salucci 1990). This is also why the companion's own disc/halo split is more halo-heavy than the primary's: the paper sets `Mdisc` = 0.13 and `Mhalo` = 0.42 (in units of the primary total), i.e. `M_halo/M_baryon` = 3.23 rather than the primary's 2.0. The 211 km/s implied rotation curve is higher than M51b's observed stellar kinematics suggest — intentional in the source model.

The larger halo core (40 arcsec vs the primary's 8) reflects "the lesser degree of halo concentration expected for smaller galaxies."

**Disc tilt vs orbital inclination** are independent quantities: the 32.5 deg disc tilt comes from the companion's observed position angle and inclination, while the 80 deg orbital inclination is a property of the trajectory.

### 3.4 Mass Ratio Verification

`Mp` = M_tot(companion) / M_tot(primary) within 4Re = Rd:

```
M51a total within Rd:  5,005,858 (baryon) + 10,011,716 (halo) = 15,017,574
M51b total within Rd:  1,952,285 (baryon) +  6,307,381 (halo) =  8,259,666

Mp = 8,259,666 / 15,017,574 = 0.550      <- paper nominal, range 0.5-0.7
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

where `M_baryon` = 5,005,858 + 1,952,285 = 6,958,143.

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

position = r_apo * e2 = (0, 458.4, 80.8)    apocenter, max height off plane
velocity = +v_t * e1  = (+200.0, 0, 0)      purely tangential, in the disc plane
```

The **sign of the velocity** sets the orbit's circulation sense and is not free. `iorb` is defined as the angle between the orbital angular momentum `L_orb = pos × vel` and the primary's disc angular momentum `L_disc`. With `+e1`:

```
L_orb  = (0, +16160, -91680)
L_disc = (0, +1, 0)              M51a's disc normal
angle  = 80.0 deg                matches iorb
```

Using `−e1` mirrors this to `180 − iorb` = 100 deg — the same near-polar geometry but counter-circulating, and outside the paper's 75–85 deg range.

The script asserts three things: `|pos|` equals the apocenter the velocity was solved for, `pos · vel` ≈ 0 (apocenter has no radial velocity), and the realized `iorb` recomputed from the emitted state vector equals the requested value. The last check is what catches a sign or basis error, which the inclination input alone cannot.

#### Three independent rotation senses

These are easy to conflate but are set separately:

| | Sense | Set by | Value |
|---|---|---|---|
| 1 | M51a's stars about M51a's centre | M51a disc normal | `L` along +y |
| 2 | M51b's stars about M51b's centre | M51b disc normal | `L` 32.5 deg from #1, same sense |
| 3 | M51b's orbit about M51a | orbital state vector | `L_orb` at `iorb` = 80 deg to #1 |

`GalaxyDisc` gives each disc a spin angular momentum parallel to its own normal (the generator uses `v_tan = -vm * cross(r_hat, n)`, which yields `L ∥ +n`). So #1 and #2 follow directly from the normal vectors in the script, and the 32.5 deg between them is the paper's observed relative disc inclination. Only #3 depends on the orbital velocity sign.

### 4.5 Computed Parameters

```
target Rcross = 1.20 Rd                <- input
  -> pericenter = 18.62 kpc            <- solved by bisection
  -> apocenter  = r_peri * (1+e)/(1-e) = 27.93 kpc
```

| Parameter | Value |
|-----------|-------|
| Target Rcross | 1.20 Rd (22.34 kpc) |
| Halo truncation Rh | Rd for each galaxy (paper sect 2.2) |
| Pericenter (solved) | 310.3 code units (18.62 kpc) |
| Apocenter (derived) | 465.4 code units (27.93 kpc) |
| Eccentricity | 0.200 |
| Orbital inclination | 80 deg |
| Argument of pericenter | 90 deg |
| Tidal strength S = Mp*(Rd/r)^3 | 0.318 |
| Tangential velocity at apocenter | 200.0 km/s |
| Velocity at pericenter | 300.0 km/s |
| v_circular at apocenter | 286.2 km/s |
| v_t / v_circ | 0.699 |
| Orbital period | ~9.95 code units (~584 Myr radial) |

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
| **crossing #1** | 3.12 | 183 | 22.34 | **1.200** |
| pericenter | 4.97 | 292 | 18.62 | 1.000 |
| **crossing** (`Rcross`) | 6.83 | 401 | 22.34 | **1.200** |
| apocenter | 9.95 | 584 | 27.93 | 1.500 |
| **crossing** (`Rdown`) | 13.06 | 767 | 22.34 | **1.200** |
| pericenter | 14.92 | 876 | 18.62 | 1.000 |
| **crossing** | 16.78 | 985 | 22.34 | **1.200** |

Every close crossing lands at exactly 1.200 Rd, inside the paper's `Rcross` = 1.2–1.4 and `Rdown` = 1.2–1.3. This is an independent check: pericenter was the targeted quantity and the crossing radius follows.

**Why every crossing is at the same radius.** Solving `Rcross` = 1.20 Rd with e = 0.2 puts the pericenter at exactly 1.000 Rd = 18.62 kpc — the primary's halo truncation radius `Rh`. The whole orbit therefore lies in the truncated regime, where both halos act as point masses and the potential is Keplerian. A Kepler ellipse does not precess, so with the argument of pericenter fixed at 90 deg the two nodes sit symmetrically at true anomaly ±90 deg from pericenter, both at the semi-latus rectum `r = r_peri(1+e)` = 18.62 × 1.2 = **22.34 kpc = 1.200 Rd**. Successive crossings repeat that radius exactly rather than alternating close/grazing as they would in a precessing (untruncated, logarithmic) potential.

### 4.7 Best-Morphology Epoch

The paper anchors its observation epoch `Tobs` to the **principal disc-plane crossing**: "the time Tobs elapsed since the disc plane crossing at T = 0." Section 3.4.3 brackets it from both sides:

- `Tobs` ≤ 4–5 — spiral structure and the long HI tail too weak.
- `Tobs` ~ 5–6 — produces the observed **"kink" on the northern arm** toward the companion, from the interplay between the most recent crossing and arms excited by the first passage.
- `Tobs` ~ 7 — that arm "would be opened into a new bridge during the second crossing," overshooting.
- `Tobs` > 6 — the spiral becomes too tightly wound.

Applying `Tobs - Tdown` = 0.5–1.0 (40–80 Myr, i.e. 0.68–1.36 of this project's 58.7 Myr unit) after the most recent crossing `Rdown` at t = 13.06:

**Observe at t = 13.74 to 14.43, centre ~14.09.**

**Crossing-to-crossing timing.** The target is `Tdown` = 4.5–6.0 of the paper's 80 Myr unit = **360–480 Myr**. Note `Tobs` (440–520 Myr) and `Tobs − Tdown` (40–80 Myr) are both crossing-to-*observation* intervals, not this quantity.

Identifying the right pair of crossings matters. A near-polar orbit crosses the disc plane twice per orbit, splitting it into two unequal arcs — one containing pericenter, one containing apocenter — which sum to the full radial period:

| Crossing pair | Gap | Apsis between |
|---|---|---|
| t = 3.12 → 6.83 | 218 Myr | pericenter |
| **t = 6.83 → 13.06** | **366 Myr** | **apocenter** ← the paper's pair |
| t = 13.06 → 16.78 | 218 Myr | pericenter |

The paper's Fig. 1 caption fixes which pair it means: *"in the multiple-passage model the apocentre is between the two disc crossings."* So `Rcross` = t 6.83 and `Rdown` = t 13.06, giving **366 Myr — inside the 360–480 Myr target.**

### 4.8 Coordinate System

- M51a disc lies in the **x-z plane**, disc normal **+y**
- M51a disc rotates **counter-clockwise** about its normal: `LoadGalaxyDiscState` uses `v_tan = -vm` with `t_hat = cross(r_hat, n_hat)`, which puts the disc angular momentum `L_disc` along **+y**, parallel to the normal. With the camera at +y looking toward the origin this also reads counter-clockwise on screen.
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
| M51a | 179.7 | 8 | 0.372 | 6.2 | 2.00 |
| M51b | 186.6 | 40 | 1.862 | 31.0 | 3.23 |

`haloVc` is the halo-only contribution, not the total observed rotation velocity. The total rotation curve is the quadrature sum of baryonic and halo terms.

---

## 6. Velocity Dispersion and Disc Stability

Both galaxies use Toomre Q = 1.5, the paper's standard value. This sets the radial and tangential velocity dispersion of disc particles:

```
sigma_r(r)   = Q * 3.36 * G * Sigma(r) / kappa(r)
sigma_phi(r) = sigma_r * kappa / (2 * Omega)
sigma_z(r)   = 0.7 * sigma_r                      (Salo & Laurikainen sect 2.2)
```

where `Sigma(r)` is the local exponential surface density and `kappa(r)` the epicyclic frequency from the full rotation curve. `Sigma(r)` is normalized to the **truncated** disc mass — `Sigma_0 = M_disc / (2*pi*h_r^2 * enc_denom)`, with `enc_denom = 1 - (1 + Rd/h_r) exp(-Rd/h_r)` — so it integrates to `M_disc` over `[0, Rd]`. (Omitting `enc_denom` would under-normalize `Sigma` by ~9% for the primary and cool the disc to an effective Q ~1.36.)

**Asymmetric drift.** The mean streaming (tangential) velocity is not the full circular speed: pressure support from the random motions lowers it, per the Binney & Tremaine asymmetric drift equation

```
v_c^2 - <v_phi>^2 = sigma_r^2 [ sigma_phi^2/sigma_r^2 - 1 - d ln(Sigma sigma_r^2)/d ln R ]
```

`LoadGalaxyDiscState` evaluates the log-derivative numerically and draws the Gaussian `sigma_phi` scatter about this reduced mean, so the disc starts closer to equilibrium.

**Vertical structure.** Particle heights follow the self-gravitating isothermal sheet `rho(z) = rho0 sech^2(z/z0)` with `z0 = sigma_z^2/(2*pi*G*Sigma)` (Spitzer 1942), sampled by inverting its CDF, and each particle is given a Gaussian vertical velocity of dispersion `sigma_z`.

Q = 1.5 gives a "warm" disc: stable against spontaneous fragmentation and particle-noise-driven multi-arm structure, while remaining responsive to the strong m=2 tidal perturbation. Below about Q = 1.2 particle noise begins driving incoherent modes (m = 2, 3, 4 comparable) that create spurious spiral structure before the encounter begins.

Note the warmup phase interacts with this: an isolated disc heats through two-body relaxation, so the effective Q at t = 0 is somewhat above the value set here.

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
- **Circular velocity** from enclosed mass plus halo: `v_c^2 = G*M_enc(r)/r + haloVc^2 * r^2/(r^2 + haloRc^2)`, reduced to the mean streaming velocity by the asymmetric drift correction
- **In-plane velocity dispersion** from Toomre Q, varying with radius (`sigma_r`, `sigma_phi`)
- **Vertical structure** as a sech^2 isothermal sheet with `sigma_z = 0.7 sigma_r`
- **Orientation** from the disc normal vector

M51b's disc normal is tilted 32.5 deg from +y:

```
normal = (sin(32.5deg), cos(32.5deg), 0) = (0.5373, 0.8434, 0.0)
```

---

## 8. Expected Behavior

### 8.1 Timeline

From the orbit integration (section 4.6). `Rd` = 18.62 kpc, `Re` = 4.65 kpc, halos truncated at each galaxy's own `Rd`.

| Phase | Sim time | Physical time | What happens |
|-------|----------|---------------|--------------|
| Approach | t = 0–3 | 0–176 Myr | M51b falls in from apocenter, descending toward the disc plane from above (y = +458.4). |
| Crossing | t ~ 3.12 | ~183 Myr | Crosses at 22.34 kpc = 1.200 Rd. |
| Pericenter | t ~ 4.97 | ~292 Myr | 18.62 kpc, below the plane. |
| **`Rcross`** | t ~ 6.83 | ~401 Myr | Principal crossing, 1.200 Rd — the paper's spiral-inducing passage. |
| Apocenter | t ~ 9.95 | ~584 Myr | 27.93 kpc, **between** the `Rcross`/`Rdown` pair. |
| **`Rdown`** | t ~ 13.06 | ~767 Myr | Most recent crossing, 1.200 Rd. 366 Myr after `Rcross`. |
| **Best morphology** | **t = 13.74–14.43** | **807–847 Myr** | 40–80 Myr after `Rdown`, per `Tobs − Tdown`. |
| Pericenter | t ~ 14.92 | ~876 Myr | |
| Crossing | t ~ 16.78 | ~985 Myr | 1.200 Rd. |

With the halo truncated every crossing occurs at the same radius, 1.200 Rd.

`End_Time` in `M51.sim` is 16.0, covering the best-morphology window.

**Choosing the pair matters.** The orbit crosses the plane twice per revolution, and the two arcs (218 Myr containing pericenter, 366 Myr containing apocenter) sum to the radial period. Only the apocenter-bracketing pair is the paper's `Rcross`/`Rdown`, so the earlier crossing at t = 3.12 and the pair t = 3.12 → 6.83 are not the constrained quantity.

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
| Orbital inclination `iorb` | 75–85 deg | 80.0 deg (verified) | yes |
| Crossing distance `Rcross` | 1.2–1.4 Rd | 1.200 Rd (all crossings) | yes |
| Apocenter between crossings | required | confirmed | yes |
| Mass ratio `Mp` | 0.5–0.7 (nom. 0.55) | 0.550 | yes |
| Rotation velocity (primary) | 220 km/s | 220 km/s | yes |
| Companion `M_halo/M_disc` | 3.23 (Mdisc=0.13/Mhalo=0.42) | 3.23 | yes |
| Toomre `Q_T` | 1.5 | 1.5 | yes |
| `M_disc/M_tot` within 4Re (primary) | 1/3 | 1/3 | yes |
| Primary Re / Rd | 100" / 400" | 100" / 400" | yes |
| Companion Re / Rd | 33" / 240" | 33" / 240" | yes |
| Halo Rc (primary / companion) | 8" / 40" | 8" / 40" | yes |
| Companion disc tilt | 32.5 deg | 32.5 deg | yes |
| Distance assumed | 9.6 Mpc | 9.6 Mpc | yes |
| Crossing-to-crossing (`Tdown`) | 360–480 Myr | 366 Myr | yes |
| Halo truncation `Rh` | Rd (400" primary) | Rd (400" primary) | yes |
| Vertical structure | sech^2, σz/σr = 0.7 | sech^2, σz/σr = 0.7 | yes |
| Asymmetric drift correction | yes | yes | yes |
| Gas component | 65,000 sticky particles | none | **no** |
| Dynamical friction | modelled (sect 4) | none | **no** |

All the mass, geometry and orbital rows agree by construction. The remaining differences are missing code features, covered in section 9.3.

### 9.2 Dobbs et al. 2010

Their SPH + N-body model adds gas physics. They find `M_halo/M_disc` ~ 2.46, spiral arms visible ~200 Myr before pericenter, strongest arm-interarm contrast ~100–200 Myr after, with pericenter 25 kpc and mass ratio 1:3.

This collisionless simulation should reproduce the stellar arm morphology but not gas features (HII regions, dust lanes).

### 9.3 Limitations

Ordered by how much the paper suggests each matters.

1. **No dynamical friction.** In reality M51b raises a trailing wake in M51a's halo and loses orbital energy (Chandrasekhar friction), so earlier passages were *more distant and weaker*. The analytic halos here are rigid spherical potentials comoving with their own galaxy's barycenter, symmetric by construction, so they exert no net drag; the only friction present comes from the live baryonic particles, a minority of the mass. The orbit therefore does not shrink between passages.

   This is **not** a departure from the paper's nominal model, which is likewise frictionless: *"only the discs were self-gravitating, [so] the orbital energy could not be transformed to the halo deformation and the amount of friction was thus strongly underestimated."* Section 4 adds friction separately — via the Chandrasekhar formula and via live-halo runs — to reconstruct the *past* orbital history, using a companion *"modelled only by a halo"* and much larger halo extents (`Rh` = 2–4 length units versus 1 here). Its conclusion is that with a sufficiently extended halo the previous passages are distant enough that *"the effects of the previous passages on the final morphology are almost negligible"* — which is what justifies the frictionless nominal model going forward.

   What it does cost: all passages occur at the same distance rather than progressively closer, so the cumulative disc heating is higher than in reality and the far tail comes out more dispersed. The paper finds the observed well-defined far tail requires previous passages "at least 30 per cent more distant than the latest two crossings."

2. **No gas.** The paper runs 50,000 + 15,000 dissipative "sticky" gas particles with fully inelastic collisions (α = 0), reaching 5–10 km/s radial dispersion. Gas produces sharp arm contrast, and the paper's morphological comparisons are largely made on the *gas* particles. Arms here will look smoother and more diffuse than the published figures.

3. **Rigid analytic halos.** Both this code and the paper's main runs use inert analytic halos, so this is not a divergence — but the paper *validates* against live self-consistent halo runs in section 4, which is how they confirm the friction result. Additionally, each halo center here tracks the mass-weighted centroid of all its member particles, so a long tidal tail drags the halo center off the nucleus. Momentum is conserved (section 10.4) but the field shape is affected.

4. **No separate bulge in the paper's sense.** The bulge/disc split here (17% primary, 40% companion) is this project's own choice; the paper works with disc + halo.

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
- **t ~ 3.1** — first close crossing. Arm excitation begins.
- **t ~ 3–5** — good morphology on a clean disc.
- **t ~ 13.1** — most recent close crossing (`Rdown`).
- **t = 13.7–14.4** — the paper's nominal observation epoch, the closest analogue to M51 as seen today.

For comparison against observations, note that M51a is seen from Earth at an inclination of ~20 deg from face-on with position angle ~170 deg (Tully 1974), which the paper uses for its projected figures. The default camera here is exactly face-on.
