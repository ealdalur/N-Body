# compute_apollo_free_return.py

## 1. Overview

This script computes initial conditions for a three-body simulation of an Apollo-style free-return trajectory in the Earth-Moon system. The spacecraft departs Earth after Trans-Lunar Injection (TLI), coasts to the Moon on a ballistic trajectory, swings around the lunar far side via gravitational slingshot, and returns to Earth — all without additional propulsive maneuvers.

The free-return trajectory was used by Apollo 8, 10, 11, 12, and 13. Apollo 13 famously relied on it for safe return after the service module explosion.

## 2. Unit System (G=1)

Unlike the galactic simulations (which use 60 pc / 10⁴ M☉ / 1 km/s), this simulation uses units natural to Earth-centered orbital mechanics:

| Dimension | 1 code unit = | Physical meaning |
|-----------|---------------|------------------|
| Distance (du) | 6,371 km | Earth mean radius (R⊕) |
| Mass (mu) | 5.972 × 10²⁴ kg | Earth mass (M⊕) |
| Velocity (vu) | 7,909.5 m/s = 7.91 km/s | Circular orbit speed at R⊕ |
| Time (tu) | 805.5 s = 13.4 min | Derived: √(R⊕³ / (G·M⊕)) |

### 2.1 Derivation

The time unit comes from Kepler's equation at the Earth's surface:

```
tu = √(R⊕³ / (G · M⊕))
   = √((6.371×10⁶)³ / (6.674×10⁻¹¹ × 5.972×10²⁴))
   = √(2.586×10²⁰ / 3.986×10¹⁴)
   = √(6.489×10⁵)
   = 805.5 s
```

The velocity unit follows:
```
vu = R⊕ / tu = 6.371×10⁶ m / 805.5 s = 7,909.5 m/s
```

### 2.2 Verification

- Circular orbit speed at r = 1 du: v = √(GM/r) = √(1) = 1.0 vu = 7.91 km/s ✓
- Escape velocity at r = 1 du: v = √(2) = 1.414 vu = 11.18 km/s ✓
- Moon orbital distance: 384,400 / 6,371 = 60.3 du (matches "60 Earth radii") ✓

### 2.3 Why G=1?

Setting G=1 eliminates one parameter from the simulation and makes orbital mechanics equations dimensionally clean: `v_circular = √(M/r)` and `v_escape = √(2M/r)` with no prefactors.

## 3. Coordinate Convention

The orbital plane is **x-z**. The **+y** axis points out of the plane (toward Earth's north pole).

Counter-clockwise (CCW) motion viewed from +y obeys:

```
angular velocity:  ω = (0, +ω, 0)

position at angle θ (CCW from +x):  (r·cos θ, 0, -r·sin θ)
velocity at angle θ:                 v·(-sin θ, 0, -cos θ)
```

This follows directly from the cross product ω×r in a right-handed coordinate system:

```
At θ=0 (on +x axis):  position = (r, 0, 0),  velocity = (0, 0, -v)
At θ=90° (on -z axis): position = (0, 0, -r), velocity = (-v, 0, 0)
```

The key sign to remember: **at the +x axis, CCW velocity is in the -z direction.**

## 4. Physical Data

### 4.1 Earth

| Parameter | Value | Source |
|-----------|-------|--------|
| Mass | 5.972 × 10²⁴ kg | IAU 2009 |
| Mean radius | 6,371 km | WGS-84 |

### 4.2 Moon

| Parameter | Value | Code units | Source |
|-----------|-------|------------|--------|
| Mass | 7.342 × 10²² kg | 0.01229 mu | IAU 2009 |
| Radius | 1,737 km | 0.2726 du | IAU |
| Mean distance | 384,400 km | 60.34 du | Lunar Laser Ranging (Williams et al. 2014) |
| Orbital period | 27.32 days | 2,927 tu | Sidereal |
| Orbital velocity | 1.012 km/s | 0.1280 vu | Derived (see §5.1) |
| Sphere of influence | 66,170 km | 10.4 du | r·(m/M)^(2/5) |
| Mass ratio M⊕/M☽ | 81.3 | — | IAU |

### 4.3 Apollo TLI Parameters

| Parameter | Value | Code units | Source |
|-----------|-------|------------|--------|
| Burnout altitude | 334 km | — | Apollo 11 Mission Report (NASA SP-238) |
| Burnout radius | 6,705 km | 1.0524 du | altitude + R⊕ |
| Burnout velocity | 10,834 m/s | 1.3697 vu | Apollo 11 S-IVB cutoff |
| Escape velocity (at burnout) | 10,899 m/s | 1.3779 vu | √(2/r_tli) |
| v_TLI / v_escape | 0.994 | — | — |
| Transit time to Moon | ~72 hours | ~320 tu | Apollo 8, 10, 11, 13 |
| Free-return perilune | 254 km | — | Apollo 13 (NASA MSC-02680) |

The TLI velocity is 99.4% of escape velocity at that altitude. The spacecraft is on a bound, highly eccentric ellipse whose apogee overshoots the Moon. The Moon's gravity intervenes before apogee.

## 5. Orbital Mechanics Calculations

### 5.1 Moon's Orbital Velocity

The Moon and Earth orbit their common barycenter. Using the two-body problem:

```
Relative velocity: v_rel = √(G·(M⊕ + M☽) / r_moon)
                         = √((1 + 0.01229) / 60.34)
                         = 0.12955 vu = 1.025 km/s

Moon velocity (COM frame): v_moon = v_rel × M⊕/(M⊕ + M☽)
                                  = 0.12955 × 1.0/1.01229
                                  = 0.12798 vu = 1.012 km/s

Earth velocity (COM):      v_earth = v_rel × M☽/(M⊕ + M☽)
                                   = 0.001573 vu = 12.4 m/s
```

The Earth's tiny counter-velocity maintains the center of mass at rest.

### 5.2 Transfer Orbit Elements

With TLI as the perigee of a Keplerian ellipse around Earth:

```
Specific energy: E = v²/2 - 1/r = (1.3697)²/2 - 1/1.0524 = -0.0122
Semi-major axis: a = -1/(2E) = 41.36 du = 263,500 km
Eccentricity:    e = 1 - r_perigee/a = 1 - 1.0524/41.36 = 0.9746
Apogee distance: r_a = a(1+e) = 81.7 du = 520,400 km
Semi-latus rectum: p = a(1-e²) = 2.056 du
```

The apogee (81.7 du) is ~35% beyond the Moon's orbit (60.3 du). In an unperturbed two-body orbit, the spacecraft would overshoot the Moon and return on the same ellipse. The Moon's gravity provides the essential third-body perturbation.

### 5.3 True Anomaly at Moon's Distance

The spacecraft crosses the Moon's orbital distance (60.34 du) BEFORE reaching apogee. The true anomaly at that crossing is found from the orbit equation:

```
r = p / (1 + e·cos ν)

cos ν = (p/r - 1) / e = (2.056/60.34 - 1) / 0.9746 = -0.9917

ν = 172.6°
```

**This is 7.4° short of apogee (180°).** This distinction is critical for targeting — at 60 Earth radii, a 7.4° angular error corresponds to a ~52,000 km miss.

### 5.4 Transit Time (Kepler's Equation)

Converting true anomaly to time requires the eccentric anomaly and mean anomaly:

```
tan(E/2) = √((1-e)/(1+e)) · tan(ν/2)

E = eccentric anomaly → M = E - e·sin(E) → t = M/(2π) × T_orbit

Transfer orbit period: T = 2π√(a³) = 2π√(41.36³) = 1,670 tu = 15.6 days
Transit time to ν=172.6°: t_transit ≈ 321 tu = 71.8 hours ≈ 3.0 days
```

This matches the observed Apollo transit times of 72-76 hours.

### 5.5 Moon Lead Angle

During the 321 tu transit, the Moon moves in its orbit:

```
Moon angular velocity: ω = 2π/T_orbit = 2π/2,927 = 0.002146 rad/tu
Moon displacement: ω × 321 = 0.689 rad = 39.5° CCW
```

The targeting must account for this 39.5° lead angle — aim at where the Moon WILL BE, not where it is now.

### 5.6 Targeting Logic (Critical Subtlety)

The naive approach is: "aim the orbit's apogee at the Moon's future position." This is WRONG for a highly eccentric orbit.

The correct logic:

1. The spacecraft crosses r_moon at true anomaly ν = 172.6° from perigee (not 180°).
2. The spacecraft's angular position at crossing is: `perigee_angle + 172.6°`
3. This must equal the Moon's angular position at the time of crossing.

Therefore:
```
perigee_angle + ν_target = moon_angle_at_encounter
perigee_angle = moon_angle_at_encounter - ν_target
             = 39.5° - 172.6° = -133.1° (mod 360° = 226.9°)
```

If you mistakenly use 180° (apogee) instead of 172.6° (actual crossing):
```
Wrong: perigee_angle = 39.5° - 180° = -140.5°
Error: 7.4° × (60.34 du × 6371 km/du × π/180) ≈ 50,000 km miss
```

This 50,000 km error is enough to miss the Moon's sphere of influence entirely.

### 5.7 Flyby Offset (Impact Parameter)

A spacecraft aimed directly at the Moon's center would collide. Adding a small angular offset (0.7°) shifts the aim point behind the Moon in its orbit, creating an impact parameter that determines perilune altitude:

```
flyby_offset_deg = 0.7°
perigee_angle += 0.7° (shift aim behind Moon)

Impact parameter: b ≈ r_moon × sin(0.7°) = 60.34 × 0.01222 = 0.737 du ≈ 4,700 km
```

This is the most sensitive parameter in the entire calculation. A 0.1° change shifts the perilune by hundreds of kilometers.

The value 0.7° was chosen to approximate the Apollo 13 free-return perilune of 254 km above the lunar surface.

## 6. Encounter Analysis

### 6.1 Why the Spacecraft "Stalls" at the Moon

The spacecraft's speed at the Moon's distance (vis-viva equation):

```
v_spacecraft = √(2/r - 1/a) = √(2/60.34 - 1/41.36) = 0.095 vu = 0.75 km/s
```

The Moon's orbital speed: **1.01 km/s**

The spacecraft is SLOWER than the Moon. This is physically correct — the spacecraft is near the apogee of its transfer ellipse, where orbital mechanics dictates minimum speed. In the simulation, this manifests as the spacecraft appearing to "stall" while the Moon catches up.

### 6.2 The Slingshot Mechanism

In the Moon's reference frame, the spacecraft approaches with a relative velocity:

```
Radial component: v_r = √(1/p) × e × sin(ν) = 0.014 vu = 0.11 km/s
Tangential component: v_t = √(1/p) × (1 + e·cos(ν)) = 0.006 vu = 0.05 km/s
Moon's tangential velocity: v_moon = 0.128 vu = 1.01 km/s

Relative speed = √(v_r² + (v_t - v_moon)²) ≈ 0.14 vu ≈ 1.1 km/s
```

This relative velocity determines the hyperbolic flyby trajectory within the Moon's sphere of influence:

```
Estimated perilune: ~2,100 km from Moon center (~370 km altitude)
Hyperbolic eccentricity: e_hyp ≈ 7.5
Deflection angle: δ = 2·arcsin(1/e_hyp) ≈ 15°
```

**Note:** These are rough two-body analytical estimates. The actual three-body dynamics determine the true flyby geometry, which can differ significantly due to Earth's continuing gravitational influence within the Moon's SOI.

### 6.3 Return Trajectory

After the flyby redirects the velocity vector, the spacecraft is on a new Earth-centered ellipse that returns to low Earth altitude. The "dramatic acceleration" on the return is simply falling back into Earth's gravity well — by energy conservation, the return speed at any altitude equals the departure speed at that altitude (~11 km/s near Earth's surface).

The full slingshot process takes ~30 hours (while traversing the Moon's SOI of ~66,000 km at ~1 km/s relative velocity). This extended timescale is why the flyby looks gradual rather than instantaneous.

## 7. Initial Conditions (Final Configuration)

### 7.1 Earth

```
Position: (0, 0, 0)
Velocity: (0, 0, +0.001573)  [COM balance: +z because Moon goes -z]
Mass: 1.0
```

### 7.2 Moon

```
Position: (60.3359, 0, 0)     = 384,400 km along +x axis
Velocity: (0, 0, -0.127955)   = 1.012 km/s in -z direction (CCW orbit)
Mass: 0.012294
```

The Moon starts along +x. Its velocity is in -z because at θ=0 (on +x axis), CCW motion gives velocity in the -z direction (see §3).

### 7.3 Spacecraft (Apollo CSM)

```
Position: (-0.70710, 0, 0.77949)   = (-4,505, 0, 4,966) km
Velocity: (1.01451, 0, 0.92030)    = (8.024, 0, 7.279) km/s
Mass: 1.0e-20                       (test particle)
|v| = 1.3697 vu = 10,834 m/s
```

The position is at 334 km altitude (r = 1.0524 du) at the perigee angle computed in §5.6. The velocity is purely tangential (perpendicular to radius) at TLI burnout speed, directed to send the spacecraft toward the Moon's future position via the transfer orbit.

## 8. Simulation Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| dt | 0.01 tu (8.1 s) | Resolves the lunar flyby (perilune passage takes ~10-50 tu) |
| r_soft | 0.01 du (64 km) | Prevents singularities; smaller than Moon radius (0.27 du) |
| Solver | LeapFrog | Symplectic — conserves energy over the multi-day trajectory |
| Gravity | P2P | Only 3 bodies — exact pairwise force computation |
| Sim time | 1200 tu (11.2 days) | Covers outbound trip, flyby, and return |
| Camera | (0, 90.5, 0) | Above orbital plane, covers full Earth-Moon system |

### 8.1 Why Not Include the Sun?

The Sun's gravity is approximately uniform over the Earth-Moon distance scale (tidal acceleration is ~0.5% of Earth's gravity at the Moon's distance). Over the 6-day timescale of a free-return trajectory, solar perturbations shift the trajectory by <100 km — negligible compared to the sensitivity of the three-body problem itself. The real Apollo guidance system accounted for solar and planetary perturbations, but this simplified three-body simulation demonstrates the free-return mechanism clearly without them.

## 9. Sensitivity and Tuning

The free-return trajectory sits at the intersection of several qualitatively different outcomes in the three-body problem:

| Problem | Adjustment |
|---------|-----------|
| Spacecraft hits Moon | Increase `flyby_offset_deg` by 0.05° |
| Spacecraft misses Moon entirely | Decrease `flyby_offset_deg` by 0.05° |
| Spacecraft orbits Moon (captured) | Increase `v_tli` by ~5 m/s |
| Spacecraft doesn't reach Moon | Increase `v_tli` by ~20 m/s |
| Perilune too low (<100 km) | Increase `flyby_offset_deg` by 0.03° |
| Perilune too high (>1000 km) | Decrease `flyby_offset_deg` by 0.1° |
| No return to Earth after flyby | Try larger changes to `flyby_offset_deg` (±0.3°) |

The `flyby_offset_deg` parameter is the primary tuning knob. In the real Apollo program, mid-course corrections of 5-20 m/s were applied during transit to compensate for exactly this kind of sensitivity. Our simulation uses no corrections — the trajectory must be correct from TLI.

## 10. Comparison to Real Apollo Missions

| Parameter | This simulation | Apollo 11 | Apollo 13 |
|-----------|----------------|-----------|-----------|
| TLI velocity | 10,834 m/s | 10,834 m/s | 10,834 m/s |
| TLI altitude | 334 km | 334 km | 330 km |
| Transit time (to Moon) | ~72 hr (Kepler est.) | 75.3 hr | 76.0 hr |
| Perilune | ~200-400 km (estimated) | 115 km (lunar orbit, not free-return) | 254 km |
| Mid-course corrections | None (ballistic) | 4 (small) | 2 (critical) |

Apollo 11 was NOT on a free-return after TLI — it performed Lunar Orbit Insertion (LOI) to enter lunar orbit. Apollo 13 used the free-return trajectory for emergency recovery after the service module explosion, with perilune at 254 km — close to our target.

## 11. Limitations

1. **No Sun**: Solar gravitational tides are neglected (~0.5% perturbation over 6 days). Acceptable for a qualitative trajectory demonstration.

2. **Circular Moon orbit**: The real Moon's orbit is slightly eccentric (e = 0.055) and inclined 5.1° to the ecliptic. We use a circular, coplanar approximation.

3. **Point masses**: Earth and Moon are treated as point masses. Atmospheric drag (<100 km altitude), Earth's J2 oblateness, and lunar mascons are not modeled.

4. **No mid-course corrections**: Real Apollo missions applied small burns (~5-20 m/s) during transit. Our simulation is purely ballistic after TLI.

5. **2D trajectory**: The orbit is confined to the x-z plane. Real Apollo trajectories had out-of-plane components due to the Moon's orbital inclination and launch site latitude.

6. **Two-body analytical estimates**: The encounter analysis (§6) uses patched-conic approximations. The actual N-body simulation handles the full three-body gravitational field continuously, so actual perilune and deflection may differ from the analytical estimates.

## 12. Sources

| Source | Used for |
|--------|----------|
| Apollo 11 Mission Report (NASA SP-238, 1971) | TLI burnout state (velocity, altitude) |
| Apollo 13 Mission Report (NASA MSC-02680, 1970) | Free-return perilune (254 km), total trip time |
| IAU 2009 System of Astronomical Constants | Earth mass, Moon mass, Moon radius |
| Lunar Laser Ranging (Williams et al. 2014) | Earth-Moon distance (384,400 km) |
| Battin, R.H. "Introduction to Astrodynamics" (MIT, 1987) | Kepler's equation, patched conics, three-body problem |
| Bate, Mueller & White "Fundamentals of Astrodynamics" (Dover, 1971) | Transfer orbit calculations, sphere of influence |
