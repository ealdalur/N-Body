# compute_moon.py

## Goal

Compute the Moon's heliocentric Cartesian state vector at J2000.0 for use in the `Solar_System.sim` script. Unlike the planets, which orbit the Sun, the Moon orbits Earth, so its state vector must be computed as a geocentric orbit and then transformed to heliocentric coordinates by adding Earth's position and velocity.

## Approach

1. Start from tabulated Keplerian orbital elements of the Moon at J2000.0, referenced to the ecliptic plane.
2. Solve Kepler's equation for the Moon's eccentric anomaly.
3. Compute geocentric position and velocity in the Moon's orbital plane.
4. Rotate from the orbital plane to the ecliptic frame using the standard three-rotation sequence.
5. Apply the ecliptic-to-simulator frame rotation (Rx(-90 deg)).
6. Add Earth's heliocentric state vector (already computed by `compute_solar_system.py`) to obtain the Moon's heliocentric state.

## Source Data

### Lunar Orbital Elements

Sourced from Meeus "Astronomical Algorithms" (2nd edition), based on the ELP 2000-82 lunar theory (Chapront-Touze & Chapront, 1988):

| Symbol | Name | Value | Units |
|--------|------|-------|-------|
| a | Semi-major axis | 384,400 | km |
| e | Eccentricity | 0.0549006 | dimensionless |
| i | Inclination (to ecliptic) | 5.145 | degrees |
| Omega | Longitude of ascending node | 125.08 | degrees |
| omega | Argument of perigee | 318.15 | degrees |
| M | Mean anomaly | 134.96 | degrees |

These are osculating elements at the J2000.0 epoch.

#### Derivation of omega

The argument of perigee is not directly tabulated in many references. It is derived from:
- Mean longitude: L = 218.32 deg
- Mean anomaly: M' = 134.96 deg
- Longitude of perigee: omega_bar = L - M' = 83.36 deg
- Argument of perigee: omega = omega_bar - Omega = 83.36 - 125.08 = -41.72 = 318.28 deg

The value 318.15 deg used in the script reflects a slightly different source precision for L.

### Physical Constants

| Constant | Value | Source |
|----------|-------|--------|
| AU | 1.495978707e8 km | IAU 2012 (exact) |
| M_sun | 1.98892e30 kg | IAU nominal |
| M_earth | 5.9724e24 kg | IAU 2015 |
| M_moon | 7.342e22 kg | IAU 2015 |
| G (code) | 4*pi^2 = 39.4784 | AU^3 / (Msun * yr^2) |

### Earth's State Vector

The Moon script depends on Earth's heliocentric state vector, which is computed independently by `compute_solar_system.py` and hardcoded:

```
position: (-0.1772, 0.0000, -0.9672) AU  [sim frame]
velocity: (-6.2836, 0.0000, 1.1558) AU/yr [sim frame]
```

## Equations and Calculations

### Semi-major Axis Conversion

The Moon's semi-major axis must be converted from km to AU:

```
a = 384400 km / 149597870.7 km/AU = 2.5696e-3 AU
```

This is about 1/389 of Earth's orbital radius, so the Moon appears very close to Earth at the scale of the solar system.

### Gravitational Parameter

The Moon orbits Earth, not the Sun. The gravitational parameter for the Earth-Moon two-body problem, expressed in code units (AU, solar masses, years), is:

```
mu = G_code * (M_earth + M_moon) / M_sun
   = 39.4784 * (5.9724e24 + 7.342e22) / 1.98892e30
   = 1.2000e-4
```

This is much smaller than the Sun's mu (= G_code * 1.0 = 39.4784), reflecting that the Moon's orbital velocity is governed by Earth's mass, not the Sun's.

### Kepler's Equation

Identical to the planetary case:

```
E - e*sin(E) = M
```

Solved by Newton-Raphson. With e = 0.0549 (low eccentricity), convergence is rapid (~4 iterations).

### Position and Velocity in the Orbital Plane

```
x' = a * (cos(E) - e)
y' = a * sqrt(1 - e^2) * sin(E)
```

```
n = sqrt(mu / a^3)        [mean motion, rad/yr]
vx' = -a * n * sin(E) / (1 - e*cos(E))
vy' = a * n * sqrt(1 - e^2) * cos(E) / (1 - e*cos(E))
```

The mean motion n corresponds to:
```
P = 2*pi / n = 0.0747 yr = 27.29 days
```

This is very close to the observed sidereal period of 27.32 days.

### Rotation to Ecliptic Frame

The same three-rotation sequence as for planets (omega, i, Omega). The Moon's orbital plane is inclined 5.145 deg to the ecliptic (not to Earth's equator), so the rotation is directly to ecliptic coordinates.

```
x_ecl = R_z(Omega) * R_x(i) * R_z(omega) * x'
```

The result is the Moon's position relative to Earth in the ecliptic frame.

### Frame Rotation (Ecliptic to Simulator)

Same mapping as for planets:

```
sim_x = ecl_x
sim_y = ecl_z
sim_z = -ecl_y
```

### Heliocentric Conversion

The final step adds Earth's heliocentric state:

```
Moon_helio = Earth_helio + Moon_geo
```

This is valid because the simulator uses heliocentric positions for all bodies (including the Moon), not Earth-centered coordinates. The N-body integrator naturally handles the three-body Sun-Earth-Moon dynamics from these initial conditions.

## Limitations and Assumptions

1. **Two-body Keplerian orbit**: The Moon's initial conditions are computed assuming a pure Keplerian orbit around Earth. Solar perturbations (which produce the 18.6-year nodal precession and the evection) are not included in the initial state. The N-body simulation will naturally generate these perturbations during integration.

2. **Osculating vs. mean elements**: The elements used are mean elements (averaged over short-period perturbations). Osculating elements at J2000.0 would give a slightly different position (~100 km difference). For the simulation's purposes, this is negligible compared to the softening length.

3. **Static node and perigee**: The Moon's node regresses at ~19.35 deg/year and its perigee advances at ~40.7 deg/year. The elements used are snapshots at J2000.0; the simulation's own dynamics will drive nodal and apsidal motion.

4. **Earth's state vector dependency**: The Moon's heliocentric position inherits any errors in Earth's computed state vector from `compute_solar_system.py`. Since Earth's elements are accurate to ~10 km, this contribution is negligible.

5. **No libration or tidal effects**: The Moon is treated as a point mass. No tidal deformation, no spin-orbit coupling, and no libration physics.

6. **Simulation timestep adequacy**: With dt = 0.0001 years = 0.88 hours, the Moon's 27.3-day orbit is resolved by ~750 timesteps per orbit. This is adequate for stable leapfrog integration at the Moon's eccentricity (0.055), though fine orbital details (apsidal precession rate) may accumulate small phase errors over many orbits.

7. **Mass precision**: The Moon's mass (7.342e22 kg) is known to better than 0.1%. The Earth/Moon mass ratio of 81.3 is well-established. Neither contributes meaningful error.

8. **Element source precision**: Different references give slightly different values for omega and Omega at J2000.0 (varying by ~0.1-0.2 deg). This translates to ~700 km positional uncertainty, which is below the simulation's resolution threshold.

## Code Implementation

### Structure

- **Constants block**: AU conversion, masses, G in code units.
- **Elements block**: Moon's six Keplerian elements at J2000.0.
- **Earth state**: Hardcoded from `compute_solar_system.py` output.
- **`solve_kepler()`**: Same Newton-Raphson solver as the planetary script.
- **`moon_state_vector()`**: Full pipeline from elements to heliocentric sim-frame state.
- **Main block**: Prints diagnostic information and the formatted Body line.

### `moon_state_vector()`

This function encapsulates the complete computation:

1. Converts semi-major axis from km to AU.
2. Computes mu using Earth+Moon mass (not Sun mass).
3. Solves Kepler, computes orbital plane position/velocity.
4. Applies the three-rotation sequence to ecliptic.
5. Applies the ecliptic-to-sim frame rotation.
6. Adds Earth's state vector.
7. Returns the full state plus diagnostic values (r, a, n).

The key difference from the planetary script is step 2 (mu based on Earth mass) and step 6 (geocentric-to-heliocentric addition).

### Verification Output

The main block prints several cross-checks:
- **Orbital period**: Should match 27.32 days (computed: 27.29 days, within 0.1%).
- **Earth-Moon distance**: Should be ~400,000 km at M=135 deg (perigee passage past, approaching apogee). Computed: 399,860 km.
- **Mass ratio**: M_moon/M_earth = 0.01229, matching the known value.

### Output Formatting

The script outputs a ready-to-use Body line for the .sim file with 10-digit scientific notation, matching the format of the planetary entries produced by `compute_solar_system.py`.
