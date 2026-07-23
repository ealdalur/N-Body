# compute_solar_system.py

## Goal

Compute Cartesian state vectors (position and velocity) for all eight planets at the J2000.0 epoch (January 1, 2000, 12:00 TT), formatted for direct insertion into the N-Body simulation's `Solar_System.sim` script. The output is in the simulation's natural unit system (AU, solar masses, years) with the ecliptic plane mapped into the simulator's coordinate frame.

## Approach

1. Start from tabulated Keplerian orbital elements at J2000.0.
2. Solve Kepler's equation to find each planet's eccentric anomaly.
3. Compute position and velocity in the orbital plane.
4. Rotate from the orbital plane to the ecliptic frame using the standard three-rotation sequence (argument of perihelion, inclination, longitude of ascending node).
5. Apply a coordinate frame rotation to map the ecliptic frame into the simulator's convention (ecliptic x-z plane with y as the out-of-plane axis).

## Source Data

### Orbital Elements

Sourced from JPL's "Approximate Positions of the Major Planets" (Standish, 1992). Each planet has six Keplerian elements at J2000.0:

| Symbol | Name | Units | Meaning |
|--------|------|-------|---------|
| a | Semi-major axis | AU | Size of the orbit ellipse |
| e | Eccentricity | dimensionless | Shape of the ellipse (0 = circle, 1 = parabola) |
| i | Inclination | degrees | Tilt of orbital plane relative to ecliptic |
| Omega | Longitude of ascending node | degrees | Where the orbit crosses the ecliptic going north |
| omega | Argument of perihelion | degrees | Angle from ascending node to closest approach point |
| M | Mean anomaly | degrees | Parametrizes position along orbit at epoch |

### Planet Masses

Masses in kg from IAU-recommended values, converted to solar masses using M_sun = 1.98892e30 kg.

### Physical Constants

| Constant | Value | Source |
|----------|-------|--------|
| G (SI) | 6.67430e-11 m^3 kg^-1 s^-2 | CODATA 2018 |
| AU | 1.495978707e11 m | IAU 2012 (exact by definition) |
| M_sun | 1.98892e30 kg | IAU nominal |
| Julian year | 365.25 * 86400 s | IAU definition |
| G (code) | 4*pi^2 = 39.4784176044 | Derived: ensures Earth's period = 1 year at 1 AU |

## Equations and Calculations

### Kepler's Equation

The mean anomaly M advances linearly with time, but the actual angular position (true anomaly) does not. Kepler's equation relates mean anomaly to eccentric anomaly E:

```
E - e*sin(E) = M
```

This transcendental equation is solved iteratively via Newton-Raphson:

```
E_{n+1} = E_n - (E_n - e*sin(E_n) - M) / (1 - e*cos(E_n))
```

Convergence is typically achieved in 4-6 iterations for planetary eccentricities (e < 0.21). Tolerance is set to 1e-12 radians.

### Position in the Orbital Plane

Once E is known, the position in the orbital plane (with x-axis toward perihelion) is:

```
x' = a * (cos(E) - e)
y' = a * sqrt(1 - e^2) * sin(E)
```

### Velocity in the Orbital Plane

The mean motion n and the velocity components are:

```
n = sqrt(mu / a^3)
vx' = -a * n * sin(E) / (1 - e*cos(E))
vy' = a * n * sqrt(1 - e^2) * cos(E) / (1 - e*cos(E))
```

where mu = G*(M_sun + M_planet). Including the planet mass in mu gives the correct two-body period.

### Rotation to Ecliptic Frame

The orbital plane is rotated to the ecliptic reference frame by applying three successive rotations:

1. R_z(omega) - rotate by argument of perihelion about the orbit normal
2. R_x(i) - rotate by inclination about the node line
3. R_z(Omega) - rotate by longitude of ascending node about the ecliptic pole

The combined rotation matrix elements are applied directly to both position and velocity vectors. The explicit formulas avoid constructing intermediate 3x3 matrices:

```
x = (cos(O)*cos(o) - sin(O)*sin(o)*cos(i)) * x' + (-cos(O)*sin(o) - sin(O)*cos(o)*cos(i)) * y'
y = (sin(O)*cos(o) + cos(O)*sin(o)*cos(i)) * x' + (-sin(O)*sin(o) + cos(O)*cos(o)*cos(i)) * y'
z = sin(o)*sin(i) * x' + cos(o)*sin(i) * y'
```

(Same for velocity with vx', vy'.)

### Frame Rotation to Simulator Coordinates

The simulator uses y as the out-of-ecliptic axis (so the ecliptic disc lies in x-z). The standard ecliptic frame has z as the pole. The mapping is a rotation Rx(-90 deg):

```
sim_x = ecl_x
sim_y = ecl_z
sim_z = -ecl_y
```

This is a proper rotation (preserves handedness) and is applied to both position and velocity.

## Limitations and Assumptions

1. **Static elements**: The orbital elements are valid at J2000.0 only. The script does not account for secular perturbations (rate-of-change terms in a, e, i, etc.). For positions within a few centuries of J2000.0, propagation via the N-Body simulation itself is more accurate than extrapolating elements.

2. **Two-body problem**: Each planet's orbit is computed independently as a Keplerian orbit around the Sun. Mutual planetary perturbations are not included in the initial conditions (they emerge naturally during N-Body integration).

3. **Point masses**: No oblateness (J2), no relativistic precession. Mercury's perihelion advance (43 arcsec/century from GR) is absent.

4. **No Moon**: The Moon's initial conditions are computed separately by `compute_moon.py` (documented in [compute_moon.md](compute_moon.md)), which uses the same Keplerian formalism but with Earth as the central body and adds Earth's heliocentric state to produce the Moon's heliocentric position and velocity.

5. **Ecliptic-of-date vs J2000 ecliptic**: The elements are referenced to the J2000 ecliptic, which differs from the ecliptic of any other date by precession. This is correct for the simulation epoch.

6. **Mass precision**: Planet masses are taken to 4-5 significant figures. The Sun's mass dominates mu, so planet mass errors have negligible effect on orbital mechanics.

## Code Implementation

### Structure

The script is organized into:
- **Constants block**: Physical constants and unit conversions.
- **Data table**: A list of tuples, one per planet, containing all six orbital elements plus mass.
- **`solve_kepler()`**: Newton-Raphson Kepler solver.
- **`orbital_to_cartesian()`**: Full pipeline from elements to ecliptic Cartesian state.
- **Main block**: Iterates over planets, applies the frame rotation, and prints formatted output.

### `solve_kepler(M_rad, e, tol=1e-12, max_iter=100)`

Starts with E = M as the initial guess (good for low eccentricity). Newton-Raphson converges quadratically. The tolerance of 1e-12 radians corresponds to sub-meter positional accuracy at planetary distances.

### `orbital_to_cartesian(a, e, i_deg, Omega_deg, omega_deg, M_deg, mu)`

Converts degrees to radians, solves Kepler, computes orbital-plane position/velocity, then applies the full rotation to ecliptic coordinates. Returns a 6-element tuple (x, y, z, vx, vy, vz).

The gravitational parameter `mu = G*(1 + m_planet)` includes the planet mass in solar-mass units. This gives the correct mean motion n = sqrt(mu/a^3) for the two-body problem.

### Output Formatting

The main block prints a table with 10-digit scientific notation for positions and velocities, suitable for direct copy into a .sim file's Body lines. The frame rotation (ecliptic to simulator) is applied inline before printing.
