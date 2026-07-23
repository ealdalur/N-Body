# compute_oppositions.py

## Goal

Compute tables of Jupiter and Saturn opposition dates from J2000.0 (t=0) through 30 years (t=30), reporting both the simulation time and the corresponding calendar date. Opposition is when an outer planet is directly opposite the Sun as seen from Earth, making it the brightest and closest it will be during that synodic cycle. These tables serve as validation targets for the solar system simulation.

## Approach

1. Propagate orbital elements forward in time using mean motion (Keplerian two-body).
2. At each timestep, compute the heliocentric ecliptic longitude of Earth and the target planet.
3. Detect sign changes in the longitude difference (indicating the planet and Earth have the same heliocentric longitude).
4. Refine each crossing to high precision using bisection.
5. Convert the simulation time to a calendar date.

## Source Data

### Orbital Elements

The same JPL J2000.0 elements used in `compute_solar_system.py`, but only for Earth, Jupiter, and Saturn:

| Planet | a (AU) | e | i (deg) | Omega (deg) | omega (deg) | M0 (deg) | Mass (kg) |
|--------|--------|------|---------|-------------|-------------|-----------|-----------|
| Earth | 1.00000261 | 0.01671 | 0.00002 | 0.0 | 102.938 | 357.527 | 5.9724e24 |
| Jupiter | 5.20289 | 0.04839 | 1.304 | 100.474 | 274.255 | 19.668 | 1.8982e27 |
| Saturn | 9.53668 | 0.05386 | 2.486 | 113.662 | 338.936 | 317.355 | 5.6834e26 |

### Unit System

- Distance: AU
- Time: Julian years (365.25 days)
- Mass: solar masses
- G = 4*pi^2 (ensures P_Earth = 1 year)

## Equations and Calculations

### Mean Anomaly Propagation

For any time t (years since J2000.0), the mean anomaly is:

```
M(t) = M_0 + 360 * t / P
```

where the orbital period P is:

```
P = 2*pi * sqrt(a^3 / mu)
```

and mu = G*(1 + m_planet/M_sun). This gives P_Jupiter ~ 11.86 yr and P_Saturn ~ 29.46 yr.

### Heliocentric Ecliptic Longitude

At each time t:
1. Compute M(t) for the planet.
2. Solve Kepler's equation for E.
3. Compute position in the orbital plane (x', y').
4. Rotate to ecliptic frame (x, y) using Omega, omega, i.
5. Compute longitude: lambda = atan2(y, x).

Only the x and y ecliptic coordinates are needed (not z) since longitude is measured in the ecliptic plane.

### Opposition Detection

Opposition occurs when the planet has the same heliocentric longitude as Earth:

```
lambda_planet - lambda_earth = 0 (mod 2*pi)
```

The script computes the wrapped difference:

```
diff = (lambda_planet - lambda_earth + pi) mod 2*pi - pi
```

This maps the difference to [-pi, pi]. A zero crossing where |diff| < pi/2 corresponds to opposition (planet on the same side as Earth). The pi/2 guard rejects the conjunction case (planet behind the Sun) where the difference also passes through zero but as a discontinuity from +pi to -pi.

### Bisection Refinement

When a sign change is detected between consecutive timesteps, bisection narrows the interval:

```
for 50 iterations:
    t_mid = (t_a + t_b) / 2
    if sign(diff(t_mid)) == sign(diff(t_a)):
        t_a = t_mid
    else:
        t_b = t_mid
```

50 iterations of bisection gives precision of dt_coarse / 2^50, which is approximately 10^-18 years (sub-nanosecond).

### Calendar Date Conversion

Simulation time t (years from J2000.0) is converted to a calendar date:

1. Compute total days: `days = t * 365.25`
2. Add 0.5 (J2000.0 is noon Jan 1, 2000)
3. Subtract full years (accounting for leap years)
4. Subtract full months from the remaining days
5. Format as "YYYY Mon DD"

## Limitations and Assumptions

1. **Keplerian propagation only**: The planets follow fixed unperturbed orbits. Real Jupiter and Saturn have orbital perturbations (especially the Jupiter-Saturn 2:5 near-resonance) that can shift opposition dates by several days. Accuracy is typically within 3-5 days of actual opposition over 30 years.

2. **Static orbital elements**: No secular drift in a, e, i, Omega, omega. Over 30 years this introduces small errors but not enough to miss an opposition or produce a spurious one.

3. **Ecliptic longitude only**: True opposition is defined in right ascension (equatorial frame). Using ecliptic longitude introduces errors of up to ~1 day due to non-zero inclinations. Jupiter and Saturn have small inclinations (1.3 and 2.5 deg) so this is minor.

4. **No light-time correction**: The geometric alignment is computed, not the apparent position as seen from Earth.

5. **Julian year assumption**: The calendar conversion uses 365.25 days/year uniformly. The Gregorian calendar's actual irregular year lengths mean the day-of-month could be off by +/- 1 day in some cases.

6. **Two-body periods**: The computed periods assume no perturbations. Jupiter's actual synodic period varies by about 1 month over decades due to Saturn's influence.

## Code Implementation

### Structure

- **`solve_kepler()`**: Same Newton-Raphson solver as in `compute_solar_system.py`.
- **`get_ecliptic_xy()`**: Computes only the x, y ecliptic position (z is not needed for longitude).
- **`ecliptic_lon()`**: Wraps the full pipeline: propagate M(t), solve Kepler, rotate to ecliptic, return atan2(y, x).
- **`find_oppositions()`**: Coarse scan with bisection refinement.
- **`t_to_date()`**: Simulation-time-to-calendar converter.

### `ecliptic_lon(a, e, inc, Omega, omega, M0, mass_kg, t)`

Propagates the mean anomaly forward by t years using the Keplerian period, then returns the ecliptic longitude. The period is computed each call (cheap, avoids global state).

### `find_oppositions(planet_name, t_start, t_end, dt_coarse=0.001)`

Scans from t_start to t_end in steps of dt_coarse (0.001 years ~ 8.8 hours). At each step, computes the longitude difference and checks for sign reversals. The coarse step size of 0.001 years ensures no opposition can be missed (the fastest longitude change is Earth's at ~1 deg/day, and 8.8 hours corresponds to ~0.37 deg, well below the pi/2 guard window).

The bisection uses 50 iterations, which is overkill for the precision needed but costs negligible computation time for so few events.

### `t_to_date(t)`

Handles leap years via the Gregorian rule (divisible by 4, except centuries unless divisible by 400). Accumulates days year-by-year, then month-by-month. An edge-case guard prevents month overflow if floating-point accumulation places the day exactly at a year boundary.

### Output

Prints two tables (Jupiter and Saturn) with simulation time to 4 decimal places and calendar date in human-readable format. The synodic period is visible in the spacing: ~1.09 years for Jupiter, ~1.04 years for Saturn.
