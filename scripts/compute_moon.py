"""
Compute the Moon's Cartesian state vector at J2000.0 epoch.

J2000.0 = January 1, 2000, 12:00 TT (JD 2451545.0)

The Moon's orbital elements are referenced to the ecliptic plane and
use the same Keplerian formalism as the planets, but with Earth as
the central body. The result is added to Earth's heliocentric position
to give the Moon's heliocentric state vector.

Lunar orbital elements from Meeus "Astronomical Algorithms" (2nd ed.),
based on Chapront-Touze & Chapront (1988) ELP 2000-82 theory.

Output in AU / solar masses / years units where G = 4*pi^2.
Ecliptic plane mapped to x-z; y is out-of-plane (sim frame).
"""

import math

# === Physical constants ===
AU_km = 1.495978707e8       # km per AU (IAU 2012)
M_sun_kg = 1.98892e30       # kg
M_earth_kg = 5.9724e24      # kg
M_moon_kg = 7.342e22        # kg (IAU)

# G in AU, solar mass, year units
G_code = 4.0 * math.pi**2

# === Moon orbital elements at J2000.0 ===
# Source: Meeus "Astronomical Algorithms" Ch. 47, Table 47.A
# Based on ELP 2000-82 (Chapront-Touze & Chapront 1988)
a_km = 384400.0             # semi-major axis (km)
e = 0.0549006               # eccentricity
i_deg = 5.145               # inclination to ecliptic (deg)
Omega_deg = 125.08          # longitude of ascending node (deg)

# Mean anomaly and argument of perigee at J2000.0
# L (mean longitude) = 218.32 deg
# M' (mean anomaly) = 134.96 deg
# omega_bar (longitude of perigee) = L - M' = 83.36 deg
# omega (argument of perigee) = omega_bar - Omega = 83.36 - 125.08 = -41.72 = 318.28 deg
M_deg = 134.96              # mean anomaly
omega_deg = 318.15          # argument of perigee

# === Earth's heliocentric state vector (from compute_solar_system.py) ===
# Already in sim frame (ecliptic mapped to x-z)
earth_pos = (-1.7717142365e-01, 2.5844928543e-07, -9.6721445286e-01)
earth_vel = (-6.2835665286e+00, -3.0883283534e-07, 1.1557686506e+00)


def solve_kepler(M_rad, e, tol=1e-12, max_iter=100):
    """Solve Kepler's equation E - e*sin(E) = M via Newton-Raphson."""
    E = M_rad
    for _ in range(max_iter):
        dE = (E - e * math.sin(E) - M_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def moon_state_vector():
    """
    Compute Moon's position and velocity relative to Earth in ecliptic frame,
    then convert to heliocentric sim frame.

    Returns:
        (x, y, z, vx, vy, vz, mass) in sim coordinates (AU, AU/yr, solar masses)
    """
    # Convert to AU
    a = a_km / AU_km

    # Gravitational parameter: Moon orbits Earth
    # mu = G * (M_earth + M_moon) in code units (solar masses)
    mu = G_code * (M_earth_kg + M_moon_kg) / M_sun_kg

    # Solve Kepler's equation
    M_rad = math.radians(M_deg)
    E = solve_kepler(M_rad, e)

    # Position in orbital plane (perigee along x-axis)
    xp = a * (math.cos(E) - e)
    yp = a * math.sqrt(1.0 - e * e) * math.sin(E)

    # Radial distance (for verification)
    r = a * (1.0 - e * math.cos(E))

    # Velocity in orbital plane
    n = math.sqrt(mu / a**3)    # mean motion (rad/yr)
    denom = 1.0 - e * math.cos(E)
    vxp = -a * n * math.sin(E) / denom
    vyp = a * n * math.sqrt(1.0 - e * e) * math.cos(E) / denom

    # Rotation matrices: orbital plane -> ecliptic frame
    inc = math.radians(i_deg)
    Om = math.radians(Omega_deg)
    w = math.radians(omega_deg)

    cos_o = math.cos(w);  sin_o = math.sin(w)
    cos_O = math.cos(Om); sin_O = math.sin(Om)
    cos_i = math.cos(inc); sin_i = math.sin(inc)

    # Position in ecliptic frame (geocentric)
    x_ecl = (cos_O * cos_o - sin_O * sin_o * cos_i) * xp + \
             (-cos_O * sin_o - sin_O * cos_o * cos_i) * yp
    y_ecl = (sin_O * cos_o + cos_O * sin_o * cos_i) * xp + \
             (-sin_O * sin_o + cos_O * cos_o * cos_i) * yp
    z_ecl = (sin_o * sin_i) * xp + (cos_o * sin_i) * yp

    # Velocity in ecliptic frame (geocentric)
    vx_ecl = (cos_O * cos_o - sin_O * sin_o * cos_i) * vxp + \
              (-cos_O * sin_o - sin_O * cos_o * cos_i) * vyp
    vy_ecl = (sin_O * cos_o + cos_O * sin_o * cos_i) * vxp + \
              (-sin_O * sin_o + cos_O * cos_o * cos_i) * vyp
    vz_ecl = (sin_o * sin_i) * vxp + (cos_o * sin_i) * vyp

    # Frame rotation: ecliptic -> sim (Rx(-90 deg))
    # sim_x = ecl_x, sim_y = ecl_z, sim_z = -ecl_y
    dx_sim = x_ecl
    dy_sim = z_ecl
    dz_sim = -y_ecl

    dvx_sim = vx_ecl
    dvy_sim = vz_ecl
    dvz_sim = -vy_ecl

    # Heliocentric = Earth + geocentric offset
    x = earth_pos[0] + dx_sim
    y = earth_pos[1] + dy_sim
    z = earth_pos[2] + dz_sim

    vx = earth_vel[0] + dvx_sim
    vy = earth_vel[1] + dvy_sim
    vz = earth_vel[2] + dvz_sim

    mass = M_moon_kg / M_sun_kg

    return x, y, z, vx, vy, vz, mass, r, a, n


if __name__ == "__main__":
    x, y, z, vx, vy, vz, mass, r, a, n = moon_state_vector()

    print("=" * 70)
    print("MOON STATE VECTOR AT J2000.0")
    print("=" * 70)

    print(f"\n--- Orbital Elements ---")
    print(f"  Semi-major axis:         {a_km:.0f} km = {a_km/AU_km:.10e} AU")
    print(f"  Eccentricity:            {e}")
    print(f"  Inclination:             {i_deg} deg (to ecliptic)")
    print(f"  Long. ascending node:    {Omega_deg} deg")
    print(f"  Argument of perigee:     {omega_deg} deg")
    print(f"  Mean anomaly:            {M_deg} deg")

    print(f"\n--- Derived Quantities ---")
    P_days = 2 * math.pi / n * 365.25
    mu = G_code * (M_earth_kg + M_moon_kg) / M_sun_kg
    v_orbital_kms = math.sqrt(mu / a) * AU_km * 1e3 / (365.25 * 86400)  # rough
    print(f"  Orbital period:          {P_days:.2f} days (observed: 27.32 days)")
    print(f"  Distance at epoch:       {r*AU_km:.0f} km (a*(1-e*cos(E)))")
    print(f"  Mean motion:             {math.degrees(n):.4f} deg/yr")
    print(f"  Gravitational parameter: {mu:.10e} (G*(M_earth+M_moon)/M_sun, code units)")

    print(f"\n--- Geocentric Offset (sim frame) ---")
    dx = x - earth_pos[0]
    dy = y - earth_pos[1]
    dz = z - earth_pos[2]
    dist_km = math.sqrt(dx**2 + dy**2 + dz**2) * AU_km
    print(f"  dx = {dx:+.10e} AU")
    print(f"  dy = {dy:+.10e} AU")
    print(f"  dz = {dz:+.10e} AU")
    print(f"  |offset| = {dist_km:.0f} km")

    print(f"\n--- Heliocentric State Vector (sim frame) ---")
    print(f"  x  = {x:+.10e} AU")
    print(f"  y  = {y:+.10e} AU")
    print(f"  z  = {z:+.10e} AU")
    print(f"  vx = {vx:+.10e} AU/yr")
    print(f"  vy = {vy:+.10e} AU/yr")
    print(f"  vz = {vz:+.10e} AU/yr")
    print(f"  mass = {mass:.10e} solar masses")

    print(f"\n--- Script Line ---")
    print(f"# Moon")
    print(f"Body  0   {x:.10e} {y:.10e} {z:.10e}   {vx:.10e} {vy:.10e} {vz:.10e}   {mass:.10e}")

    print(f"\n--- Verification ---")
    print(f"  Moon mass / Earth mass = {M_moon_kg/M_earth_kg:.6f} (expected: 0.01230)")
    print(f"  Earth-Moon distance:   {dist_km:.0f} km (expected: ~400,000 km at M=135 deg)")
    print(f"  Orbital period:        {P_days:.2f} days (expected: 27.32 days)")
