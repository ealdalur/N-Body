"""
Compute Cartesian state vectors for all 8 planets at J2000.0 epoch.

J2000.0 = January 1, 2000, 12:00 TT (JD 2451545.0)

Orbital elements from JPL "Approximate Positions of the Major Planets"
(Standish, 1992). Output in AU / solar masses / years units where
G = 4*pi^2 = 39.4784176044.

Ecliptic plane mapped to x-z; y is out-of-plane (normal to ecliptic).
"""

import math

# Physical constants
G_SI = 6.67430e-11       # m^3 kg^-1 s^-2
AU = 1.495978707e11      # m (IAU 2012 exact)
M_sun_kg = 1.98892e30    # kg
yr = 365.25 * 86400      # seconds in a Julian year

# G in AU, solar mass, year units
G_code = 4.0 * math.pi**2

# Planet data: (name, a[AU], e, i[deg], Omega[deg], omega[deg], M[deg], mass[kg])
planets = [
    ("Mercury",  0.38709927, 0.20563593,  7.00497902,  48.33076593,  29.12703035, 174.79252722, 3.3011e23),
    ("Venus",    0.72333566, 0.00677672,  3.39467605,  76.67984255,  54.92262463,  50.37663232, 4.8675e24),
    ("Earth",    1.00000261, 0.01671123,  0.00001531,   0.0,        102.93768193, 357.52689973, 5.9724e24),
    ("Mars",     1.52371034, 0.09339410,  1.84969142,  49.55953891, 286.49683150,  19.39019754, 6.4171e23),
    ("Jupiter",  5.20288700, 0.04838624,  1.30439695, 100.47390909, 274.25457074,  19.66796068, 1.8982e27),
    ("Saturn",   9.53667594, 0.05386179,  2.48599187, 113.66242448, 338.93645383, 317.35536592, 5.6834e26),
    ("Uranus",  19.18916464, 0.04725744,  0.77263783,  74.01692503,  96.93735127, 142.28382821, 8.6810e25),
    ("Neptune", 30.06992276, 0.00859048,  1.77004347, 131.78422574, 273.18053653, 259.91520804, 1.02413e26),
]


def solve_kepler(M_rad, e, tol=1e-12, max_iter=100):
    """Solve Kepler's equation E - e*sin(E) = M via Newton-Raphson."""
    E = M_rad
    for _ in range(max_iter):
        dE = (E - e * math.sin(E) - M_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def orbital_to_cartesian(a, e, i_deg, Omega_deg, omega_deg, M_deg, mu):
    """
    Convert Keplerian elements to Cartesian state vector (ecliptic frame).

    Parameters:
        a       - semi-major axis (AU)
        e       - eccentricity
        i_deg   - inclination (degrees)
        Omega_deg - longitude of ascending node (degrees)
        omega_deg - argument of perihelion (degrees)
        M_deg   - mean anomaly (degrees)
        mu      - gravitational parameter G*(M_sun + M_planet)

    Returns:
        (x, y, z, vx, vy, vz) in ecliptic coordinates
    """
    i = math.radians(i_deg)
    Omega = math.radians(Omega_deg)
    omega = math.radians(omega_deg)
    M_rad = math.radians(M_deg)

    E = solve_kepler(M_rad, e)

    # Position in orbital plane
    xp = a * (math.cos(E) - e)
    yp = a * math.sqrt(1.0 - e * e) * math.sin(E)

    # Velocity in orbital plane
    n = math.sqrt(mu / (a**3))
    denom = 1.0 - e * math.cos(E)
    vxp = -a * n * math.sin(E) / denom
    vyp = a * n * math.sqrt(1.0 - e * e) * math.cos(E) / denom

    # Rotation to ecliptic frame
    cos_o = math.cos(omega); sin_o = math.sin(omega)
    cos_O = math.cos(Omega); sin_O = math.sin(Omega)
    cos_i = math.cos(i);     sin_i = math.sin(i)

    # Position
    x = (cos_O * cos_o - sin_O * sin_o * cos_i) * xp + (-cos_O * sin_o - sin_O * cos_o * cos_i) * yp
    y = (sin_O * cos_o + cos_O * sin_o * cos_i) * xp + (-sin_O * sin_o + cos_O * cos_o * cos_i) * yp
    z = (sin_o * sin_i) * xp + (cos_o * sin_i) * yp

    # Velocity
    vx = (cos_O * cos_o - sin_O * sin_o * cos_i) * vxp + (-cos_O * sin_o - sin_O * cos_o * cos_i) * vyp
    vy = (sin_O * cos_o + cos_O * sin_o * cos_i) * vxp + (-sin_O * sin_o + cos_O * cos_o * cos_i) * vyp
    vz = (sin_o * sin_i) * vxp + (cos_o * sin_i) * vyp

    return (x, y, z, vx, vy, vz)


if __name__ == "__main__":
    print("# Solar System state vectors at J2000.0 (2000-Jan-1.5 TT)")
    print(f"# Units: AU, solar masses, years. G = 4*pi^2 = {G_code:.10f}")
    print("# Ecliptic mapped to x-z plane (y = out-of-plane)")
    print()
    print(f"{'Planet':10s}  {'mass (M_sun)':>14s}  {'x (AU)':>16s} {'y (AU)':>16s} {'z (AU)':>16s}  {'vx (AU/yr)':>14s} {'vy (AU/yr)':>14s} {'vz (AU/yr)':>14s}")
    print("-" * 130)

    for name, a, e, inc, Omega, omega, M_anom, mass_kg in planets:
        m = mass_kg / M_sun_kg
        mu = G_code * (1.0 + m)
        x, y, z, vx, vy, vz = orbital_to_cartesian(a, e, inc, Omega, omega, M_anom, mu)

        # Rotate ecliptic to sim: ecl_x -> sim_x, ecl_z -> sim_y, -ecl_y -> sim_z
        # This is Rx(-90deg), a proper rotation preserving handedness.
        sx, sy, sz = x, z, -y
        svx, svy, svz = vx, vz, -vy

        print(f"{name:10s}  {m:>14.10e}  {sx:>+16.10e} {sy:>+16.10e} {sz:>+16.10e}  {svx:>+14.10e} {svy:>+14.10e} {svz:>+14.10e}")

    print()
    print("# To use in Solar_System.sim:")
    print("# Body  system  posX posY posZ  velX velY velZ  mass")
