"""
Compute Jupiter and Saturn oppositions from J2000.0 (t=0) to t=30 years.

Opposition occurs when an outer planet has the same heliocentric ecliptic
longitude as Earth (Sun-Earth-Planet aligned, planet opposite Sun in sky).

Uses the same Keplerian orbital elements as the N-Body simulation.
Output: simulation time (years since J2000.0) and calendar date.
"""

import math

M_sun_kg = 1.98892e30
G_code = 4.0 * math.pi**2

planets = {
    'Earth':   (1.00000261, 0.01671123, 0.00001531,   0.0,        102.93768193, 357.52689973, 5.9724e24),
    'Jupiter': (5.20288700, 0.04838624, 1.30439695, 100.47390909, 274.25457074,  19.66796068, 1.8982e27),
    'Saturn':  (9.53667594, 0.05386179, 2.48599187, 113.66242448, 338.93645383, 317.35536592, 5.6834e26),
}


def solve_kepler(M_rad, e, tol=1e-12):
    """Solve Kepler's equation E - e*sin(E) = M via Newton-Raphson."""
    E = M_rad
    for _ in range(100):
        dE = (E - e * math.sin(E) - M_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def get_ecliptic_xy(a, e, inc_deg, Omega_deg, omega_deg, M_deg, mu):
    """Get x,y position in ecliptic plane from orbital elements."""
    i = math.radians(inc_deg)
    Om = math.radians(Omega_deg)
    w = math.radians(omega_deg)
    M_rad = math.radians(M_deg)
    E = solve_kepler(M_rad, e)
    xp = a * (math.cos(E) - e)
    yp = a * math.sqrt(1.0 - e * e) * math.sin(E)
    cos_o = math.cos(w); sin_o = math.sin(w)
    cos_O = math.cos(Om); sin_O = math.sin(Om)
    cos_i = math.cos(i)
    x = (cos_O * cos_o - sin_O * sin_o * cos_i) * xp + (-cos_O * sin_o - sin_O * cos_o * cos_i) * yp
    y = (sin_O * cos_o + cos_O * sin_o * cos_i) * xp + (-sin_O * sin_o + cos_O * cos_o * cos_i) * yp
    return x, y


def ecliptic_lon(a, e, inc, Omega, omega, M0, mass_kg, t):
    """Compute heliocentric ecliptic longitude at time t (years from J2000.0)."""
    m = mass_kg / M_sun_kg
    mu = G_code * (1.0 + m)
    P = 2 * math.pi * math.sqrt(a**3 / mu)
    M = M0 + 360.0 * t / P
    M = M % 360.0
    x, y = get_ecliptic_xy(a, e, inc, Omega, omega, M, mu)
    return math.atan2(y, x)


def find_oppositions(planet_name, t_start, t_end, dt_coarse=0.001):
    """Find times when planet is at opposition (same heliocentric longitude as Earth)."""
    p = planets[planet_name]
    e_data = planets['Earth']
    oppositions = []
    prev_diff = None
    t = t_start
    while t < t_end:
        lon_e = ecliptic_lon(*e_data, t)
        lon_p = ecliptic_lon(*p, t)
        diff = lon_p - lon_e
        diff = (diff + math.pi) % (2 * math.pi) - math.pi
        if prev_diff is not None:
            if prev_diff * diff < 0 and abs(diff) < math.pi / 2:
                ta, tb = t - dt_coarse, t
                da = prev_diff
                for _ in range(50):
                    tm = (ta + tb) / 2
                    lon_em = ecliptic_lon(*e_data, tm)
                    lon_pm = ecliptic_lon(*p, tm)
                    dm = lon_pm - lon_em
                    dm = (dm + math.pi) % (2 * math.pi) - math.pi
                    if dm * da > 0:
                        ta = tm
                        da = dm
                    else:
                        tb = tm
                oppositions.append((ta + tb) / 2)
        prev_diff = diff
        t += dt_coarse
    return oppositions


def t_to_date(t):
    """Convert simulation time (years from J2000.0) to calendar date string."""
    # J2000.0 = 2000 Jan 1.5 (noon Jan 1)
    total_days = t * 365.25
    year = 2000
    day_accum = total_days + 0.5

    while True:
        days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
        if day_accum >= days_in_year:
            day_accum -= days_in_year
            year += 1
        else:
            break

    is_leap = (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0))
    months = [31, 29 if is_leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    m = 0
    while m < 12 and day_accum >= months[m]:
        day_accum -= months[m]
        m += 1
    if m >= 12:
        m = 11
        day_accum = 30
    return f'{year} {month_names[m]} {int(day_accum) + 1:2d}'


if __name__ == "__main__":
    print('=== Jupiter Oppositions (t=0 to t=30) ===')
    print(f'{"Sim time (t)":>12s}   {"Calendar date":s}')
    print('-' * 32)
    for t in find_oppositions('Jupiter', 0.0, 30.0):
        print(f'{t:12.4f}   {t_to_date(t)}')

    print()
    print('=== Saturn Oppositions (t=0 to t=30) ===')
    print(f'{"Sim time (t)":>12s}   {"Calendar date":s}')
    print('-' * 32)
    for t in find_oppositions('Saturn', 0.0, 30.0):
        print(f'{t:12.4f}   {t_to_date(t)}')
