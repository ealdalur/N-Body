"""
Compute initial conditions for an Apollo-style free-return trajectory.

Three-body system: Earth + Moon + spacecraft (test particle).
The spacecraft starts at the instant of Trans-Lunar Injection (TLI)
burnout and coasts ballistically to the Moon, swings around the far
side via gravitational slingshot, and returns to Earth.

=== Unit system ===

Chosen so that G=1 with Earth-radius and Earth-mass as base units:

  1 distance unit (du) = 6,371 km = 1 Earth radius
  1 mass unit (mu)     = 5.972e24 kg = 1 Earth mass
  1 time unit (tu)     = sqrt(R_e^3 / (G * M_e)) = 805.5 s = 13.4 min
  1 velocity unit (vu) = R_e / tu = 7,909.5 m/s = 7.91 km/s

Verification:
  - Circular orbit speed at r=1 du: v = sqrt(G*M/r) = sqrt(1) = 1.0 vu = 7.91 km/s  (correct)
  - Escape velocity at r=1 du: v = sqrt(2) = 1.414 vu = 11.18 km/s  (correct)
  - Moon orbital distance: 384,400 / 6,371 = 60.3 du  (matches "60 Earth radii")

=== Coordinate convention ===

The orbital plane is x-z. The +y axis points out of the plane (north pole).
Counter-clockwise (CCW) motion viewed from +y means:

  omega = (0, +omega, 0)
  position at angle theta CCW from +x: (r*cos(theta), 0, -r*sin(theta))
  velocity (tangent) at angle theta:   v * (-sin(theta), 0, -cos(theta))

This follows from omega x r in a right-handed coordinate system.
At theta=0 (i.e. along +x): velocity is in the -z direction.

=== Physical data sources ===

Earth & Moon:
  - Mass ratio: M_moon/M_earth = 1/81.3 (IAU 2009)
  - Earth-Moon distance: 384,400 km (Lunar Laser Ranging, Williams et al. 2014)
  - Moon radius: 1,737 km (IAU)
  - Moon orbital period: 27.32 days (sidereal)

Apollo TLI:
  - Burnout altitude: 334 km (Apollo 11 Mission Report, NASA SP-238)
  - Burnout velocity: 10,834 m/s inertial (Apollo 11 S-IVB cutoff)
  - Transit time to Moon: 72 hours typical (Apollo 8, 10, 11, 13)
  - Free-return perilune: 254 km (Apollo 13, NASA MSC-02680)

=== Algorithm ===

1. Compute Moon's circular orbital velocity from two-body mechanics.
2. Compute spacecraft transfer orbit elements (ellipse with TLI at perigee).
3. Find the true anomaly at which the spacecraft crosses the Moon's orbital
   distance (r_moon = 60.34 du). This is 172.2 deg, NOT 180 deg (apogee).
4. Compute the transit time from perigee to that true anomaly using
   Kepler's equation.
5. Determine how far the Moon moves during transit (39.3 deg CCW).
6. Aim the spacecraft's crossing point at the Moon's FUTURE position:
   perigee_angle = moon_future_angle - nu_at_crossing
7. Add a small flyby offset (0.7 deg) to ensure the spacecraft passes
   behind the Moon (far side) rather than colliding with it.
"""

import math

# ============================================================
# Unit system
# ============================================================
R_earth_km = 6371.0          # Earth mean radius (km)
M_earth_kg = 5.972e24        # Earth mass (kg)
G_SI = 6.674e-11             # gravitational constant (SI)

du = R_earth_km              # 1 code distance unit = 6,371 km
mu_kg = M_earth_kg           # 1 code mass unit = 5.972e24 kg
tu = math.sqrt((R_earth_km * 1e3)**3 / (G_SI * M_earth_kg))  # seconds
vu = R_earth_km * 1e3 / tu   # m/s per code velocity unit

print("=" * 65)
print("APOLLO FREE-RETURN TRAJECTORY PARAMETERS")
print("=" * 65)
print(f"\n--- Unit system (G=1) ---")
print(f"  1 du = {R_earth_km:.0f} km (Earth radius)")
print(f"  1 mu = {M_earth_kg:.3e} kg (Earth mass)")
print(f"  1 tu = {tu:.1f} s = {tu/60:.2f} min")
print(f"  1 vu = {vu:.1f} m/s = {vu/1000:.4f} km/s")
print(f"  G = 1.0 in these units")

# ============================================================
# Moon parameters
# ============================================================
M_moon_kg = 7.342e22         # Moon mass (IAU)
R_moon_km = 1737.0           # Moon mean radius (IAU)
r_moon_km = 384400.0         # Earth-Moon distance (LLR)

M_moon = M_moon_kg / M_earth_kg           # 0.01229 code mass units
R_moon = R_moon_km / du                    # 0.2726 code distance units
r_moon = r_moon_km / du                    # 60.34 code distance units
M_total = 1.0 + M_moon                    # Earth + Moon

# Moon orbital velocity (circular two-body problem)
# The relative velocity satisfies: v_rel = sqrt(G*(M1+M2)/r)
# Each body's velocity in the COM frame is proportional to the other's mass.
v_rel = math.sqrt(M_total / r_moon)
v_moon = v_rel * 1.0 / M_total            # Moon's speed in COM frame
v_earth = v_rel * M_moon / M_total        # Earth's speed (opposite direction)

T_orbit = 2 * math.pi * r_moon / v_rel    # orbital period in code units
omega_moon = 2 * math.pi / T_orbit        # angular velocity (rad/tu)

print(f"\n--- Moon ---")
print(f"  Mass:              {M_moon:.5f} mu ({M_moon_kg:.3e} kg)")
print(f"  Mass ratio E/M:    {1/M_moon:.1f}")
print(f"  Radius:            {R_moon:.4f} du ({R_moon_km:.0f} km)")
print(f"  Orbital distance:  {r_moon:.2f} du ({r_moon_km:.0f} km)")
print(f"  Orbital period:    {T_orbit:.1f} tu ({T_orbit*tu/86400:.2f} days)")
print(f"  Relative velocity: {v_rel:.5f} vu ({v_rel*vu/1000:.4f} km/s)")
print(f"  Moon velocity:     {v_moon:.5f} vu ({v_moon*vu/1000:.4f} km/s)")
print(f"  Earth velocity:    {v_earth:.6f} vu ({v_earth*vu/1000:.5f} km/s)")
print(f"  Angular velocity:  {omega_moon:.6f} rad/tu")
print(f"  Sphere of influence: {r_moon*(M_moon)**(2.0/5.0):.1f} du ({r_moon_km*(M_moon_kg/M_earth_kg)**(2.0/5.0):.0f} km)")

# ============================================================
# Trans-Lunar Injection (TLI) parameters
# ============================================================
# From Apollo 11 Mission Report (NASA SP-238, 1971):
#   S-IVB TLI cutoff conditions:
#     Altitude: 334 km above Earth surface
#     Inertial velocity: 10,834 m/s
#     Flight path angle: ~0 deg (tangential burnout)
#
# This is 99.4% of escape velocity at that altitude, placing the
# spacecraft on a highly eccentric ellipse (e=0.975) whose apogee
# extends well beyond the Moon's orbit.

alt_tli_km = 334.0           # TLI burnout altitude (Apollo 11)
v_tli_ms = 10834.0           # TLI burnout speed (Apollo 11)

r_tli = (R_earth_km + alt_tli_km) / du    # burnout radius in code units
v_tli = v_tli_ms / vu                     # burnout speed in code units

# Escape velocity at burnout altitude: v_esc = sqrt(2*G*M/r) = sqrt(2/r)
v_esc_tli = math.sqrt(2.0 / r_tli)

# Two-body (Earth-only) Keplerian orbit elements
E_specific = 0.5 * v_tli**2 - 1.0 / r_tli    # specific orbital energy
a_transfer = -1.0 / (2.0 * E_specific)         # semi-major axis
e_transfer = 1.0 - r_tli / a_transfer          # eccentricity (perigee at TLI)
r_apogee = 2 * a_transfer - r_tli              # apogee distance
p = a_transfer * (1 - e_transfer**2)           # semi-latus rectum

print(f"\n--- Trans-Lunar Injection (Apollo 11 data) ---")
print(f"  TLI altitude:      {alt_tli_km:.0f} km")
print(f"  TLI radius:        {r_tli:.4f} du ({R_earth_km + alt_tli_km:.0f} km)")
print(f"  TLI speed:         {v_tli:.5f} vu ({v_tli_ms:.0f} m/s = {v_tli_ms/1000:.3f} km/s)")
print(f"  Escape velocity:   {v_esc_tli:.5f} vu ({v_esc_tli*vu/1000:.3f} km/s)")
print(f"  v_TLI / v_escape:  {v_tli/v_esc_tli:.4f} (< 1 = bound orbit)")
print(f"  Specific energy:   {E_specific:.6f} (negative = bound)")
print(f"  Semi-major axis:   {a_transfer:.2f} du ({a_transfer*du:.0f} km)")
print(f"  Eccentricity:      {e_transfer:.5f}")
print(f"  Apogee distance:   {r_apogee:.2f} du ({r_apogee*du:.0f} km)")
print(f"  Note: apogee ({r_apogee:.0f} du) > Moon distance ({r_moon:.0f} du)")
print(f"  The Moon's gravity intervenes before apogee is reached.")

# ============================================================
# Transfer timing via Kepler's equation
# ============================================================
# Find the true anomaly at which the spacecraft crosses r_moon.
# The orbit equation: r = p / (1 + e*cos(nu))
# Solving for nu at r = r_moon:

r_target = r_moon            # target distance = Moon's orbit
cos_nu = (p / r_target - 1) / e_transfer
cos_nu = max(-1.0, min(1.0, cos_nu))
nu_target = math.acos(cos_nu)

# Convert true anomaly to time using Kepler's equation:
#   E = 2*atan(sqrt((1-e)/(1+e)) * tan(nu/2))   (eccentric anomaly)
#   M = E - e*sin(E)                              (mean anomaly)
#   t = M / (2*pi) * T                           (time from perigee)
tan_E_half = math.sqrt((1 - e_transfer) / (1 + e_transfer)) * math.tan(nu_target / 2)
E_target = 2.0 * math.atan(tan_E_half)
if E_target < 0:
    E_target += 2 * math.pi

M_target = E_target - e_transfer * math.sin(E_target)
T_transfer = 2 * math.pi * math.sqrt(a_transfer**3)
t_to_target = M_target / (2 * math.pi) * T_transfer

# How far the Moon moves during the spacecraft's transit
moon_angle_at_encounter = omega_moon * t_to_target  # radians CCW

print(f"\n--- Transfer timing (Kepler's equation) ---")
print(f"  True anomaly at r_moon: {math.degrees(nu_target):.1f} deg (apogee is at 180 deg)")
print(f"  Transfer orbit period: {T_transfer:.1f} tu ({T_transfer*tu/86400:.2f} days)")
print(f"  Transit time to r_moon: {t_to_target:.1f} tu ({t_to_target*tu/3600:.1f} hrs = {t_to_target*tu/86400:.2f} days)")
print(f"  Moon moves during transit: {math.degrees(moon_angle_at_encounter):.1f} deg CCW")

# ============================================================
# Launch geometry (targeting)
# ============================================================
# Key insight: the spacecraft crosses r_moon at true anomaly nu_target
# (172.2 deg from perigee), NOT at apogee (180 deg). The angular
# position of the spacecraft at that moment is:
#
#   spacecraft_angle = perigee_angle + nu_target
#
# For the spacecraft to arrive at the Moon:
#   perigee_angle + nu_target = moon_angle_at_encounter
#   perigee_angle = moon_angle_at_encounter - nu_target
#
# A small flyby offset ensures passage behind the Moon (far side)
# rather than collision. This offset creates the impact parameter
# that determines perilune altitude.

perigee_angle = moon_angle_at_encounter - nu_target

# Flyby offset: shifts aim point slightly behind Moon in its orbit.
# 0.7 deg produces ~200-400 km perilune altitude (Apollo 13 was 254 km).
# This is the most sensitive parameter in the entire calculation.
flyby_offset_deg = 0.7
perigee_angle += math.radians(flyby_offset_deg)

apogee_angle = perigee_angle + math.pi

print(f"\n--- Launch geometry ---")
print(f"  Moon position at encounter: {math.degrees(moon_angle_at_encounter):.1f} deg CCW from +x")
print(f"  Apogee direction:  {math.degrees(apogee_angle):.1f} deg CCW from +x")
print(f"  Perigee angle:     {math.degrees(perigee_angle):.1f} deg CCW from +x")
print(f"  Flyby offset:      {flyby_offset_deg} deg (for far-side passage)")

# ============================================================
# Compute spacecraft state vector at TLI
# ============================================================
# Position at perigee_angle, radius = r_tli:
#   CCW convention: pos = (r*cos(theta), 0, -r*sin(theta))
sc_x = r_tli * math.cos(perigee_angle)
sc_y = 0.0
sc_z = -r_tli * math.sin(perigee_angle)

# Velocity: tangential (CCW), speed = v_tli, flight path angle = 0:
#   CCW tangent at theta: v * (-sin(theta), 0, -cos(theta))
sc_vx = v_tli * (-math.sin(perigee_angle))
sc_vy = 0.0
sc_vz = v_tli * (-math.cos(perigee_angle))

print(f"\n--- Spacecraft initial state (TLI burnout) ---")
print(f"  Position: ({sc_x:.5f}, {sc_y:.5f}, {sc_z:.5f}) du")
print(f"            ({sc_x*du:.1f}, {sc_y*du:.1f}, {sc_z*du:.1f}) km")
print(f"  Velocity: ({sc_vx:.5f}, {sc_vy:.5f}, {sc_vz:.5f}) vu")
print(f"            ({sc_vx*vu/1000:.4f}, {sc_vy*vu/1000:.4f}, {sc_vz*vu/1000:.4f}) km/s")
print(f"  |v|:      {v_tli:.5f} vu ({v_tli*vu/1000:.4f} km/s)")

# ============================================================
# Encounter analysis (analytical estimate)
# ============================================================
# Spacecraft speed at Moon's distance (vis-viva equation, Earth-only):
v_at_moon = math.sqrt(max(0, 2.0/r_moon - 1.0/a_transfer))

# Decompose into radial and tangential at nu_target:
v_r_at_moon = math.sqrt(1.0/p) * e_transfer * math.sin(nu_target)
v_t_at_moon = math.sqrt(1.0/p) * (1 + e_transfer * math.cos(nu_target))

# Relative velocity between spacecraft and Moon at encounter:
# (determines hyperbolic flyby speed in Moon's frame)
delta_v_r = v_r_at_moon           # Moon has no radial velocity
delta_v_t = v_t_at_moon - v_moon  # tangential difference
v_relative = math.sqrt(delta_v_r**2 + delta_v_t**2)

# Impact parameter from the flyby offset angle
b_estimate = r_moon * math.sin(math.radians(flyby_offset_deg))

# Hyperbolic perilune from impact parameter and v_infinity:
#   r_p = (-c + sqrt(c^2 + 4*b^2)) / 2  where c = 2*mu_moon/v_inf^2
c_hyp = 2 * M_moon / v_relative**2
r_perilune = (-c_hyp + math.sqrt(c_hyp**2 + 4 * b_estimate**2)) / 2.0
alt_perilune = r_perilune - R_moon

# Hyperbolic deflection angle
e_hyp = 1 + r_perilune * v_relative**2 / M_moon
delta_angle = 2 * math.asin(1.0 / e_hyp)

print(f"\n--- Encounter analysis (two-body estimate) ---")
print(f"  Spacecraft speed at Moon distance: {v_at_moon:.5f} vu ({v_at_moon*vu/1000:.4f} km/s)")
print(f"    (Near apogee — very slow. This is why the flyby looks 'stalled'.)")
print(f"  Moon orbital speed: {v_moon:.5f} vu ({v_moon*vu/1000:.4f} km/s)")
print(f"  Relative velocity:  {v_relative:.5f} vu ({v_relative*vu/1000:.4f} km/s)")
print(f"    (This is the hyperbolic excess velocity in Moon's frame)")
print(f"  Impact parameter:   {b_estimate:.3f} du ({b_estimate*du:.0f} km)")
print(f"  Estimated perilune: {r_perilune:.4f} du ({r_perilune*du:.0f} km from center)")
print(f"  Estimated altitude: {alt_perilune:.4f} du ({alt_perilune*du:.0f} km above surface)")
print(f"  Deflection angle:   {math.degrees(delta_angle):.1f} deg")
print(f"  NOTE: This is a rough analytical estimate. The actual three-body")
print(f"  simulation determines the true perilune and deflection.")

# ============================================================
# Simulation parameters
# ============================================================
dt = 0.01           # 8.1 seconds — resolves lunar flyby dynamics
t_total_tu = 1200   # ~11 days — covers full outbound + return trip
n_steps = int(t_total_tu / dt)
cam_dist = r_moon * 1.5

print(f"\n--- Simulation parameters ---")
print(f"  dt:          {dt} tu ({dt*tu:.1f} s)")
print(f"  Total time:  {t_total_tu} tu ({t_total_tu*tu/86400:.1f} days)")
print(f"  Steps:       {n_steps:,}")
print(f"  r_soft:      0.01 du ({0.01*du:.0f} km)")
print(f"  Gravity:     P2P (exact, 3 bodies)")
print(f"  Camera:      (0, {cam_dist:.1f}, 0) — above orbital plane")

# ============================================================
# Script output
# ============================================================
print(f"\n{'='*65}")
print(f"SCRIPT FILE CONTENTS")
print(f"{'='*65}")
print(f"G               1.0")
print(f"FDE             0.0")
print(f"dt              {dt}")
print(f"r_soft          0.01")
print(f"Solver          LeapFrog")
print(f"Gravity         P2P")
print(f"N_SystemBodies  3")
print(f"Camera          0.0  {cam_dist:.1f}  0.0")
print(f"")
print(f"# Earth at origin (slight +z velocity for COM balance with Moon)")
print(f"Body  0   0.0 0.0 0.0   0.0 0.0 {v_earth:.6f}   1.0")
print(f"# Moon along +x, velocity in -z (CCW orbit viewed from +y)")
print(f"Body  0   {r_moon:.4f} 0.0 0.0   0.0 0.0 {-v_moon:.6f}   {M_moon:.6f}")
print(f"# Apollo CSM post-TLI (test particle)")
print(f"Body  0   {sc_x:.5f} {sc_y:.5f} {sc_z:.5f}   {sc_vx:.5f} {sc_vy:.5f} {sc_vz:.5f}   1.0e-20")

# ============================================================
# Tuning notes
# ============================================================
print(f"\n{'='*65}")
print(f"TUNING NOTES")
print(f"{'='*65}")
print(f"  The free-return trajectory is EXTREMELY sensitive to initial")
print(f"  conditions. A 0.1 deg change in flyby_offset_deg shifts the")
print(f"  perilune by hundreds of kilometers.")
print(f"")
print(f"  If spacecraft hits Moon:   increase flyby_offset_deg by 0.05")
print(f"  If spacecraft misses Moon: decrease flyby_offset_deg by 0.05")
print(f"  If no return to Earth:     try larger flyby_offset_deg changes")
print(f"")
print(f"  In the real Apollo program, mid-course corrections of 5-20 m/s")
print(f"  were applied during transit. This sensitivity is inherent to")
print(f"  the three-body problem.")
