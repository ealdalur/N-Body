"""
Compute initial conditions for an M51 (Whirlpool Galaxy) interaction simulation.

The M51 system consists of NGC 5194 (M51a, a grand-design Sbc spiral) and
NGC 5195 (M51b, a compact SB0-pec lenticular companion). M51b has already
made at least one close prograde passage through M51a's disc, exciting the
prominent two-arm spiral structure that makes this system famous.

This script computes pre-encounter initial conditions that will produce
M51-like tidal spiral arms when the simulation is run forward.

=== Key Observational Constraints ===

From Salo & Laurikainen 2000, Dobbs et al. 2010, and others:
  - Encounter is PROGRADE (M51b orbits in same sense as M51a disc rotation)
  - Pericenter distance: ~15-25 kpc (250-420 code units)
  - Time since first passage: ~300-400 Myr
  - M51b is currently behind M51a's disc plane (farther from Earth)
  - Orbital inclination to M51a disc: ~10-20 degrees
  - Mass ratio (total dynamical) M51b:M51a ~ 1:3 to 1:4
  - The prograde nature is critical for strong spiral arm excitation

=== Unit system (G=1) ===
  1 distance unit = 60 pc = 0.060 kpc
  1 mass unit     = 10^4 solar masses
  1 velocity unit = 1 km/s
  1 time unit     = 60 pc / (1 km/s) = 58.7 Myr

=== Sources ===
  - Salo & Laurikainen 2000, A&A 359, 895 (orbit-constrained N-body model)
  - Dobbs et al. 2010, MNRAS 403, 625 (SPH + N-body with gas)
  - Querejeta et al. 2015 (mass decomposition from Spitzer 3.6um)
  - Mentuch Cooper et al. 2012 (stellar masses)
  - McQuinn et al. 2016 (TRGB distance 8.58 Mpc)
  - Sofue et al. 1999; Meidt et al. 2013 (rotation curve)
  - Toomre & Toomre 1972 (original tidal interaction model)
"""

import math

# === Unit system ===
du = 0.060   # 1 code distance unit = 60 pc = 0.060 kpc
mu = 1.0e4   # 1 code mass unit = 10^4 Msun
vu = 1.0     # 1 code velocity unit = 1 km/s
tu = 58.7    # 1 code time unit = 60 pc / (1 km/s) = 58.7 Myr

# === NGC 5194 (M51a) — Grand-design spiral ===
# Source: Querejeta et al. 2015; Mentuch Cooper et al. 2012
m51a_bulge_msun = 1.0e10        # bulge mass (central concentration)
m51a_disc_msun = 4.0e10         # stellar disc mass
m51a_gas_msun = 0.9e10          # HI + H2 gas mass
m51a_radius_kpc = 11.2          # stellar disc radius (R25)
m51a_inner_kpc = 0.3            # inner bulge region (no disc particles)
m51a_vc_kms = 210               # flat rotation curve velocity
m51a_haloRc_kpc = 5.0           # DM halo core radius

# === NGC 5195 (M51b) — Compact lenticular companion ===
# Source: Mentuch Cooper et al. 2012; Salo & Laurikainen 2000
m51b_mass_msun = 1.0e10         # total baryonic mass (bulge-dominated)
m51b_disc_msun = 0.5e10         # small "disc" component (disturbed)
m51b_radius_kpc = 2.1           # optical extent
m51b_inner_kpc = 0.2            # inner region
m51b_vc_kms = 130               # estimated from mass-Vc scaling
m51b_haloRc_kpc = 2.5           # DM halo core radius

# === Interaction geometry ===
# Source: Salo & Laurikainen 2000; Dobbs et al. 2010
pericenter_kpc = 20.0           # closest approach distance
orbital_inclination_deg = 15.0  # M51b orbit tilted slightly from disc plane
distance_mpc = 8.0              # adopted distance (McQuinn et al. 2016: 8.58)

# ============================================================
print("=" * 65)
print("M51 WHIRLPOOL GALAXY INTERACTION PARAMETERS")
print("=" * 65)

# === Convert to code units ===
print("\n--- NGC 5194 (M51a) ---")
m51a_R = m51a_radius_kpc / du
m51a_Ri = m51a_inner_kpc / du
m51a_M = m51a_bulge_msun / mu
m51a_Mfrac = (m51a_disc_msun + m51a_gas_msun) / m51a_bulge_msun
m51a_haloVc = m51a_vc_kms / vu
m51a_haloRc = m51a_haloRc_kpc / du
m51a_baryonic = m51a_M * (1 + m51a_Mfrac)

print(f"  Disc radius:       {m51a_R:.1f} code units ({m51a_radius_kpc} kpc)")
print(f"  Inner radius:      {m51a_Ri:.1f} code units ({m51a_inner_kpc} kpc)")
print(f"  Bulge mass (M):    {m51a_M:.0f} code units ({m51a_bulge_msun:.1e} Msun)")
print(f"  Disc+gas mass:     {(m51a_disc_msun+m51a_gas_msun)/mu:.0f} code units ({m51a_disc_msun+m51a_gas_msun:.1e} Msun)")
print(f"  Mfrac:             {m51a_Mfrac:.2f}")
print(f"  Halo Vc:           {m51a_haloVc:.1f} ({m51a_vc_kms} km/s)")
print(f"  Halo Rc:           {m51a_haloRc:.1f} code units ({m51a_haloRc_kpc} kpc)")
print(f"  Total baryonic:    {m51a_baryonic:.0f} code units ({m51a_baryonic*mu:.1e} Msun)")

print("\n--- NGC 5195 (M51b) ---")
m51b_R = m51b_radius_kpc / du
m51b_Ri = m51b_inner_kpc / du
m51b_M = m51b_mass_msun / mu
m51b_Mfrac = m51b_disc_msun / m51b_mass_msun
m51b_haloVc = m51b_vc_kms / vu
m51b_haloRc = m51b_haloRc_kpc / du
m51b_baryonic = m51b_M * (1 + m51b_Mfrac)

print(f"  Radius:            {m51b_R:.1f} code units ({m51b_radius_kpc} kpc)")
print(f"  Inner radius:      {m51b_Ri:.1f} code units ({m51b_inner_kpc} kpc)")
print(f"  Central mass (M):  {m51b_M:.0f} code units ({m51b_mass_msun:.1e} Msun)")
print(f"  Disc mass:         {m51b_disc_msun/mu:.0f} code units ({m51b_disc_msun:.1e} Msun)")
print(f"  Mfrac:             {m51b_Mfrac:.2f}")
print(f"  Halo Vc:           {m51b_haloVc:.1f} ({m51b_vc_kms} km/s)")
print(f"  Halo Rc:           {m51b_haloRc:.1f} code units ({m51b_haloRc_kpc} kpc)")
print(f"  Total baryonic:    {m51b_baryonic:.0f} code units ({m51b_baryonic*mu:.1e} Msun)")

# === Mass ratio ===
mass_ratio = m51b_baryonic / m51a_baryonic
print(f"\n--- Mass ratio ---")
print(f"  M51b / M51a (baryonic): 1:{1/mass_ratio:.1f}")

# === Orbital setup ===
# We want M51b to approach M51a on a prograde, slightly inclined orbit
# and reach pericenter at ~20 kpc (333 code units).
#
# Strategy: place M51b at ~40 kpc (667 code units) from M51a on an
# infalling trajectory. Compute the velocity from energy conservation
# for a bound (elliptical) orbit with the desired pericenter.
#
# For two-body orbit:
#   E = (1/2)*v^2 - M_eff/r  (specific orbital energy, per unit reduced mass)
#   L = r * v_perp  (specific angular momentum)
#   At pericenter: L = r_p * v_p, E = v_p^2/2 - M_eff/r_p
#
# Total effective gravitating mass at separation r (halo-dominated):
# M_eff(r) ~ M51a_halo(r) + M51b_halo(r) + baryonic
# Using isothermal halo: M_halo(r) = Vc^2 * r [G=1]

print("\n--- Orbital computation ---")

r_peri = pericenter_kpc / du  # pericenter in code units
r_start_kpc = 40.0           # initial separation
r_start = r_start_kpc / du   # in code units

print(f"  Pericenter:        {r_peri:.1f} code units ({pericenter_kpc} kpc)")
print(f"  Initial separation: {r_start:.1f} code units ({r_start_kpc} kpc)")

# Effective mass at initial separation (halo + baryonic)
# Each halo contributes Vc^2 * r at distance r from its center
# At separation r_start, M51a's halo encloses Vc_a^2 * r_start
# M51b's halo contribution is smaller (its mass is inside r_start)
M_eff_start = (m51a_haloVc**2 * r_start + m51b_haloVc**2 * r_start * 0.5
               + m51a_baryonic + m51b_baryonic)

# Effective mass at pericenter
M_eff_peri = (m51a_haloVc**2 * r_peri + m51b_haloVc**2 * r_peri * 0.5
              + m51a_baryonic + m51b_baryonic)

print(f"  M_eff at r_start:  {M_eff_start:.2e} code units")
print(f"  M_eff at r_peri:   {M_eff_peri:.2e} code units")

# For a bound orbit, we need to specify either eccentricity or apocenter.
# Use energy conservation + angular momentum conservation:
#   E = v^2/2 - M_eff/r = v_p^2/2 - M_eff/r_p
#   L = r * v_perp = r_p * v_p (at pericenter, all velocity is tangential)
#
# At initial position r_start, the velocity has radial and tangential components:
#   v_r = radial (inward) component
#   v_t = tangential component
#   L = r_start * v_t = r_p * v_p
#   v_t = r_p * v_p / r_start
#
# From energy conservation (using average M_eff for simplicity):
# Use the M_eff at the geometric mean distance for a rough orbit average:
M_eff_avg = (M_eff_start + M_eff_peri) / 2.0

# More carefully: for a Keplerian orbit with M_eff constant:
# E = -M_eff / (r_p + r_a)  where r_a is apocenter
# L = sqrt(M_eff * r_p * r_a * 2 / (r_p + r_a)) -- wrong, let me redo
#
# Actually, for a Keplerian orbit:
# r_p = a(1-e), r_a = a(1+e)
# E = -M/(2a), L^2 = M*a*(1-e^2) = M * r_p * (1+e) = M * r_p * (r_p + r_a) / a...
#
# Simpler approach: treat as two-point energy/angular momentum conservation
# with M_eff evaluated at each radius.
#
# At r_start: E = (1/2)*(v_r^2 + v_t^2) - M_eff_start / r_start
# At r_peri:  E = (1/2)*v_p^2 - M_eff_peri / r_peri
# Angular momentum: r_start * v_t = r_peri * v_p
#
# From L: v_p = r_start * v_t / r_peri
# Substituting into energy equation:
# (1/2)*(v_r^2 + v_t^2) - M_eff_start/r_start = (1/2)*(r_start*v_t/r_peri)^2 - M_eff_peri/r_peri
#
# For a "just arriving" scenario, let's set v_r such that the orbit has
# apocenter at r_start (i.e., r_start IS the apocenter, v_r = 0 there).
# Then all velocity at r_start is tangential.
#
# With v_r = 0 at r_start (apocenter):
# L = r_start * v_t
# v_p = r_start * v_t / r_peri
#
# Energy conservation:
# (1/2)*v_t^2 - M_eff_start/r_start = (1/2)*(r_start/r_peri)^2*v_t^2 - M_eff_peri/r_peri
# v_t^2 * [1 - (r_start/r_peri)^2] = 2*(M_eff_start/r_start - M_eff_peri/r_peri)
# v_t^2 = 2*(M_eff_start/r_start - M_eff_peri/r_peri) / [1 - (r_start/r_peri)^2]

# But (r_start/r_peri)^2 > 1 since r_start > r_peri, so denominator is negative.
# And M_eff_start/r_start < M_eff_peri/r_peri (since M_eff grows slower than r),
# wait no: M_eff = Vc^2*r + const, so M_eff/r ~ Vc^2 + const/r, which DECREASES with r.
# So M_eff_peri/r_peri > M_eff_start/r_start. Numerator is negative. Negative/negative = positive. Good.

numerator = 2.0 * (M_eff_start / r_start - M_eff_peri / r_peri)
denominator = 1.0 - (r_start / r_peri)**2
v_t_sq = numerator / denominator
v_t = math.sqrt(v_t_sq)

v_p = r_start * v_t / r_peri
E_specific = 0.5 * v_t**2 - M_eff_start / r_start

print(f"\n  Apocenter approach (v_r = 0 at r_start):")
print(f"  v_tangential at r_start: {v_t:.1f} km/s")
print(f"  v_pericenter:            {v_p:.1f} km/s")
print(f"  Specific energy:         {E_specific:.1f} (< 0 = bound)")
print(f"  Angular momentum L:      {r_start * v_t:.0f}")

# Time estimate: for an elliptical orbit in a logarithmic potential
# (isothermal halo), the orbital period is roughly:
# T ~ 2*pi*r / v_circular ~ 2*pi * r_mean / Vc
r_mean = (r_start + r_peri) / 2.0
T_orbit_approx = 2 * math.pi * r_mean / m51a_haloVc
t_to_peri = T_orbit_approx / 2.0  # half orbit from apocenter to pericenter

print(f"\n  Approximate half-orbit time to pericenter:")
print(f"    {t_to_peri:.0f} code units = {t_to_peri*tu:.0f} Myr")
print(f"    At dt=0.0005: {int(t_to_peri/0.0005)} steps")

# === Set up coordinate system ===
# M51a at origin with disc in x-z plane (normal along +y).
# M51b approaches in a prograde orbit (same direction as disc rotation).
# M51a rotates counter-clockwise viewed from +y.
# Prograde means M51b orbits counter-clockwise viewed from +y.
#
# Place M51b at apocenter along +x axis.
# Tangential velocity is in the +z direction (CCW when viewed from +y).
#
# Slight orbital inclination: tilt M51b's orbit 15 deg from the disc plane.
# This means M51b's velocity has a small +y component.

inc = math.radians(orbital_inclination_deg)

# Position: along +x at r_start, slightly above disc plane
pos_x = r_start * math.cos(inc)  # ~= r_start (inc is small)
pos_y = 0.0                       # start in-plane, orbit will carry it above/below
pos_z = 0.0

# Velocity: tangential (prograde = +z direction when at +x position)
# with small y-component for orbital inclination
vel_x = 0.0  # at apocenter, no radial velocity
vel_y = v_t * math.sin(inc)   # small out-of-plane component
vel_z = v_t * math.cos(inc)   # main tangential component

print(f"\n--- Initial conditions (code units) ---")
print(f"  M51b position: ({pos_x:.1f}, {pos_y:.1f}, {pos_z:.1f})")
print(f"  M51b velocity: ({vel_x:.1f}, {vel_y:.1f}, {vel_z:.1f})")
print(f"  |velocity|:    {math.sqrt(vel_x**2+vel_y**2+vel_z**2):.1f} km/s")

# === Camera placement ===
# M51a is nearly face-on (20 deg inclination to line of sight).
# With disc in x-z plane and normal along +y, the observer sees it from +y.
# Camera above, slightly offset to see depth.
cam_height = m51a_R * 4.0  # enough to see the full interaction
print(f"\n--- Camera ---")
print(f"  Position: (0, {cam_height:.0f}, 0)")

# === Particle counts ===
# M51a is the dominant galaxy — give it more particles for spiral resolution.
# M51b is compact — fewer particles suffice.
n_m51a = 40000
n_m51b = 10000
n_total = n_m51a + n_m51b
m51a_particle_msun = (m51a_disc_msun + m51a_gas_msun) / (n_m51a - 1)
m51b_particle_msun = m51b_disc_msun / (n_m51b - 1)

print(f"\n--- Particle counts ---")
print(f"  M51a: {n_m51a} ({n_m51a-1} disc + 1 central)")
print(f"  M51b: {n_m51b} ({n_m51b-1} disc + 1 central)")
print(f"  Total: {n_total}")
print(f"  M51a particle mass: ~{m51a_particle_msun:.0f} Msun")
print(f"  M51b particle mass: ~{m51b_particle_msun:.0f} Msun")

# === Disc orientation for M51b ===
# M51b is a disturbed lenticular. Its disc (what remains of it) is likely
# roughly aligned with its orbital plane. Use a normal tilted 15 deg from +y.
# Normal = (sin(inc), cos(inc), 0) = (0.259, 0.966, 0.0)
m51b_normal_x = math.sin(inc)
m51b_normal_y = math.cos(inc)
m51b_normal_z = 0.0

print(f"\n--- M51b disc normal ---")
print(f"  ({m51b_normal_x:.4f}, {m51b_normal_y:.4f}, {m51b_normal_z:.4f})")
print(f"  (tilted {orbital_inclination_deg} deg from M51a disc normal)")

# === Summary ===
print(f"\n{'='*65}")
print(f"SCRIPT LINES")
print(f"{'='*65}")
print(f"N_SystemBodies  {n_m51a}  {n_m51b}")
print(f"Camera          0.0  {cam_height:.0f}  0.0")
print(f"")
print(f"# M51a (NGC 5194) at origin, disc in x-z plane")
print(f"GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   {m51a_M:.1f} {m51a_Mfrac:.2f} {m51a_R:.1f} {m51a_Ri:.1f} 0.1  {m51a_haloVc:.1f} {m51a_haloRc:.1f}")
print(f"")
print(f"# M51b (NGC 5195) approaching on prograde orbit")
print(f"GalaxyDisc  1   {pos_x:.1f} {pos_y:.1f} {pos_z:.1f}   {vel_x:.1f} {vel_y:.1f} {vel_z:.1f}   {m51b_normal_x:.4f} {m51b_normal_y:.4f} {m51b_normal_z:.4f}   {m51b_M:.1f} {m51b_Mfrac:.2f} {m51b_R:.1f} {m51b_Ri:.1f} 0.1  {m51b_haloVc:.1f} {m51b_haloRc:.1f}")

# === Timing expectations ===
print(f"\n{'='*65}")
print(f"EXPECTED TIMELINE")
print(f"{'='*65}")
print(f"  t = 0:           M51b at apocenter ({r_start_kpc:.0f} kpc from M51a)")
print(f"  t ~ {t_to_peri/2:.0f}:       Spiral arms begin developing (M51b approaching)")
print(f"  t ~ {t_to_peri:.0f}:       First pericenter passage ({pericenter_kpc:.0f} kpc)")
print(f"  t ~ {t_to_peri*1.5:.0f}:      Strong tidal tails, M51b receding")
print(f"  t ~ {T_orbit_approx:.0f}:      M51b returns for second passage")
print(f"")
print(f"  In physical time:")
print(f"  t_pericenter ~ {t_to_peri*tu:.0f} Myr")
print(f"  t_full_orbit ~ {T_orbit_approx*tu:.0f} Myr")
print(f"  Comparable to Salo & Laurikainen 2000 model (~400 Myr between passages)")
