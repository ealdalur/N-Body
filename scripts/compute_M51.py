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

=== Mass Decomposition Methodology ===

Following Salo & Laurikainen (2000): the disc contributes ~1/3 of the total
rotational support at the disc edge, and the halo contributes ~2/3. This gives
M_halo/M_disc ~ 2 within the disc radius. The baryonic (disc+bulge) masses
used here are rotation-curve-decomposed values, not photometric stellar masses.
This is standard for N-body tidal interaction studies — the dynamically cold
disc mass that participates in spiral arm formation is lower than the total
photometric mass (which includes dynamically hot thick-disc and bulge stars).

  V_total^2 = V_baryon^2 + V_halo^2
  V_baryon = V_total / sqrt(3)    [disc provides 1/3 of V^2]
  V_halo   = V_total * sqrt(2/3)  [halo provides 2/3 of V^2]

=== Unit system (G=1) ===
  1 distance unit = 60 pc = 0.060 kpc
  1 mass unit     = 10^4 solar masses
  1 velocity unit = 1 km/s
  1 time unit     = 60 pc / (1 km/s) = 58.7 Myr

=== Sources ===
  - Salo & Laurikainen 2000, MNRAS 319, 377 (orbit-constrained N-body model)
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
# Observed rotation curve: V_total = 210 km/s (flat)
# Rotation curve decomposition (Salo & Laurikainen methodology):
#   V_baryon = 210/sqrt(3) = 121 km/s  -> M_baryon = V_b^2 * R / G
#   V_halo = 210*sqrt(2/3) = 172 km/s  -> haloVc parameter
# This gives M_halo/M_baryon ~ 2 within the disc radius.
m51a_Vtotal_kms = 210               # observed flat rotation velocity
m51a_radius_kpc = 11.2              # stellar disc radius (R25)
m51a_inner_kpc = 0.3                # inner bulge region (no disc particles)
m51a_haloRc_kpc = 1.0               # DM halo core radius (Salo: ~0.5-1 kpc)

# Rotation curve decomposition
m51a_Vbaryon = m51a_Vtotal_kms / math.sqrt(3)
m51a_Vhalo = m51a_Vtotal_kms * math.sqrt(2.0/3.0)

# Baryonic mass from V_baryon^2 * R / G (G=1 in code units)
m51a_R_code = m51a_radius_kpc / du
m51a_baryon_code = m51a_Vbaryon**2 * m51a_R_code
m51a_baryon_msun = m51a_baryon_code * mu

# Bulge/disc decomposition: bulge ~ 17% of baryonic (typical for Sbc)
m51a_bulge_frac = 0.17
m51a_bulge_code = round(m51a_baryon_code * m51a_bulge_frac / 10000) * 10000
m51a_disc_code = m51a_baryon_code - m51a_bulge_code
m51a_Mfrac = m51a_disc_code / m51a_bulge_code

# Halo: cored isothermal V_halo(r) = haloVc * r / sqrt(r^2 + Rc^2)
# At disc edge, V_halo = haloVc * R / sqrt(R^2 + Rc^2) should equal m51a_Vhalo
m51a_haloRc_code = m51a_haloRc_kpc / du
m51a_haloVc = m51a_Vhalo / math.sqrt(m51a_R_code**2 / (m51a_R_code**2 + m51a_haloRc_code**2))

# === NGC 5195 (M51b) — Compact lenticular companion ===
# Observed/estimated rotation curve: V_total = 130 km/s
# Same decomposition: M_halo/M_baryon ~ 2 within disc radius.
m51b_Vtotal_kms = 130               # estimated from mass-Vc scaling
m51b_radius_kpc = 2.1               # optical extent
m51b_inner_kpc = 0.2                # inner region
m51b_haloRc_kpc = 0.6               # DM halo core radius (compact galaxy)

m51b_Vbaryon = m51b_Vtotal_kms / math.sqrt(3)
m51b_Vhalo = m51b_Vtotal_kms * math.sqrt(2.0/3.0)

m51b_R_code = m51b_radius_kpc / du
m51b_baryon_code = m51b_Vbaryon**2 * m51b_R_code
m51b_baryon_msun = m51b_baryon_code * mu

# Bulge/disc decomposition: bulge ~ 40% for SB0 lenticular
m51b_bulge_frac = 0.40
m51b_bulge_code = round(m51b_baryon_code * m51b_bulge_frac / 10000) * 10000
m51b_disc_code = m51b_baryon_code - m51b_bulge_code
m51b_Mfrac = m51b_disc_code / m51b_bulge_code

m51b_haloRc_code = m51b_haloRc_kpc / du
m51b_haloVc = m51b_Vhalo / math.sqrt(m51b_R_code**2 / (m51b_R_code**2 + m51b_haloRc_code**2))

# === Interaction geometry ===
# Source: Salo & Laurikainen 2000; Dobbs et al. 2010
pericenter_kpc = 15.0               # closest approach distance
orbital_inclination_deg = 15.0      # M51b orbit tilted slightly from disc plane

# ============================================================
print("=" * 65)
print("M51 WHIRLPOOL GALAXY INTERACTION PARAMETERS")
print("=" * 65)

# === Print M51a parameters ===
print("\n--- NGC 5194 (M51a) ---")
print(f"  Observed V_total:  {m51a_Vtotal_kms} km/s")
print(f"  V_baryon:          {m51a_Vbaryon:.1f} km/s (1/3 of V^2)")
print(f"  V_halo:            {m51a_Vhalo:.1f} km/s (2/3 of V^2)")
print(f"  Disc radius:       {m51a_R_code:.1f} code units ({m51a_radius_kpc} kpc)")
print(f"  Inner radius:      {m51a_inner_kpc/du:.1f} code units ({m51a_inner_kpc} kpc)")
print(f"  Bulge mass (M):    {m51a_bulge_code:.0f} code units ({m51a_bulge_code*mu:.2e} Msun)")
print(f"  Disc mass:         {m51a_disc_code:.0f} code units ({m51a_disc_code*mu:.2e} Msun)")
print(f"  Total baryonic:    {m51a_baryon_code:.0f} code units ({m51a_baryon_msun:.2e} Msun)")
print(f"  Mfrac:             {m51a_Mfrac:.2f}")
print(f"  Halo Vc:           {m51a_haloVc:.1f} km/s")
print(f"  Halo Rc:           {m51a_haloRc_code:.1f} code units ({m51a_haloRc_kpc} kpc)")

# Verify rotation curve at disc edge
v_b_check = math.sqrt(m51a_baryon_code / m51a_R_code)
v_h_check = m51a_haloVc * m51a_R_code / math.sqrt(m51a_R_code**2 + m51a_haloRc_code**2)
v_total_check = math.sqrt(v_b_check**2 + v_h_check**2)
print(f"  V_total check at R: {v_total_check:.1f} km/s (target: {m51a_Vtotal_kms})")

# M_halo/M_baryon within disc
m51a_Mhalo_disc = m51a_haloVc**2 * m51a_R_code**3 / (m51a_R_code**2 + m51a_haloRc_code**2)
print(f"  M_halo within disc: {m51a_Mhalo_disc:.0f} code units")
print(f"  M_halo/M_baryon:   {m51a_Mhalo_disc/m51a_baryon_code:.2f} (target: ~2.0)")

# === Print M51b parameters ===
print(f"\n--- NGC 5195 (M51b) ---")
print(f"  Observed V_total:  {m51b_Vtotal_kms} km/s")
print(f"  V_baryon:          {m51b_Vbaryon:.1f} km/s")
print(f"  V_halo:            {m51b_Vhalo:.1f} km/s")
print(f"  Disc radius:       {m51b_R_code:.1f} code units ({m51b_radius_kpc} kpc)")
print(f"  Inner radius:      {m51b_inner_kpc/du:.1f} code units ({m51b_inner_kpc} kpc)")
print(f"  Bulge mass (M):    {m51b_bulge_code:.0f} code units ({m51b_bulge_code*mu:.2e} Msun)")
print(f"  Disc mass:         {m51b_disc_code:.0f} code units ({m51b_disc_code*mu:.2e} Msun)")
print(f"  Total baryonic:    {m51b_baryon_code:.0f} code units ({m51b_baryon_msun:.2e} Msun)")
print(f"  Mfrac:             {m51b_Mfrac:.2f}")
print(f"  Halo Vc:           {m51b_haloVc:.1f} km/s")
print(f"  Halo Rc:           {m51b_haloRc_code:.1f} code units ({m51b_haloRc_kpc} kpc)")

v_b_check = math.sqrt(m51b_baryon_code / m51b_R_code)
v_h_check = m51b_haloVc * m51b_R_code / math.sqrt(m51b_R_code**2 + m51b_haloRc_code**2)
v_total_check = math.sqrt(v_b_check**2 + v_h_check**2)
print(f"  V_total check at R: {v_total_check:.1f} km/s (target: {m51b_Vtotal_kms})")

m51b_Mhalo_disc = m51b_haloVc**2 * m51b_R_code**3 / (m51b_R_code**2 + m51b_haloRc_code**2)
print(f"  M_halo within disc: {m51b_Mhalo_disc:.0f} code units")
print(f"  M_halo/M_baryon:   {m51b_Mhalo_disc/m51b_baryon_code:.2f} (target: ~2.0)")

# === Mass ratio at pericenter ===
r_peri_code = pericenter_kpc / du
m51a_halo_peri = m51a_haloVc**2 * r_peri_code**3 / (r_peri_code**2 + m51a_haloRc_code**2)
m51b_halo_peri = m51b_haloVc**2 * r_peri_code**3 / (r_peri_code**2 + m51b_haloRc_code**2)
m51a_dyn_peri = m51a_baryon_code + m51a_halo_peri
m51b_dyn_peri = m51b_baryon_code + m51b_halo_peri
dyn_ratio = m51b_dyn_peri / m51a_dyn_peri

print(f"\n--- Mass ratio (dynamical, at pericenter = {pericenter_kpc} kpc) ---")
print(f"  M51a enclosed: {m51a_dyn_peri:.0f} (baryonic {m51a_baryon_code:.0f} + halo {m51a_halo_peri:.0f})")
print(f"  M51b enclosed: {m51b_dyn_peri:.0f} (baryonic {m51b_baryon_code:.0f} + halo {m51b_halo_peri:.0f})")
print(f"  M51b / M51a (dynamical): 1:{1/dyn_ratio:.1f}")
print(f"  Observational constraint: 1:3 to 1:5 (Salo & Laurikainen 2000)")

# === Orbital setup ===
# Place M51b at apocenter (~40 kpc) with purely tangential velocity.
# Energy and angular momentum conservation gives the velocity needed
# to reach the desired pericenter.
#
# The simulation applies halo acceleration toward the halo center:
#   a_halo = Vc^2 * r / (r^2 + Rc^2)  (attractive, directed inward)
#
# The corresponding gravitational potential (with a_r = -dPhi/dr):
#   Phi_halo(r) = +0.5 * Vc^2 * ln(r^2 + Rc^2)
#
# This INCREASES outward (particle at larger r has higher potential energy
# and falls inward). For baryonic point masses: Phi_baryon(r) = -G*M/r
# (also increases outward from -inf to 0).
#
# The mutual force between the two systems includes BOTH halos (by Newton's
# 3rd law — M51a's particles feel M51b's halo and vice versa), plus baryons:
#   Phi(r) = +0.5*Vc_a^2*ln(r^2+Rc_a^2) + 0.5*Vc_b^2*ln(r^2+Rc_b^2) - M_baryon/r

print("\n--- Orbital computation ---")

r_peri = pericenter_kpc / du
r_start_kpc = 40.0
r_start = r_start_kpc / du

print(f"  Pericenter:        {r_peri:.1f} code units ({pericenter_kpc} kpc)")
print(f"  Initial separation: {r_start:.1f} code units ({r_start_kpc} kpc)")

M_baryon_total = m51a_baryon_code + m51b_baryon_code

def Phi(r):
    phi_halo_a = 0.5 * m51a_haloVc**2 * math.log(r**2 + m51a_haloRc_code**2)
    phi_halo_b = 0.5 * m51b_haloVc**2 * math.log(r**2 + m51b_haloRc_code**2)
    phi_baryon = -M_baryon_total / r
    return phi_halo_a + phi_halo_b + phi_baryon

Phi_start = Phi(r_start)
Phi_peri = Phi(r_peri)

print(f"  Phi(r_start):      {Phi_start:.1f} (higher = further from potential minimum)")
print(f"  Phi(r_peri):       {Phi_peri:.1f}")
print(f"  Delta Phi:         {Phi_peri - Phi_start:.1f} (negative = peri is lower)")

# Energy and angular momentum conservation:
#   At apocenter (r_start), v_r = 0:  E = 0.5*v_t^2 + Phi(r_start)
#   At pericenter (r_peri), v_r = 0:  E = 0.5*v_p^2 + Phi(r_peri)
#   Angular momentum: L = r_start * v_t = r_peri * v_p
#
# Solving (same algebra, now dPhi < 0 and denominator < 0, so v_t^2 > 0):
#   v_t^2 = 2*(Phi_peri - Phi_start) / (1 - (r_start/r_peri)^2)

numerator = 2.0 * (Phi_peri - Phi_start)
denominator = 1.0 - (r_start / r_peri)**2
v_t_sq = numerator / denominator
v_t = math.sqrt(v_t_sq)

v_p = r_start * v_t / r_peri
E_specific = 0.5 * v_t**2 + Phi_start

print(f"\n  Apocenter approach (v_r = 0 at r_start):")
print(f"  v_tangential at r_start: {v_t:.1f} km/s")
print(f"  v_pericenter:            {v_p:.1f} km/s")
print(f"  Specific energy:         {E_specific:.1f}")
print(f"  Angular momentum L:      {r_start * v_t:.0f}")

# Circular velocity at r_start for reference
v_circ_start = math.sqrt(M_baryon_total/r_start
               + m51a_haloVc**2 * r_start**2 / (r_start**2 + m51a_haloRc_code**2)
               + m51b_haloVc**2 * r_start**2 / (r_start**2 + m51b_haloRc_code**2))
print(f"  v_circular at r_start:   {v_circ_start:.1f} km/s")
print(f"  v_t / v_circ:            {v_t/v_circ_start:.3f}")

# Time estimate
r_mean = (r_start + r_peri) / 2.0
v_circ_mean = math.sqrt(M_baryon_total/r_mean
              + m51a_haloVc**2 * r_mean**2 / (r_mean**2 + m51a_haloRc_code**2)
              + m51b_haloVc**2 * r_mean**2 / (r_mean**2 + m51b_haloRc_code**2))
T_orbit_approx = 2 * math.pi * r_mean / v_circ_mean
t_to_peri = T_orbit_approx / 2.0

print(f"\n  Approximate half-orbit time to pericenter:")
print(f"    {t_to_peri:.0f} code units = {t_to_peri*tu:.0f} Myr")
print(f"    At dt=0.0005: {int(t_to_peri/0.0005)} steps")

# === Coordinate system ===
# M51a at origin, disc in x-z plane (normal +y).
# M51a disc rotates CLOCKWISE viewed from +y (because LoadGalaxyDiscState
# uses v_tan = -vm, and t_hat = cross(r_hat, n) = +z at +x position,
# so particles at +x move in -z direction = CW viewed from +y).
#
# For a PROGRADE encounter, M51b must orbit in the same sense (CW from +y),
# so its tangential velocity at +x must be in the -z direction.
# Orbital inclination tilts M51b velocity slightly out of disc plane.

inc = math.radians(orbital_inclination_deg)

pos_x = r_start * math.cos(inc)
pos_y = 0.0
pos_z = 0.0

vel_x = 0.0
vel_y = -v_t * math.sin(inc)
vel_z = -v_t * math.cos(inc)

print(f"\n--- Initial conditions (code units) ---")
print(f"  M51b position: ({pos_x:.1f}, {pos_y:.1f}, {pos_z:.1f})")
print(f"  M51b velocity: ({vel_x:.1f}, {vel_y:.1f}, {vel_z:.1f})")
print(f"  |velocity|:    {math.sqrt(vel_x**2+vel_y**2+vel_z**2):.1f} km/s")

# === Camera ===
cam_height = m51a_R_code * 4.0
print(f"\n--- Camera ---")
print(f"  Position: (0, {cam_height:.0f}, 0)")

# === Particle counts ===
n_m51a = 40000
n_m51b = 10000
n_total = n_m51a + n_m51b

print(f"\n--- Particle counts ---")
print(f"  M51a: {n_m51a} ({n_m51a-1} disc + 1 central)")
print(f"  M51b: {n_m51b} ({n_m51b-1} disc + 1 central)")
print(f"  Total: {n_total}")
print(f"  M51a particle mass: ~{m51a_disc_code*mu/(n_m51a-1):.0f} Msun")
print(f"  M51b particle mass: ~{m51b_disc_code*mu/(n_m51b-1):.0f} Msun")

# === M51b disc orientation ===
m51b_normal_x = math.sin(inc)
m51b_normal_y = math.cos(inc)
m51b_normal_z = 0.0

print(f"\n--- M51b disc normal ---")
print(f"  ({m51b_normal_x:.4f}, {m51b_normal_y:.4f}, {m51b_normal_z:.4f})")
print(f"  (tilted {orbital_inclination_deg} deg from M51a disc normal)")

# === Script output ===
m51a_Ri_code = m51a_inner_kpc / du
m51b_Ri_code = m51b_inner_kpc / du

print(f"\n{'='*65}")
print(f"SCRIPT LINES")
print(f"{'='*65}")
print(f"N_SystemBodies  {n_m51a}  {n_m51b}")
print(f"Camera          0.0  {cam_height:.0f}  0.0")
print(f"")
print(f"# M51a (NGC 5194) at origin, disc in x-z plane")
print(f"GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   {m51a_bulge_code:.1f} {m51a_Mfrac:.2f} {m51a_R_code:.1f} {m51a_Ri_code:.1f} 1.5  {m51a_haloVc:.1f} {m51a_haloRc_code:.1f}")
print(f"")
print(f"# M51b (NGC 5195) approaching on prograde orbit")
print(f"GalaxyDisc  1   {pos_x:.1f} {pos_y:.1f} {pos_z:.1f}   {vel_x:.1f} {vel_y:.1f} {vel_z:.1f}   {m51b_normal_x:.4f} {m51b_normal_y:.4f} {m51b_normal_z:.4f}   {m51b_bulge_code:.1f} {m51b_Mfrac:.2f} {m51b_R_code:.1f} {m51b_Ri_code:.1f} 1.5  {m51b_haloVc:.1f} {m51b_haloRc_code:.1f}")

# === Timeline ===
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
