"""
Compute physical parameters for Milky Way - Andromeda simulation script.

Converts real astrophysical measurements into the N-Body code unit system:
  1 distance unit = 60 pc = 0.060 kpc
  1 mass unit     = 10^4 solar masses
  1 velocity unit = 1 km/s
  1 time unit     = 58.7 Myr  (= 60 pc / 1 km/s)

Sources:
  - MW rotation curve: Eilers et al. 2019 (~220 km/s)
  - MW baryonic mass: Bland-Hawthorn & Gerhard 2016 (~6e10 Msun)
  - Sgr A* mass: GRAVITY Collaboration 2019 (~4e6 Msun)
  - M31 mass: Tamm et al. 2012 (~1e11 Msun baryonic)
  - M31* mass: Bender et al. 2005 (~1.4e8 Msun)
  - MW-M31 distance: Riess et al. 2012 (~780 kpc)
  - M31 approach velocity: van der Marel et al. 2012 (~110 km/s radial)
  - M31 inclination: de Vaucouleurs 1958 (~77 deg)
"""

import math

# === Unit system ===
du = 0.060   # 1 code distance unit = 60 pc = 0.060 kpc
mu = 1.0e4   # 1 code mass unit = 10^4 solar masses
vu = 1.0     # 1 code velocity unit = 1 km/s

# === Milky Way parameters ===
mw_radius_kpc = 26.8        # stellar disc radius
mw_inner_kpc = 1.0          # inner hole (bulge region)
mw_bulge_msun = 1.5e10      # bulge + central mass
mw_disc_msun = 4.5e10       # disc mass (stars + gas)
mw_smbh_msun = 4.0e6        # Sgr A*
mw_vc_kms = 220             # flat rotation curve velocity
mw_haloRc_kpc = 10.0        # dark matter halo core radius

# === Andromeda parameters ===
and_radius_kpc = 33.5        # stellar disc radius
and_inner_kpc = 1.5          # inner hole
and_bulge_msun = 3.0e10      # bulge + central mass
and_disc_msun = 7.0e10       # disc mass
and_smbh_msun = 1.4e8        # M31*
and_vc_kms = 225             # flat rotation curve velocity
and_haloRc_kpc = 12.0        # dark matter halo core radius

# === Relative geometry ===
distance_kpc = 780           # current MW-M31 distance
v_radial_kms = -110          # approaching (negative = toward MW)
v_transverse_kms = 17        # tangential velocity

# === Andromeda disc orientation ===
# Inclination 77 deg to line of sight (nearly edge-on from Earth)
# Position angle ~35 deg on sky
inclination_deg = 77
position_angle_deg = 35

# === Convert to code units ===
print("=" * 60)
print("MILKY WAY - ANDROMEDA PARAMETERS (code units)")
print("=" * 60)

print("\n--- Milky Way ---")
mw_R = mw_radius_kpc / du
mw_Ri = mw_inner_kpc / du
mw_M_central = mw_bulge_msun / mu
mw_Mfrac = mw_disc_msun / mw_bulge_msun
mw_haloVc = mw_vc_kms / vu
mw_haloRc = mw_haloRc_kpc / du
mw_particle_msun = mw_disc_msun / 39999

print(f"  Disc radius:     {mw_R:.1f} code units ({mw_radius_kpc} kpc)")
print(f"  Inner radius:    {mw_Ri:.1f} code units ({mw_inner_kpc} kpc)")
print(f"  Central mass:    {mw_M_central:.1f} code units ({mw_bulge_msun:.1e} Msun)")
print(f"  Mass fraction:   {mw_Mfrac:.2f}")
print(f"  Particle mass:   ~{mw_particle_msun:.0f} Msun each")
print(f"  Halo Vc:         {mw_haloVc:.1f} code units ({mw_vc_kms} km/s)")
print(f"  Halo Rc:         {mw_haloRc:.1f} code units ({mw_haloRc_kpc} kpc)")
print(f"  Disc normal:     (0, 1, 0)")

print("\n--- Andromeda ---")
and_R = and_radius_kpc / du
and_Ri = and_inner_kpc / du
and_M_central = and_bulge_msun / mu
and_Mfrac = and_disc_msun / and_bulge_msun
and_haloVc = and_vc_kms / vu
and_haloRc = and_haloRc_kpc / du
and_particle_msun = and_disc_msun / 39999

print(f"  Disc radius:     {and_R:.1f} code units ({and_radius_kpc} kpc)")
print(f"  Inner radius:    {and_Ri:.1f} code units ({and_inner_kpc} kpc)")
print(f"  Central mass:    {and_M_central:.1f} code units ({and_bulge_msun:.1e} Msun)")
print(f"  Mass fraction:   {and_Mfrac:.2f}")
print(f"  Particle mass:   ~{and_particle_msun:.0f} Msun each")
print(f"  Halo Vc:         {and_haloVc:.1f} code units ({and_vc_kms} km/s)")
print(f"  Halo Rc:         {and_haloRc:.1f} code units ({and_haloRc_kpc} kpc)")

# Compute Andromeda disc normal
inc = math.radians(inclination_deg)
pa = math.radians(position_angle_deg)
nx = math.cos(inc)
ny = math.sin(inc) * math.cos(pa)
nz = math.sin(inc) * math.sin(pa)
print(f"  Disc normal:     ({nx:.4f}, {ny:.4f}, {nz:.4f})")
print(f"    (inclination={inclination_deg} deg, PA={position_angle_deg} deg)")

print("\n--- Separation & Velocity ---")
sep = distance_kpc / du
print(f"  Distance:        {sep:.1f} code units ({distance_kpc} kpc)")
print(f"  Radial vel:      {v_radial_kms:.1f} code units ({v_radial_kms} km/s, approaching)")
print(f"  Transverse vel:  {v_transverse_kms:.1f} code units ({v_transverse_kms} km/s)")

print("\n--- Timing ---")
# Time to collision (rough): distance / approach speed
t_collision_myr = distance_kpc / (abs(v_radial_kms) * 1.022)  # 1 km/s ~ 1.022 kpc/Gyr
t_collision_code = t_collision_myr * 1000 / 58.7  # Myr to code units
print(f"  Estimated collision time: ~{t_collision_myr:.0f} Gyr = ~{t_collision_code:.0f} code time units")
print(f"  (Very long! For faster interaction, reduce initial separation)")

print("\n--- Script lines ---")
print(f"N_SystemBodies  40000  40000")
print(f"Camera          0.0  8000.0  8000.0")
print(f"# Milky Way at origin, disc in x-z plane")
print(f"GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   {mw_M_central:.1f} {mw_Mfrac:.2f} {mw_R:.1f} {mw_Ri:.1f} 0.1  {mw_haloVc:.1f} {mw_haloRc:.1f}")
print(f"# Andromeda along +x axis")
print(f"GalaxyDisc  1   {sep:.1f} 0.0 0.0   {v_radial_kms:.1f} {v_transverse_kms:.1f} 0.0   {nx:.4f} {ny:.4f} {nz:.4f}   {and_M_central:.1f} {and_Mfrac:.2f} {and_R:.1f} {and_Ri:.1f} 0.1  {and_haloVc:.1f} {and_haloRc:.1f}")

# Rotation directions:
print("\n--- Rotation directions ---")
print("  MW: CCW when viewed from +y (disc normal direction)")
print("  This matches the real MW rotation (clockwise from south galactic pole)")
print("  Andromeda: CCW relative to its disc normal")
