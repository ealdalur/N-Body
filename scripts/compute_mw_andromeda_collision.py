"""
Compute initial conditions for MW-Andromeda collision scenario.

Places the galaxies close enough (~150 kpc apart) with an estimated
infall velocity so the collision begins within a reasonable simulation time.

Uses energy conservation to estimate approach velocity at close range
given the starting conditions at 780 kpc with 110 km/s radial velocity.
"""

import math

# Unit system
du = 0.060   # 1 code unit = 60 pc = 0.060 kpc
mu = 1.0e4   # 1 code mass unit = 10^4 Msun
vu = 1.0     # 1 code velocity unit = 1 km/s

# Galaxy parameters (same as Milky_Way_Andromeda.sim)
mw_M = 1500000.0       # MW central body (code units)
mw_Mfrac = 3.0
mw_R = 446.7
and_M = 3000000.0      # Andromeda central body (code units)
and_Mfrac = 2.33
and_R = 558.3

# Halo parameters
mw_haloVc = 220.0      # km/s
and_haloVc = 225.0

# Starting conditions (current real observations)
r0 = 13000.0           # current separation (780 kpc)
v0 = 110.0             # current approach velocity (km/s)

# Desired close separation
r_close = 2500.0       # 150 kpc in code units

# Estimate total enclosed mass at close range using isothermal halo model:
# M_halo(r) ~ v_c^2 * r / G (with G=1)
mw_halo_enclosed = mw_haloVc**2 * r_close
and_halo_enclosed = and_haloVc**2 * r_close
mw_baryonic = mw_M * (1 + mw_Mfrac)
and_baryonic = and_M * (1 + and_Mfrac)
M_total = mw_halo_enclosed + and_halo_enclosed + mw_baryonic + and_baryonic

print("=== Mass estimates ===")
print(f"  MW halo enclosed at {r_close*du:.0f} kpc: {mw_halo_enclosed:.2e} code = {mw_halo_enclosed*mu:.2e} Msun")
print(f"  And halo enclosed at {r_close*du:.0f} kpc: {and_halo_enclosed:.2e} code = {and_halo_enclosed*mu:.2e} Msun")
print(f"  MW baryonic: {mw_baryonic:.2e} code = {mw_baryonic*mu:.2e} Msun")
print(f"  And baryonic: {and_baryonic:.2e} code = {and_baryonic*mu:.2e} Msun")
print(f"  Total: {M_total:.2e} code = {M_total*mu:.2e} Msun")

# Energy conservation: 0.5*v0^2 - G*M/r0 = 0.5*v^2 - G*M/r
# v^2 = v0^2 + 2*G*M*(1/r - 1/r0)  [G=1]
v_sq = v0**2 + 2 * 1.0 * M_total * (1.0/r_close - 1.0/r0)
v_approach = math.sqrt(v_sq)

print(f"\n=== Velocity at close range ===")
print(f"  v^2 = {v0}^2 + 2*{M_total:.2e}*(1/{r_close} - 1/{r0})")
print(f"  v^2 = {v0**2} + {2*M_total*(1/r_close - 1/r0):.0f}")
print(f"  v = {v_approach:.1f} km/s")

# Add transverse velocity for off-center collision (more realistic and visually interesting)
v_transverse = 50.0  # km/s

print(f"\n=== Collision scenario ===")
sep_kpc = r_close * du
print(f"  Separation: {r_close:.0f} code units = {sep_kpc:.0f} kpc")
print(f"  Approach velocity: {v_approach:.0f} km/s (from energy conservation)")
print(f"  Transverse velocity: {v_transverse:.0f} km/s (for grazing encounter)")
print(f"  Gap between disc edges: {r_close - mw_R - and_R:.0f} code units = {(r_close-mw_R-and_R)*du:.0f} kpc")
print(f"  Time to first overlap (no extra accel): {(r_close-mw_R-and_R)/v_approach:.0f} time units = {(r_close-mw_R-and_R)/v_approach*58.7:.0f} Myr")

# Camera: midpoint between galaxies, elevated
cam_x = r_close / 2.0
cam_y = 2000.0
print(f"\n=== Script values ===")
print(f"  Camera: {cam_x:.1f} {cam_y:.1f} 0.0")
print(f"  MW: pos=(0,0,0) vel=(0,0,0)")
print(f"  And: pos=({r_close:.1f}, 0, 0) vel=(-{v_approach:.1f}, {v_transverse:.1f}, 0)")

# Compare different separations for choosing best scenario
print(f"\n=== Separation comparison ===")
print(f"{'Sep (code)':>10s} {'Sep (kpc)':>9s} {'v (km/s)':>9s} {'Gap':>6s} {'t_contact':>10s} {'t (Myr)':>8s}")
print("-" * 60)

for r in [5000, 4000, 3500, 3000, 2500]:
    # Two-body problem: each galaxy's halo contributes to the potential
    # At separation r, approximate total gravitating mass as sum of
    # both halos enclosed within r (from their respective centers)
    # plus baryonic mass
    M_eff = (mw_haloVc**2 + and_haloVc**2) * r * 0.5 + mw_baryonic + and_baryonic
    v_sq_r = v0**2 + 2 * M_eff * (1.0/r - 1.0/r0)
    if v_sq_r < 0:
        v_sq_r = v0**2
    v_r = math.sqrt(v_sq_r)
    gap = r - mw_R - and_R
    t_contact = gap / v_r if v_r > 0 else 9999
    print(f"{r:10d} {r*du:9.0f} {v_r:9.1f} {gap:6.0f} {t_contact:10.0f} {t_contact*58.7:8.0f}")
