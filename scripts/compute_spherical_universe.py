"""
Compute initial conditions for a Spherical Universe script that produces
filamentary structure formation with voids.

=== Universe Classification ===

This simulation represents a CLOSED UNIVERSE with weak dark energy:
  - H/H_crit = 0.8 (sub-critical expansion rate)
  - KE/|PE| = 0.64 (total energy is negative -> gravitationally bound)
  - Without dark energy: expands, turns around, recollapses (Big Crunch)
  - With FDE=0.10: recollapse is significantly delayed but not prevented;
    the system remains formally bound. Dark energy is too weak to fully
    unbind the system, but it substantially extends the lifetime of the
    cosmic web before eventual recollapse.

=== Approach: Hubble Flow + Gravitational Instability ===

In real cosmology, the universe starts in a nearly uniform expanding state
(the Hubble flow: v = H*r). Structure forms because overdense regions
decelerate faster than the mean, eventually turning around and collapsing,
while underdense regions decelerate less and become voids.

This is implemented by giving each particle an initial velocity proportional
to its distance from the center: v_i = H * r_i (radially outward). Dark
energy (FDE) is set small enough to be negligible during structure formation
but provides gentle support against late-time recollapse.

=== Expected Simulation Timeline ===

  t ~ 0-0.3 code units (0-18 Myr):
    Sphere expands uniformly. Density perturbations begin growing via
    gravitational instability. No visible structure yet.

  t ~ 0.3-1.0 code units (18-59 Myr):
    Filaments become visible as overdense streaks collapse. Voids open
    between them. The cosmic web takes shape.

  t ~ 1-3 code units (59-176 Myr):
    Well-developed filamentary structure with nodes (intersections of
    filaments), filaments (linear collapsed structures), and voids
    (expanding underdense regions). This is the peak of structural
    complexity.

  t > 3 code units (> 176 Myr):
    Expansion slows as gravity decelerates the outer shells. FDE provides
    gentle outward support. The cosmic web persists but the bulk system
    will eventually recollapse on a much longer timescale than without FDE.

Note: The absolute timescales (Myr) are set by the code unit system but
should not be taken as cosmologically realistic. The simulation captures
the QUALITATIVE progression of structure formation, not the quantitative
timeline of the real universe (where structure formation spans ~10 Gyr).

=== Why Hubble Flow (not random velocities + FDE) ===

A previous approach using random velocities and strong FDE failed because:
  1. FDE pushes radially from origin: center (r~0) feels no force and
     collapses into a central blob while the exterior evaporates.
  2. Without coherent outward expansion, there is no differential
     deceleration between over/underdense regions — no filaments form.
  3. Random velocities provide thermal pressure that prevents small-scale
     collapse without enabling large-scale structure.

Hubble flow solves all three: every point expands proportionally, overdense
patches decelerate faster (eventually collapsing), underdense patches
decelerate slower (becoming voids). This is identical in mechanism to real
cosmological structure formation.

=== Unit system (G=1) ===
  1 distance unit = 60 pc = 0.060 kpc
  1 mass unit     = 10^4 solar masses
  1 velocity unit = 1 km/s
  1 time unit     = 58.7 Myr  (= 60 pc / 1 km/s)
"""

import math

# === System parameters ===
R = 200.0           # sphere radius (code units) = 12 kpc
N = 100000          # number of particles
M_total = 5.0e7     # total mass (code units) = 5e11 Msun

# === Derived quantities ===
rho = M_total / ((4.0 / 3.0) * math.pi * R**3)
spacing = R * ((4.0 / 3.0) * math.pi / N) ** (1.0 / 3.0)
v_virial = math.sqrt(M_total / R)
v_escape = math.sqrt(2) * v_virial
t_freefall = math.sqrt(R**3 / M_total) * math.pi / 2.0

print("=" * 60)
print("FILAMENTARY STRUCTURE FORMATION PARAMETERS")
print("  Universe type: CLOSED (with weak dark energy)")
print("=" * 60)

print(f"\n--- Base quantities ---")
print(f"  Sphere radius:       {R:.0f} code units = {R*0.06:.1f} kpc")
print(f"  Total mass:          {M_total:.2e} code units = {M_total*1e4:.2e} Msun")
print(f"  Particle count:      {N}")
print(f"  Particle mass:       {M_total/N:.1f} code units = {M_total/N*1e4:.2e} Msun")
print(f"  Mean density:        {rho:.4f} code units")
print(f"  Inter-particle spacing: {spacing:.1f} code units")
print(f"  Gaussian position scatter: sigma = 10 code units (hardcoded)")
print(f"  Scatter / spacing:   {10/spacing:.2f} (>1 = coherent overdensities)")

print(f"\n--- Velocity scales ---")
print(f"  Virial velocity:     {v_virial:.1f} km/s")
print(f"  Escape velocity:     {v_escape:.1f} km/s")
print(f"  Free-fall time:      {t_freefall:.2f} code units = {t_freefall*58.7:.1f} Myr")

# === Critical Hubble parameter ===
# For a uniform sphere with Hubble flow v = H*r:
#
# Total KE = integral_0^R (1/2) * rho * (H*r)^2 * 4*pi*r^2 dr
#          = (1/2) * rho * H^2 * (4*pi/5) * R^5
#          = (3/10) * M * H^2 * R^2
#
# Total PE = -(3/5) * G * M^2 / R  (gravitational self-energy of uniform sphere)
#
# Critical condition (KE = |PE|):
#   (3/10) * M * H^2 * R^2 = (3/5) * M^2 / R     [G=1]
#   H^2 = 2*M / R^3
#   H_crit = sqrt(2*M/R^3)
#
# H > H_crit: open (unbound, expands forever, no structure forms)
# H = H_crit: flat (marginally unbound, Einstein-de Sitter analog)
# H < H_crit: closed (bound, recollapses unless dark energy intervenes)

H_crit = math.sqrt(2.0 * M_total / R**3)
v_edge_crit = H_crit * R

print(f"\n--- Critical Hubble parameter ---")
print(f"  H_crit = sqrt(2*M/R^3) = {H_crit:.4f}")
print(f"  v at edge (H_crit * R) = {v_edge_crit:.1f} km/s")
print(f"  v_escape = {v_escape:.1f} km/s")

# === Chosen Hubble parameter ===
H_frac = 0.8
H_use = H_frac * H_crit
v_edge = H_use * R
KE_frac = H_frac**2  # KE/|PE| = (H/H_crit)^2
total_energy_sign = "NEGATIVE (bound)" if KE_frac < 1 else "POSITIVE (unbound)"

print(f"\n--- Chosen Hubble parameter ---")
print(f"  H/H_crit = {H_frac}")
print(f"  H = {H_use:.4f}")
print(f"  v at edge = {v_edge:.1f} km/s")
print(f"  v at R/2 = {H_use*R/2:.1f} km/s")
print(f"  KE/|PE| = {KE_frac:.2f}")
print(f"  Total energy: {total_energy_sign}")
print(f"  Universe classification: {'CLOSED' if KE_frac < 1 else 'OPEN' if KE_frac > 1 else 'FLAT'}")

# === Dark energy (FDE) ===
# FDE_balance = M/R^3 = (4/3)*pi*rho: where dark energy exactly balances mean gravity
# For FDE < FDE_balance: gravity still dominates at the initial scale
# The "dominance ratio" = a_FDE / a_grav at radius r (with all mass inside):
#   a_FDE = FDE * r
#   a_grav = G*M / r^2 = M / r^2   [G=1]
#   ratio = FDE * r^3 / M
# FDE dominates when r > (M/FDE)^(1/3)

FDE_balance = M_total / R**3
FDE_use = 0.10  # chosen empirically: negligible during structure formation,
                # delays recollapse but does not prevent it (closed universe)

r_crossover = (M_total / FDE_use) ** (1.0 / 3.0)
ratio_at_100 = FDE_use * 100**3 / M_total
ratio_at_200 = FDE_use * 200**3 / M_total
ratio_at_300 = FDE_use * 300**3 / M_total

print(f"\n--- Dark energy (FDE) ---")
print(f"  FDE_balance (= M/R^3):  {FDE_balance:.4f}")
print(f"  FDE chosen:             {FDE_use:.2f}")
print(f"  FDE / FDE_balance:      {FDE_use/FDE_balance:.3f}")
print(f"  Crossover radius:       {r_crossover:.0f} code units")
print(f"    (FDE dominates gravity for r > {r_crossover:.0f})")
print(f"  Dominance ratio at r=100: {ratio_at_100:.4f} ({ratio_at_100*100:.2f}%)")
print(f"  Dominance ratio at r=200: {ratio_at_200:.4f} ({ratio_at_200*100:.2f}%)")
print(f"  Dominance ratio at r=300: {ratio_at_300:.4f} ({ratio_at_300*100:.2f}%)")
print(f"  -> FDE is negligible during structure formation (r ~ 100-300)")
print(f"  -> FDE delays recollapse once sphere expands beyond r ~ {r_crossover:.0f}")

# === Jeans analysis ===
# The Jeans length is the minimum scale at which gravity overcomes thermal
# pressure (random motions) to enable collapse.
#
# In the Hubble flow context, the effective "thermal" velocity is set by
# the velocity dispersion arising from the Gaussian position scatter:
# A particle displaced by delta_r from its "ideal" position has a velocity
# error of ~H * delta_r relative to its neighbors.
#
# lambda_J ~ v_thermal / sqrt(G * rho)  [G=1]
# Structures larger than lambda_J collapse; smaller ones are pressure-supported.

v_thermal = H_use * 10.0  # velocity dispersion from sigma=10 position scatter
lambda_J = v_thermal * math.sqrt(math.pi / rho)
N_jeans = (4.0 / 3.0) * math.pi * (lambda_J / 2.0)**3 * rho / (M_total / N)
t_jeans = 1.0 / math.sqrt(rho)

print(f"\n--- Jeans analysis ---")
print(f"  Effective thermal velocity: H * sigma = {v_thermal:.1f} km/s")
print(f"  Jeans length: {lambda_J:.1f} code units")
print(f"  Particles per Jeans volume: {N_jeans:.0f}")
print(f"  Jeans collapse time: {t_jeans:.2f} code units = {t_jeans*58.7:.1f} Myr")
print(f"  -> Perturbations larger than {lambda_J:.0f} code units collapse")
print(f"  -> Perturbations smaller than {lambda_J:.0f} are pressure-supported")

# === Timescales ===
# Turnaround: outer shell reaches max radius then falls back
# Without FDE:
r_max_no_fde = R / (1.0 - KE_frac)  # from energy: E = -M/r_max -> r_max = M/|E|_per_unit
t_turnaround_no_fde = math.pi * math.sqrt(r_max_no_fde**3 / (2.0 * M_total))

print(f"\n--- Timescales ---")
print(f"  Free-fall time (whole sphere):  {t_freefall:.2f} code units")
print(f"  Jeans collapse time:            {t_jeans:.2f} code units")
print(f"  Turnaround time (no FDE):       {t_turnaround_no_fde:.1f} code units")
print(f"  Turnaround radius (no FDE):     {r_max_no_fde:.0f} code units")
print(f"  Structure visible at:           t ~ {t_jeans:.1f}-{3*t_jeans:.1f} code units")
print(f"  With FDE=0.10: turnaround delayed well beyond this")

# === Particle count ===
print(f"\n--- Particle count effect ---")
print(f"  Gaussian scatter sigma = 10 code units (hardcoded in generator)")
print(f"  For scatter to seed coherent overdensities: scatter > spacing")
for N_t in [20000, 40000, 100000, 200000]:
    sp = R * ((4.0/3.0)*math.pi/N_t)**(1.0/3.0)
    N_J_t = (4.0/3.0)*math.pi*(lambda_J/2)**3 * rho / (M_total/N_t)
    print(f"    N={N_t:6d}: spacing={sp:.1f}, scatter/spacing={10/sp:.2f}, particles/Jeans_vol={N_J_t:.0f}")
print(f"  N=100,000: good resolution (scatter/spacing=1.44, 108 particles per Jeans vol)")

# === Output ===
print(f"\n{'='*60}")
print(f"SCRIPT PARAMETERS")
print(f"{'='*60}")
print(f"  G               1.0")
print(f"  FDE             {FDE_use}")
print(f"  dt              0.001")
print(f"  r_soft          0.5")
print(f"  BH_Opening_Theta  0.7")
print(f"  N_SystemBodies  {N}")
print(f"  Camera          0.0  960.0  0.0")
print(f"  SphericalUniverse  0   0.0 0.0 0.0   0.0 0.0 0.0   {M_total:.1e}  {R:.1f}  {H_use:.2f}  0.0  1.0")
print(f"")
print(f"  H = {H_use:.2f}: Hubble parameter (each particle gets v = H*r, radial outward)")
print(f"  H_crit = {H_crit:.2f}, H/H_crit = {H_frac} -> closed universe")
print(f"  Camera on +Y axis at 960 (= R*4.8), looking down at origin")
