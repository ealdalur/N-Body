"""
Compute initial conditions for an M51 (Whirlpool Galaxy) interaction simulation.

The M51 system consists of NGC 5194 (M51a, a grand-design Sbc spiral) and
NGC 5195 (M51b, a compact SB0-pec lenticular companion). M51b has already
crossed M51a's disc plane at least once, exciting the prominent two-arm spiral
structure that makes this system famous.

This script computes pre-encounter initial conditions that will produce
M51-like tidal spiral arms when the simulation is run forward.

=== Key Observational Constraints ===

All orbital parameters below are taken from Salo & Laurikainen (2000), Table 2
and the Fig. 1 caption, for their BOUND multiple-encounter model (the model
whose morphology and kinematics they find best matches M51). Angular sizes are
converted at the paper's own assumed distance of 9.6 Mpc, so that its
dimensionless ratios carry over unchanged.

  - Encounter is a BOUND multiple-passage orbit (not a single parabolic flyby)
  - Orbital inclination iorb = 75-85 deg -- a near-POLAR orbit
  - Eccentricity e ~ 0.2 for the current (bound) orbit
  - Disc-plane crossings at Rcross = 1.2-1.4 and Rdown = 1.2-1.3 primary Rd,
    i.e. ~22-26 kpc; apocenter lies between the two crossings
  - Principal (spiral-inducing) crossing ~400-500 Myr ago, in the south
  - Most recent crossing 50-100 Myr ago
  - Mass ratio Mp = M_tot(comp)/M_tot(prim) within 4Re = 0.5-0.7 (nominal 0.55)
  - Toomre Q_T = 1.5 standard value for both discs
  - M_disc/M_tot = 1/3 within 4Re (so M_halo/M_disc ~ 2)
  - Primary disc:   Re = 100 arcsec, Rd = 4Re   = 400 arcsec; halo Rc =  8 arcsec
  - Companion disc: Re =  33 arcsec, Rd = 7.3Re = 240 arcsec; halo Rc = 40 arcsec
  - Companion disc tilted 32.5 deg relative to the primary disc

Note the paper parameterizes by DISC-PLANE CROSSING distance rather than
pericenter, "motivated by the high inclination of the relative orbit".

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
  - Sofue et al. 1999; Meidt et al. 2013 (rotation curve)
  - Toomre & Toomre 1972 (original tidal interaction model)
"""

import math

# === Unit system ===
du = 0.060   # 1 code distance unit = 60 pc = 0.060 kpc
mu = 1.0e4   # 1 code mass unit = 10^4 Msun
vu = 1.0     # 1 code velocity unit = 1 km/s
tu = 58.7    # 1 code time unit = 60 pc / (1 km/s) = 58.7 Myr

# === Distance ===
# The paper's own assumed distance, used throughout so that every angular size
# and every dimensionless ratio it quotes translates consistently. The modern
# TRGB value is 8.58 Mpc (McQuinn et al. 2016) and is the better measurement of
# the real galaxy, but it is the wrong choice here: Rcross is quoted in units of
# the primary Rd, and the paper's 1.2-1.4 range refers to ITS Rd of 18.6 kpc.
distance_Mpc = 9.6                  # Scoville & Young 1983, as adopted by S&L
arcsec_to_kpc = distance_Mpc * 1.0e3 / 206265.0

# === NGC 5194 (M51a) — Grand-design spiral ===
# Salo & Laurikainen (2000) section 2.2, verbatim:
#   Re = 100 arcsec (exponential scalelength)
#   Rd = 4*Re = 400 arcsec (disc truncation; also their unit of length)
#   halo core Rc = 8 arcsec
m51a_Re_arcsec = 100.0
m51a_Rd_arcsec = 400.0
m51a_haloRc_arcsec = 8.0

m51a_Vtotal_kms = 210               # observed flat rotation velocity
m51a_radius_kpc = m51a_Rd_arcsec * arcsec_to_kpc      # Rd, the truncation
m51a_h_r_kpc = m51a_Re_arcsec * arcsec_to_kpc         # Re, the scale length
m51a_inner_kpc = 0.3                # inner bulge region (no disc particles)
m51a_haloRc_kpc = m51a_haloRc_arcsec * arcsec_to_kpc

# Rotation curve decomposition
m51a_Vbaryon = m51a_Vtotal_kms / math.sqrt(3)
m51a_Vhalo = m51a_Vtotal_kms * math.sqrt(2.0/3.0)

# Baryonic mass from V_baryon^2 * R / G (G=1 in code units)
m51a_R_code = m51a_radius_kpc / du
m51a_h_r_code = m51a_h_r_kpc / du
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
# Salo & Laurikainen (2000) section 2.2: the companion disc is taken to be 70%
# of the primary's extent, Rd = 240 arcsec ~ 7*Re with Re = 33 arcsec, and its
# halo has Rc = 40 arcsec ("lesser degree of halo concentration expected for
# smaller galaxies").
#
# Crucially, the paper does NOT set the companion mass from an assumed rotation
# curve. It sets the total mass RATIO Mp = M_tot(comp) / M_tot(prim) within 4Re,
# with a nominal Mp = 0.55 (Table 1 caption; Table 2 range 0.5-0.7). The paper
# notes this exceeds the observed disc mass ratio (~0.4) deliberately, because
# "small galaxies are likely to have more dominant halo components".
#
# Their explicit constraint: "Passages with Mp < 0.55 ... seem to be too weak to
# account for the observations", while "for Mp > 0.7 the velocities in the north
# start to be somewhat too large".
# NOTE the companion's Rd/Re ratio is 240/33 = 7.3, NOT 4 like the primary, so
# its scale length must be passed explicitly to the simulation rather than left
# to the default Rd/4.
m51b_Re_arcsec = 33.0
m51b_Rd_arcsec = 240.0
m51b_haloRc_arcsec = 40.0

m51b_radius_kpc = m51b_Rd_arcsec * arcsec_to_kpc      # Rd, the truncation
m51b_h_r_kpc = m51b_Re_arcsec * arcsec_to_kpc         # Re, the scale length
m51b_inner_kpc = 0.2                # inner region
m51b_haloRc_kpc = m51b_haloRc_arcsec * arcsec_to_kpc
m51b_mass_ratio = 0.55              # Mp, paper nominal (Table 1 caption)

m51b_R_code = m51b_radius_kpc / du
m51b_h_r_code = m51b_h_r_kpc / du
m51b_haloRc_code = m51b_haloRc_kpc / du

# Companion total mass from the paper's mass ratio Mp, applied to the primary's
# total (baryon + halo) mass within its disc. This replaces the previous
# approach of assuming a 130 km/s rotation curve for M51b, which produced
# Mp ~ 0.30 -- below the paper's stated lower bound of 0.55.
m51a_Mhalo_disc_pre = (m51a_haloVc**2 * m51a_R_code**3
                       / (m51a_R_code**2 + m51a_haloRc_code**2))
m51a_total_code = m51a_baryon_code + m51a_Mhalo_disc_pre
m51b_total_code = m51b_mass_ratio * m51a_total_code

# Split the companion total into baryon + halo on the same 1/3 : 2/3 basis used
# for the primary, keeping M_halo/M_baryon ~ 2 (the paper holds the halo-to-disc
# ratio fixed at 2:1 within the disc radius for both components).
m51b_baryon_code = m51b_total_code / 3.0
m51b_baryon_msun = m51b_baryon_code * mu
m51b_Mhalo_target = m51b_total_code - m51b_baryon_code

# haloVc that puts M_halo_target inside the companion disc radius, inverting
# M_halo(R) = Vc^2 * R^3 / (R^2 + Rc^2)
m51b_haloVc = math.sqrt(m51b_Mhalo_target * (m51b_R_code**2 + m51b_haloRc_code**2)
                        / m51b_R_code**3)

# Implied rotation curve of the companion (diagnostic, no longer an input)
m51b_Vbaryon = math.sqrt(m51b_baryon_code / m51b_R_code)
m51b_Vhalo = m51b_haloVc * m51b_R_code / math.sqrt(m51b_R_code**2 + m51b_haloRc_code**2)
m51b_Vtotal_kms = math.sqrt(m51b_Vbaryon**2 + m51b_Vhalo**2)

# Bulge/disc decomposition: bulge ~ 40% for SB0 lenticular
m51b_bulge_frac = 0.40
m51b_bulge_code = round(m51b_baryon_code * m51b_bulge_frac / 10000) * 10000
m51b_disc_code = m51b_baryon_code - m51b_bulge_code
m51b_Mfrac = m51b_disc_code / m51b_bulge_code

# === Interaction geometry (Salo & Laurikainen 2000, Table 2 + Fig. 1 caption) ===
#
# PERICENTER / CROSSING DISTANCE. The paper parameterizes by disc-plane crossing
# distance, not pericenter, "motivated by the high inclination of the relative
# orbit". Rcross = 1.2-1.4 and Rdown = 1.2-1.3, in units of the primary Rd.
#   Abstract, rescaled from 9.6 to 8.58 Mpc: main crossing 22.3-26.8 kpc,
#   latest crossing 17.9-22.3 kpc.
#
# ECCENTRICITY. Fig. 1 caption, bound model: "for the last case e = 0.2 and
# iorb = 85 deg". Their other (single-passage) models have e = 0.67-0.83; we
# reproduce the BOUND model, so e = 0.2.
#
# APOCENTER is DERIVED from pericenter and e, not chosen independently:
#   r_apo = r_peri * (1+e)/(1-e) = 1.5 * r_peri   for e = 0.2
#
# INCLINATION. Table 2: iorb = 75-85 deg; Fig. 1 caption gives 85 deg for the
# bound model. Near-POLAR. The paper treats the high inclination as essential --
# it produces the out-of-plane velocities, the 40-50 deg tilted tail, and the
# S-shaped major-axis rotation curve.
#
# NODE GEOMETRY (argument of pericenter). Fig. 1 caption, bound model: "the
# apocentre is between the two disc crossings". So apocenter must sit at maximum
# height ABOVE the disc plane, with the line of nodes ~90 deg away in orbital
# phase -- i.e. argument of pericenter = 90 deg.
#
# This is a separate degree of freedom from pericenter, apocenter, e and iorb,
# and it was previously wrong: apocenter was placed IN the disc plane (y=0),
# which put the nodes near the apsides instead of between them. The companion
# then started on the plane, swung away to -y, and returned through the plane
# from below -- not the paper's geometry.
# TARGET CROSSING DISTANCE.
#
# Rcross is the paper's actual constraint (Table 2: 1.2-1.4 in units of the
# primary Rd), so it is the INPUT and the pericenter is solved to match it.
# Pericenter is not independently meaningful in the paper's framework.
#
# Set to 1.20, the LOW end of the range and therefore the strongest perturbation
# the paper allows: "passages with Mp < 0.55 or Rcross > 1.4 seem to be too weak
# to account for the observations".
#   Rcross = 1.20 -> S = Mp*(Rd/r)^3 = 0.318
#   Rcross = 1.30 -> S = 0.250
#   Rcross = 1.40 -> S = 0.200  (the paper's stated weak limit)
# Toomre Q for both discs: the paper's standard value. Gives a warm disc, stable
# against noise-driven multi-arm structure while still responsive to the m=2
# tidal forcing.
toomre_Q = 1.5                      # Salo & Laurikainen section 2.2

target_Rcross = 1.20                # units of primary Rd; paper range 1.2-1.4
orbital_eccentricity = 0.2          # bound model (Fig. 1 caption)
orbital_inclination_deg = 80.0      # Table 2 range 75-85, mid-range

# pericenter_kpc and apocenter_kpc are solved below, once the halo constants
# needed for the orbit integration are available (see "Orbit solve" section).

# ============================================================
print("=" * 65)
print("M51 WHIRLPOOL GALAXY INTERACTION PARAMETERS")
print("=" * 65)

# === Print M51a parameters ===
print("\n--- NGC 5194 (M51a) ---")
print(f"  Observed V_total:  {m51a_Vtotal_kms} km/s")
print(f"  V_baryon:          {m51a_Vbaryon:.1f} km/s (1/3 of V^2)")
print(f"  V_halo:            {m51a_Vhalo:.1f} km/s (2/3 of V^2)")
print(f"  Disc radius Rd:    {m51a_R_code:.1f} code units ({m51a_radius_kpc:.2f} kpc) = {m51a_Rd_arcsec:.0f} arcsec")
print(f"  Scale length Re:   {m51a_h_r_code:.1f} code units ({m51a_h_r_kpc:.2f} kpc) = {m51a_Re_arcsec:.0f} arcsec")
print(f"  Rd / Re:           {m51a_Rd_arcsec/m51a_Re_arcsec:.2f}")
print(f"  Inner radius:      {m51a_inner_kpc/du:.1f} code units ({m51a_inner_kpc} kpc)")
print(f"  Bulge mass (M):    {m51a_bulge_code:.0f} code units ({m51a_bulge_code*mu:.2e} Msun)")
print(f"  Disc mass:         {m51a_disc_code:.0f} code units ({m51a_disc_code*mu:.2e} Msun)")
print(f"  Total baryonic:    {m51a_baryon_code:.0f} code units ({m51a_baryon_msun:.2e} Msun)")
print(f"  Mfrac:             {m51a_Mfrac:.2f}")
print(f"  Halo Vc:           {m51a_haloVc:.1f} km/s")
print(f"  Halo Rc:           {m51a_haloRc_code:.1f} code units ({m51a_haloRc_kpc:.3f} kpc)")

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
print(f"  Implied V_total:   {m51b_Vtotal_kms:.1f} km/s (from Mp, not assumed)")
print(f"  V_baryon:          {m51b_Vbaryon:.1f} km/s")
print(f"  V_halo:            {m51b_Vhalo:.1f} km/s")
print(f"  Disc radius Rd:    {m51b_R_code:.1f} code units ({m51b_radius_kpc:.2f} kpc) = {m51b_Rd_arcsec:.0f} arcsec")
print(f"  Scale length Re:   {m51b_h_r_code:.1f} code units ({m51b_h_r_kpc:.2f} kpc) = {m51b_Re_arcsec:.0f} arcsec")
print(f"  Rd / Re:           {m51b_Rd_arcsec/m51b_Re_arcsec:.2f}  (NOT 4 -- must be passed explicitly)")
print(f"  Re(comp)/Re(prim): {m51b_Re_arcsec/m51a_Re_arcsec:.3f}  (paper 33/100; was 0.60 when h_r=Rd/4)")
print(f"  Inner radius:      {m51b_inner_kpc/du:.1f} code units ({m51b_inner_kpc} kpc)")
print(f"  Bulge mass (M):    {m51b_bulge_code:.0f} code units ({m51b_bulge_code*mu:.2e} Msun)")
print(f"  Disc mass:         {m51b_disc_code:.0f} code units ({m51b_disc_code*mu:.2e} Msun)")
print(f"  Total baryonic:    {m51b_baryon_code:.0f} code units ({m51b_baryon_msun:.2e} Msun)")
print(f"  Mfrac:             {m51b_Mfrac:.2f}")
print(f"  Halo Vc:           {m51b_haloVc:.1f} km/s")
print(f"  Halo Rc:           {m51b_haloRc_code:.1f} code units ({m51b_haloRc_kpc:.3f} kpc)")

v_b_check = math.sqrt(m51b_baryon_code / m51b_R_code)
v_h_check = m51b_haloVc * m51b_R_code / math.sqrt(m51b_R_code**2 + m51b_haloRc_code**2)
v_total_check = math.sqrt(v_b_check**2 + v_h_check**2)
print(f"  V_total consistency:  {v_total_check:.1f} km/s (recomputed from masses)")

m51b_Mhalo_disc = m51b_haloVc**2 * m51b_R_code**3 / (m51b_R_code**2 + m51b_haloRc_code**2)
print(f"  M_halo within disc: {m51b_Mhalo_disc:.0f} code units")
print(f"  M_halo/M_baryon:   {m51b_Mhalo_disc/m51b_baryon_code:.2f} (target: ~2.0)")

# === Mass ratio Mp (the paper's primary mass constraint) ===
# The paper defines Mp = M_tot(companion) / M_tot(primary) within 4Re, where
# 4Re = Rd is the disc truncation radius. This is the quantity constrained to
# 0.5-0.7, so it is the one to verify.
Mp_check = m51b_total_code / m51a_total_code

print(f"\n--- Mass ratio Mp (companion/primary total within 4Re = Rd) ---")
print(f"  M51a total within Rd: {m51a_total_code:.0f} "
      f"(baryon {m51a_baryon_code:.0f} + halo {m51a_Mhalo_disc_pre:.0f})")
print(f"  M51b total within Rd: {m51b_total_code:.0f} "
      f"(baryon {m51b_baryon_code:.0f} + halo {m51b_Mhalo_target:.0f})")
print(f"  Mp = {Mp_check:.3f}  (target {m51b_mass_ratio}, paper range 0.5-0.7)")
print(f"  Equivalent as a ratio: 1:{1/Mp_check:.2f}")

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

M_baryon_total = m51a_baryon_code + m51b_baryon_code

def Phi(r):
    phi_halo_a = 0.5 * m51a_haloVc**2 * math.log(r**2 + m51a_haloRc_code**2)
    phi_halo_b = 0.5 * m51b_haloVc**2 * math.log(r**2 + m51b_haloRc_code**2)
    phi_baryon = -M_baryon_total / r
    return phi_halo_a + phi_halo_b + phi_baryon

# === Solve pericenter from the target crossing distance ===
# Rcross is the paper's constraint, so it is the input; pericenter follows.
# The crossing radius is not a closed-form function of pericenter -- it depends
# on the node geometry, and the apsides precess in a logarithmic potential --
# so bisect using a short numerical integration.

def _accel_solve(p):
    r = math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
    g = (M_baryon_total / (r*r)
         + m51a_haloVc**2 * r / (r*r + m51a_haloRc_code**2)
         + m51b_haloVc**2 * r / (r*r + m51b_haloRc_code**2))
    return [-g*p[0]/r, -g*p[1]/r, -g*p[2]/r], r

def _vt_for(r_apo, r_per):
    num = 2.0 * (Phi(r_per) - Phi(r_apo))
    den = 1.0 - (r_apo / r_per)**2
    return math.sqrt(num / den)

def _first_crossing_Rd(peri_kpc, dt=5.0e-5, t_max=12.0):
    """Radius of the first disc-plane crossing, in units of the primary Rd,
    starting from apocenter with argument of pericenter = 90 deg."""
    rp = peri_kpc / du
    ra = rp * (1.0 + orbital_eccentricity) / (1.0 - orbital_eccentricity)
    vt = _vt_for(ra, rp)
    _i = math.radians(orbital_inclination_deg)
    p = [0.0, ra*math.sin(_i), ra*math.cos(_i)]
    v = [-vt, 0.0, 0.0]
    a, _ = _accel_solve(p)
    t = 0.0
    y_prev = p[1]
    while t < t_max:
        for k in range(3):
            p[k] += v[k]*dt + 0.5*a[k]*dt*dt
        a_new, r = _accel_solve(p)
        for k in range(3):
            v[k] += 0.5*(a[k] + a_new[k])*dt
        a = a_new
        t += dt
        if y_prev * p[1] < 0.0:
            return r / m51a_R_code
        y_prev = p[1]
    return None

_lo = 0.5 * target_Rcross * m51a_radius_kpc
_hi = 1.5 * target_Rcross * m51a_radius_kpc
for _ in range(50):
    _mid = 0.5*(_lo + _hi)
    _rr = _first_crossing_Rd(_mid)
    if _rr is None:
        break
    if _rr < target_Rcross:
        _lo = _mid
    else:
        _hi = _mid
pericenter_kpc = 0.5*(_lo + _hi)
apocenter_kpc = pericenter_kpc * (1.0 + orbital_eccentricity) / (1.0 - orbital_eccentricity)

r_peri = pericenter_kpc / du
r_start_kpc = apocenter_kpc
r_start = r_start_kpc / du

print("\n--- Orbital computation ---")
print(f"  Target Rcross:     {target_Rcross:.2f} Rd (paper range 1.2-1.4, using min)")
print(f"  -> Pericenter:     {r_peri:.1f} code units ({pericenter_kpc:.2f} kpc)  [solved]")
print(f"  -> Apocenter:      {r_start:.1f} code units ({apocenter_kpc:.2f} kpc)  [from e]")
print(f"  Eccentricity:      {orbital_eccentricity} (paper bound model)")
print(f"  Check e:           {(r_start-r_peri)/(r_start+r_peri):.4f}")
_S_now = m51b_mass_ratio*(1.0/target_Rcross)**3
_S_paper = m51b_mass_ratio*(1.0/1.20)**3
print(f"  Tidal strength S = Mp*(Rd/r)^3 = {_S_now:.4f}")
print(f"    paper's Rcross=1.20 Rd would give S = {_S_paper:.4f} "
      f"-> this is {_S_now/_S_paper:.1f}x stronger")
print(f"  Crossing is {target_Rcross*4.0:.2f} disc scale lengths out (h_r = Rd/4)")
print(f"  Observed projected separation of NGC 5195 ~11.5 kpc = 0.69 Rd")

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

# Enclosed-mass ratio at pericenter (different quantity: both halos are
# evaluated at the encounter separation rather than at each disc radius)
r_peri_code = pericenter_kpc / du
m51a_halo_peri = m51a_haloVc**2 * r_peri_code**3 / (r_peri_code**2 + m51a_haloRc_code**2)
m51b_halo_peri = m51b_haloVc**2 * r_peri_code**3 / (r_peri_code**2 + m51b_haloRc_code**2)
m51a_dyn_peri = m51a_baryon_code + m51a_halo_peri
m51b_dyn_peri = m51b_baryon_code + m51b_halo_peri
dyn_ratio = m51b_dyn_peri / m51a_dyn_peri

print(f"\n--- Enclosed mass ratio at pericenter ({pericenter_kpc:.2f} kpc) ---")
print(f"  M51a enclosed: {m51a_dyn_peri:.0f} (baryonic {m51a_baryon_code:.0f} + halo {m51a_halo_peri:.0f})")
print(f"  M51b enclosed: {m51b_dyn_peri:.0f} (baryonic {m51b_baryon_code:.0f} + halo {m51b_halo_peri:.0f})")
print(f"  M51b / M51a: 1:{1/dyn_ratio:.2f}")

# === Coordinate system ===
# M51a at origin, disc in the x-z plane, disc normal +y.
#
# NODES: the two points where M51b's orbit crosses M51a's disc plane (y = 0).
# The line joining them (through the origin) is the LINE OF NODES. Because
# gravity here is a central force, angular momentum is conserved and the orbit
# plane -- hence the line of nodes -- is fixed for the whole run.
#
# NODE GEOMETRY is where the pericenter sits relative to that line, i.e. the
# argument of pericenter. It is a genuinely independent degree of freedom: you
# can hold pericenter, apocenter, eccentricity and inclination all fixed and
# still slide the apsides around the orbit, changing where in the orbit the disc
# crossings happen and at what radius.
#
# The paper's bound model (Fig. 1 caption) puts "the apocentre between the two
# disc crossings", which means:
#   - apocenter at MAXIMUM height off the disc plane
#   - line of nodes ~90 deg away in orbital phase (argument of pericenter 90 deg)
#   - crossings therefore occur on the way in/out, NOT at the apsides
#
# Construction: build an orthonormal orbit-plane basis
#   e1 = (1, 0, 0)                        along the line of nodes, in the disc plane
#   e2 = (0, sin(iorb), cos(iorb))        perpendicular to e1, tilted out of the plane
# Placing M51b at apocenter along e2 puts it at maximum |y|, with its velocity
# purely tangential along -e1 (in the disc plane). The nodes are then at +-e1,
# a quarter-orbit away from each apsis -- exactly the paper's arrangement.
#
# Two earlier bugs this replaces:
#   1. pos_x = r_start*cos(inc), which made |pos| smaller than the r_start that
#      v_t was solved for (3% error at 15 deg, 83% at 80 deg).
#   2. apocenter placed at y = 0, i.e. IN the disc plane, putting the nodes near
#      the apsides. M51b then started on the plane, swung to -y, and re-crossed
#      from below -- visibly not "apocentre between the crossings".

inc = math.radians(orbital_inclination_deg)

# Orbit-plane basis (see above)
e1 = (1.0, 0.0, 0.0)                                  # line of nodes
e2 = (0.0, math.sin(inc), math.cos(inc))              # toward apocenter

# M51b at apocenter, i.e. r_start along +e2 -> maximum height off the disc plane
pos_x = r_start * e2[0]
pos_y = r_start * e2[1]
pos_z = r_start * e2[2]

# Velocity purely tangential at apocenter, along +e1 (in the disc plane).
#
# The SIGN sets the circulation sense of the orbit relative to M51a's disc spin,
# and it must be +e1 rather than -e1. The generator gives each disc a spin angular
# momentum along its own normal, so M51a's L_disc points along +y. With velocity
# along +e1 the orbital angular momentum L_orb = pos x vel has a positive
# y-component, making the angle between L_orb and L_disc exactly iorb = 80 deg.
# Using -e1 mirrors it to 180 - iorb = 100 deg: the same near-polar geometry but
# the opposite circulation, i.e. retrograde, and outside the paper's 75-85 range.
# The assertion after this block guards the sign.
vel_x = v_t * e1[0]
vel_y = v_t * e1[1]
vel_z = v_t * e1[2]
# e1 has exact zero components; v_t*0.0 can yield -0.0, which prints as "-0.0".
# Normalize so the emitted script lines are clean.
vel_x += 0.0; vel_y += 0.0; vel_z += 0.0
pos_x += 0.0; pos_y += 0.0; pos_z += 0.0

# === Verify the orbital inclination actually realized ===
# iorb is defined as the angle between the orbital angular momentum and the
# primary's disc angular momentum. Recomputing it from the emitted state vector
# catches a sign or basis error that the inclination input alone cannot.
_L_orb = (pos_y*vel_z - pos_z*vel_y,
          pos_z*vel_x - pos_x*vel_z,
          pos_x*vel_y - pos_y*vel_x)
_L_orb_mag = math.sqrt(sum(c*c for c in _L_orb))
# M51a's disc normal is +y (set in its GalaxyDisc line), and the generator makes
# disc spin angular momentum parallel to the normal.
_L_disc = (0.0, 1.0, 0.0)
_cos_i = sum(a*b for a, b in zip(_L_orb, _L_disc)) / _L_orb_mag
_iorb_actual = math.degrees(math.acos(max(-1.0, min(1.0, _cos_i))))

assert abs(_iorb_actual - orbital_inclination_deg) < 1e-6, (
    f"realized iorb {_iorb_actual:.3f} deg != requested "
    f"{orbital_inclination_deg} deg -- check the sign of the velocity basis "
    f"vector; a flipped sign gives the 180-iorb mirror (retrograde)")

print(f"\n--- Initial conditions (code units) ---")
print(f"  M51b position: ({pos_x:.1f}, {pos_y:.1f}, {pos_z:.1f})")
print(f"  |position|:    {math.sqrt(pos_x**2+pos_y**2+pos_z**2):.1f} "
      f"(must equal apocenter {r_start:.1f})")
print(f"  height off disc plane: {abs(pos_y):.1f} code = {abs(pos_y)*du:.2f} kpc")
print(f"  M51b velocity: ({vel_x:.1f}, {vel_y:.1f}, {vel_z:.1f})")
print(f"  |velocity|:    {math.sqrt(vel_x**2+vel_y**2+vel_z**2):.1f} km/s (= v_t, all tangential)")
print(f"  pos . vel = {pos_x*vel_x+pos_y*vel_y+pos_z*vel_z:.3e} (must be ~0 at apocenter)")
print(f"  Orbit inclination: {_iorb_actual:.1f} deg (requested {orbital_inclination_deg:.0f}),"
      f" verified from L_orb vs L_disc")
print(f"  Line of nodes along +/-x; orbit circulates in the same sense as the")
print(f"    M51a disc (L_orb . L_disc > 0), i.e. prograde, not the 100 deg mirror")
print(f"  Argument of pericenter: 90 deg (apocenter between the disc crossings)")

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
# Salo & Laurikainen section 2.2: adopting PA_disc = 90 deg and i_disc = 30 deg
# for the companion (Schweizer 1977; Smith et al. 1990) "yield a relative
# inclination of 32.5 deg with respect to the disc of M51".
#
# This is a property of the companion's own disc, set by observation, and is
# INDEPENDENT of the orbital inclination. The previous version tied the disc
# normal to the orbit inclination, which conflated two unrelated quantities.
m51b_disc_rel_inclination_deg = 32.5
_dinc = math.radians(m51b_disc_rel_inclination_deg)

m51b_normal_x = math.sin(_dinc)
m51b_normal_y = math.cos(_dinc)
m51b_normal_z = 0.0

print(f"\n--- M51b disc normal ---")
print(f"  ({m51b_normal_x:.4f}, {m51b_normal_y:.4f}, {m51b_normal_z:.4f})")
print(f"  (tilted {m51b_disc_rel_inclination_deg} deg from M51a disc normal;")
print(f"   observed relative disc inclination, independent of the orbit)")

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
print(f"GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   {m51a_bulge_code:.1f} {m51a_Mfrac:.2f} {m51a_R_code:.1f} {m51a_Ri_code:.1f} {m51a_h_r_code:.1f} {toomre_Q:.1f}  {m51a_haloVc:.1f} {m51a_haloRc_code:.1f}")
print(f"")
print(f"# M51b (NGC 5195) at apocenter, near-polar bound orbit (iorb={orbital_inclination_deg:.0f} deg)")
print(f"GalaxyDisc  1   {pos_x:.1f} {pos_y:.1f} {pos_z:.1f}   {vel_x:.1f} {vel_y:.1f} {vel_z:.1f}   {m51b_normal_x:.4f} {m51b_normal_y:.4f} {m51b_normal_z:.4f}   {m51b_bulge_code:.1f} {m51b_Mfrac:.2f} {m51b_R_code:.1f} {m51b_Ri_code:.1f} {m51b_h_r_code:.1f} {toomre_Q:.1f}  {m51b_haloVc:.1f} {m51b_haloRc_code:.1f}")

# === Orbit integration: verify and report the actual event sequence ===
# The analytic v_t solution guarantees the pericenter/apocenter radii, but the
# DISC-PLANE CROSSINGS -- the events the paper actually constrains -- depend on
# the node geometry and can only be found by integrating. This also verifies the
# node geometry is right: with the apocenter between the crossings, the crossing
# radii should sit near 1.3 Rd, not at the apsides.

def _accel(p):
    r = math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
    g = (M_baryon_total / (r*r)
         + m51a_haloVc**2 * r / (r*r + m51a_haloRc_code**2)
         + m51b_haloVc**2 * r / (r*r + m51b_haloRc_code**2))
    return [-g*p[0]/r, -g*p[1]/r, -g*p[2]/r], r

def integrate_events(t_max=24.0, dt=5.0e-5):
    """Velocity-Verlet the relative orbit. Returns a time-ordered list of
    (kind, t, r_code) for disc-plane crossings and apsides."""
    p = [pos_x, pos_y, pos_z]
    v = [vel_x, vel_y, vel_z]
    a, _ = _accel(p)
    t = 0.0
    events = []
    y_prev = p[1]
    r_prev = None
    dr_prev = None
    while t < t_max:
        for k in range(3):
            p[k] += v[k]*dt + 0.5*a[k]*dt*dt
        a_new, r = _accel(p)
        for k in range(3):
            v[k] += 0.5*(a[k] + a_new[k])*dt
        a = a_new
        t += dt
        # disc-plane crossing: y changes sign (M51a disc is the x-z plane)
        if y_prev * p[1] < 0.0:
            events.append(("crossing", t, r))
        y_prev = p[1]
        if r_prev is not None:
            dr = r - r_prev
            if dr_prev is not None:
                if dr_prev < 0.0 and dr >= 0.0:
                    events.append(("pericenter", t, r_prev))
                elif dr_prev > 0.0 and dr <= 0.0:
                    events.append(("apocenter", t, r_prev))
            dr_prev = dr
        r_prev = r
    return events

events = integrate_events()

print("")
print(f"{'='*65}")
print(f"ORBIT EVENTS (numerically integrated)")
print(f"{'='*65}")
print(f"  Disc-plane crossings are what Salo & Laurikainen constrain:")
print(f"    Rcross = 1.2-1.4 Rd (principal), Rdown = 1.2-1.3 Rd (most recent)")
print(f"  Rd = {m51a_R_code:.1f} code units = {m51a_radius_kpc} kpc")
print("")
print(f"  {'event':>12} {'t':>7} {'Myr':>7} {'r (kpc)':>9} {'r / Rd':>8}")
print(f"  {'-'*47}")
n_cross = 0
for kind, te, re_ in events:
    flag = ""
    if kind == "crossing":
        n_cross += 1
        in_range = 1.2 <= re_/m51a_R_code <= 1.4
        flag = "  <- in Rcross range" if in_range else ""
    print(f"  {kind:>12} {te:7.2f} {te*tu:7.0f} {re_*du:9.2f} {re_/m51a_R_code:8.3f}{flag}")

crossings = [(te, re_) for kind, te, re_ in events if kind == "crossing"]

print("")
print(f"  Node geometry check: apocenter should fall BETWEEN crossings")
print(f"  (paper Fig. 1: 'the apocentre is between the two disc crossings')")
for kind, te, re_ in events:
    if kind == "apocenter":
        before = [t for t, _ in crossings if t < te]
        after = [t for t, _ in crossings if t > te]
        if before and after:
            print(f"    apocenter t={te:.2f} lies between crossings "
                  f"t={before[-1]:.2f} and t={after[0]:.2f}  OK")
        break

# Close crossings only. Because the apsides precess in a logarithmic potential,
# the orbit alternates between crossing near the semi-latus rectum (~1.3 Rd, a
# genuine close passage) and crossing near apocenter (~1.8 Rd, a grazing pass
# that does little). Only the close ones are the paper's Rcross / Rdown.
_close_cut = (r_peri + r_start) / 2.0 / m51a_R_code   # semi-major axis, in Rd
close = [(te, re_) for te, re_ in crossings if re_ / m51a_R_code <= _close_cut]

print("")
print(f"  Close crossings (<= {_close_cut:.2f} Rd = semi-major axis; Rcross/Rdown events):")
for te, re_ in close:
    print(f"    t={te:6.2f} ({te*tu:5.0f} Myr)  {re_*du:5.2f} kpc = {re_/m51a_R_code:.2f} Rd")
print(f"  Grazing crossings near apocenter (~1.8 Rd) are skipped: the companion")
print(f"  is at its greatest distance there, so the perturbation is weak.")

if len(close) >= 2:
    t_principal = close[0][0]      # Rcross: the spiral-inducing passage
    t_down = close[1][0]           # Rdown: the most recent crossing
    # Paper: Tobs - Tdown = 0.5-1.0 in ITS time unit (80 Myr) after Rdown
    lo = t_down + 0.5*80.0/tu
    hi = t_down + 1.0*80.0/tu
    print("")
    print(f"  BEST-MORPHOLOGY WINDOW")
    print(f"  The paper's observation epoch (Tobs) is 0.5-1.0 of its own time")
    print(f"  units (80 Myr each) after the most recent close crossing Rdown.")
    print(f"    principal crossing (Rcross): t={t_principal:.2f} at {close[0][1]/m51a_R_code:.2f} Rd")
    print(f"    most recent crossing (Rdown): t={t_down:.2f} at {close[1][1]/m51a_R_code:.2f} Rd")
    print(f"  => observe at t = {lo:.2f} to {hi:.2f}  (centre ~{(lo+hi)/2:.2f})")
    print(f"  Set End_Time above {hi:.1f}.")
    print("")
    dt_principal = (t_down - t_principal)*tu
    print(f"  Tobs cross-check: the paper wants {5.5*80.0:.0f}-{6.5*80.0:.0f} Myr between")
    print(f"  successive close crossings. Here it is {dt_principal:.0f} Myr", end="")
    if 440.0 <= dt_principal <= 520.0:
        print(" -- in range.")
    else:
        print(f" -- {'above' if dt_principal > 520 else 'BELOW'} range.")
        print(f"  The untruncated halo (no Rh cutoff, see docs section 9.3)")
        print(f"  overestimates enclosed mass at large radii, which shortens the")
        print(f"  period for a given orbit shape.")

# === Timeline ===
print("")
print(f"{'='*65}")
print(f"SUMMARY")
print(f"{'='*65}")
print(f"  Pericenter:   {pericenter_kpc:.1f} kpc ({r_peri:.1f} code)")
print(f"  Apocenter:    {apocenter_kpc:.1f} kpc ({r_start:.1f} code), derived from e")
print(f"  Eccentricity: {orbital_eccentricity}")
print(f"  Inclination:  {orbital_inclination_deg:.0f} deg (near-polar)")
print(f"  Arg. of peri: 90 deg (apocenter between crossings)")
print(f"  Mass ratio Mp:{Mp_check:.3f}")
print(f"  Orbital period: {(close[1][0]-close[0][0])*2 if len(close)>1 else float('nan'):.2f} code units"
      if len(close) > 1 else "")
