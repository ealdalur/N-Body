"""
Analyse the M51 orbit-decay diagnostic.

The simulation, with `OrbitDiagnostic <N>` in the .sim file, writes
`orbit_diagnostic.csv`: for every N steps it logs the two galaxies' barycentre
positions/velocities, their separation, and the radial/tangential split of the
relative velocity.

This script compares the LIVE orbit against the CONSERVATIVE analytic orbit that
`compute_M51.py` integrates. The key quantity is the specific orbital energy

    E_orb = 0.5 * v_rel^2 + Phi(sep)

where Phi is the exact potential of the two truncated cored-isothermal halos plus
the two baryonic point masses (G = 1, code units) -- the same mass model
`compute_M51.py` uses. In the idealised rigid-potential two-body problem E_orb is
conserved. If the live run's E_orb declines steadily, the simulation is losing
orbital energy (spurious drag, e.g. the halo-recentering artifact, or genuine but
under-resolved dynamical friction from the live discs). Jumps localised at
pericentre instead point to integration/softening error at close approach.

Usage:
    python analyze_orbit_diagnostic.py [orbit_diagnostic.csv] [M51.sim]

Defaults: CSV in the current directory; M51.sim next to this script.
A PNG (orbit_diagnostic.png) is written if matplotlib is available.
"""

import csv
import math
import os
import sys

# ---- inputs -----------------------------------------------------------------
here = os.path.dirname(os.path.abspath(__file__))
csv_path = sys.argv[1] if len(sys.argv) > 1 else "orbit_diagnostic.csv"
sim_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "M51.sim")

if not os.path.exists(csv_path):
    sys.exit(f"CSV not found: {csv_path}\n"
             f"Run M51.sim (with 'OrbitDiagnostic 20') first; it writes "
             f"orbit_diagnostic.csv in the simulation's working directory.")


# ---- mass model from the .sim GalaxyDisc lines ------------------------------
def parse_galaxy_discs(path):
    """Return a list of dicts with the halo + baryon parameters of each disc."""
    discs = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s.startswith("GalaxyDisc"):
                continue
            tok = s.split()
            # tok[0]=GalaxyDisc, [1]=sys, [2:5]=pos, [5:8]=vel, [8:11]=normal,
            # [11]=centralMass, [12]=discMass, [13]=R, [14]=Ri, [15]=h_r, [16]=Q,
            # [17]=haloVc, [18]=haloRc, [19]=haloRh, [20]=sigma_z_ratio (opt)
            M_central = float(tok[11]); M_disc = float(tok[12])
            discs.append(dict(
                M_baryon=M_central + M_disc,
                R=float(tok[13]),
                normal=(float(tok[8]), float(tok[9]), float(tok[10])),
                haloVc=float(tok[17]), haloRc=float(tok[18]), haloRh=float(tok[19]),
            ))
    return discs


discs = parse_galaxy_discs(sim_path)
if len(discs) < 2:
    sys.exit(f"Expected >=2 GalaxyDisc lines in {sim_path}, found {len(discs)}")

M_baryon_total = sum(d["M_baryon"] for d in discs[:2])

# Primary disc geometry, for the disc-plane crossings the paper actually
# constrains: Rd sets the units, and the normal defines the plane the companion
# crosses (a crossing is where (companion - primary) . n_primary changes sign).
Rd_primary = discs[0]["R"]
_nn = discs[0]["normal"]
_nmag = math.sqrt(sum(c * c for c in _nn)) or 1.0
n_primary = tuple(c / _nmag for c in _nn)


def M_halo(r, Vc, Rc, Rh):
    """Enclosed halo mass of a cored isothermal sphere, frozen beyond Rh."""
    rr = r if (Rh <= 0.0 or r <= Rh) else Rh
    return Vc * Vc * rr**3 / (rr * rr + Rc * Rc)


def M_enc(r):
    """Total mass governing the relative orbit: both baryons + both halos.
    Mirrors _M_enc() in compute_M51.py and HaloScale() in Simulation.cpp."""
    return (M_baryon_total
            + M_halo(r, discs[0]["haloVc"], discs[0]["haloRc"], discs[0]["haloRh"])
            + M_halo(r, discs[1]["haloVc"], discs[1]["haloRc"], discs[1]["haloRh"]))


# Phi(r) on a grid: Phi(r) = -integral_r^ref M_enc(x)/x^2 dx, Phi(ref)=0.
# Only differences matter, so the reference offset is irrelevant.
_REF = 8000.0
_grid = [1.0 + i * (_REF - 1.0) / 400000 for i in range(400001)]
_phi = [0.0] * len(_grid)
acc = 0.0
for i in range(len(_grid) - 2, -1, -1):  # integrate downward from _REF
    x0, x1 = _grid[i], _grid[i + 1]
    g0 = M_enc(x0) / (x0 * x0)
    g1 = M_enc(x1) / (x1 * x1)
    acc += 0.5 * (g0 + g1) * (x1 - x0)
    _phi[i] = -acc


def Phi(r):
    if r <= _grid[0]:
        return _phi[0]
    if r >= _grid[-1]:
        return 0.0
    # linear interpolation on the uniform grid
    frac = (r - _grid[0]) / (_grid[-1] - _grid[0]) * (len(_grid) - 1)
    j = int(frac)
    if j >= len(_grid) - 1:
        return _phi[-1]
    w = frac - j
    return _phi[j] * (1 - w) + _phi[j + 1] * w


# ---- read the diagnostic CSV ------------------------------------------------
rows = []
with open(csv_path, newline="") as f:
    for row in csv.DictReader(f):
        rows.append({k: float(v) for k, v in row.items()})

# Keep only post-warmup rows. Use a small positive threshold rather than >= 0.0:
# the warmup->active transition logs a row at t = -0.0 (negative zero, which
# passes ">= 0") where the bulk velocity has not yet been applied, so its
# relative velocity is zero. Seeding the analytic orbit from that row would give
# a v=0 radial plunge. The first genuinely post-bulk row is one log-interval later.
post = [r for r in rows if r["t"] > 1e-6]
if len(post) < 2:
    sys.exit("No post-warmup (t > 0) rows in the CSV.")

# The halo centres are the paper's orbital coordinates and, unlike the particle
# barycentres, are not dragged by tidal debris -- so prefer them for the true
# orbit when present. Fall back to the barycentre columns for older CSVs.
have_halo = "hsep" in post[0]

for r in post:
    r["Eorb"] = 0.5 * r["vrel"] ** 2 + Phi(r["sep"])
    if have_halo:
        hvrel = math.sqrt((r["hv1x"]-r["hv0x"])**2 + (r["hv1y"]-r["hv0y"])**2
                          + (r["hv1z"]-r["hv0z"])**2)
        r["hvrel"] = hvrel
        r["Eorb_halo"] = 0.5 * hvrel ** 2 + Phi(r["hsep"])

# Which separation/energy defines "the orbit" for the report.
SEPK = "hsep" if have_halo else "sep"
EK = "Eorb_halo" if have_halo else "Eorb"
label = "halo-centre (core) orbit" if have_halo else "barycentre orbit"


# ---- analytic orbit from the first post-warmup relative state ---------------
r0 = post[0]
if have_halo:
    p = [r0["h1x"] - r0["h0x"], r0["h1y"] - r0["h0y"], r0["h1z"] - r0["h0z"]]
    v = [r0["hv1x"] - r0["hv0x"], r0["hv1y"] - r0["hv0y"], r0["hv1z"] - r0["hv0z"]]
    # radial velocity of the halo-centre orbit, for apocentre detection
    for r in post:
        hdr = [r["h1x"]-r["h0x"], r["h1y"]-r["h0y"], r["h1z"]-r["h0z"]]
        hdv = [r["hv1x"]-r["hv0x"], r["hv1y"]-r["hv0y"], r["hv1z"]-r["hv0z"]]
        r["hvrad"] = sum(hdr[k]*hdv[k] for k in range(3)) / r["hsep"] if r["hsep"] > 0 else 0.0
else:
    p = [r0["c1x"] - r0["c0x"], r0["c1y"] - r0["c0y"], r0["c1z"] - r0["c0z"]]
    v = [r0["v1x"] - r0["v0x"], r0["v1y"] - r0["v0y"], r0["v1z"] - r0["v0z"]]
VRADK = "hvrad" if have_halo else "vrad"


def accel(p):
    r = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    g = M_enc(r) / (r * r)
    return [-g * p[0] / r, -g * p[1] / r, -g * p[2] / r]


dt = 5.0e-4
t_end = post[-1]["t"]
t = r0["t"]
a = accel(p)
ana_t, ana_sep = [t], [math.sqrt(sum(c * c for c in p))]
ana_plane = [sum(p[k] * n_primary[k] for k in range(3))]
while t < t_end:
    for k in range(3):
        p[k] += v[k] * dt + 0.5 * a[k] * dt * dt
    an = accel(p)
    for k in range(3):
        v[k] += 0.5 * (a[k] + an[k]) * dt
    a = an
    t += dt
    ana_t.append(t)
    ana_sep.append(math.sqrt(sum(c * c for c in p)))
    ana_plane.append(sum(p[k] * n_primary[k] for k in range(3)))


def interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    w = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] * (1 - w) + ys[hi] * w


# ---- report -----------------------------------------------------------------
E0 = post[0][EK]
Ef = post[-1][EK]
sep0 = post[0][SEPK]
sep_min_sim = min(r[SEPK] for r in post)
sep_min_ana = min(ana_sep)
du = 0.060  # code length unit in kpc

# E_orb sampled at successive apocentres (vrad + -> -), to see whether the
# conserved quantity steps down passage by passage (real energy loss) or wobbles.
apo_E = []
for i in range(1, len(post)):
    if post[i - 1][VRADK] > 0.0 and post[i][VRADK] <= 0.0:
        apo_E.append((post[i]["t"], post[i][EK], post[i][SEPK]))

print("=" * 68)
print("M51 ORBIT-DECAY DIAGNOSTIC")
print("=" * 68)
print(f"  rows (t>=0):            {len(post)}   t = {post[0]['t']:.2f} .. {t_end:.2f}")
print(f"  orbit measured from:    {label}")
if have_halo:
    bc0, bcf = post[0]["sep"], post[-1]["sep"]
    print(f"  (barycentre separation for reference: {bc0*0.06:.2f} -> {bcf*0.06:.2f} kpc,"
          f" min {min(r['sep'] for r in post)*0.06:.2f} kpc -- includes tidal-debris shift)")
print(f"  M_baryon_total:         {M_baryon_total:.3e} code")
print(f"  E_orb (specific orbital energy, code units):")
print(f"     initial:             {E0:.2f}")
print(f"     final:               {Ef:.2f}")
print(f"     change:              {Ef - E0:+.2f}  ({100.0 * (Ef - E0) / abs(E0):+.1f}%)")
print(f"  Separation:")
print(f"     initial:             {sep0:.1f} code ({sep0 * du:.2f} kpc)")
print(f"     min (sim):           {sep_min_sim:.1f} code ({sep_min_sim * du:.2f} kpc)")
print(f"     min (analytic):      {sep_min_ana:.1f} code ({sep_min_ana * du:.2f} kpc)")
print()
if apo_E:
    print("  E_orb at successive apocentres (steady decline => real energy loss):")
    for tt, ee, ss in apo_E:
        print(f"     t={tt:6.2f}   E_orb={ee:10.2f}   sep={ss:6.1f} ({ss*du:5.1f} kpc)")
    print()

print("  sim vs analytic separation (kpc) at sampled times:")
print(f"  {'t':>7} {'sim':>9} {'analytic':>9} {'sim-ana':>9}")
n = len(post)
for i in range(0, n, max(1, n // 25)):
    r = post[i]
    sa = interp(ana_t, ana_sep, r["t"])
    print(f"  {r['t']:7.2f} {r[SEPK]*du:9.2f} {sa*du:9.2f} {(r[SEPK]-sa)*du:9.2f}")

# ---- disc-plane crossings (the quantity the paper constrains) ---------------
# A crossing is where the companion passes through the primary's disc plane,
# i.e. where (companion - primary) . n_primary changes sign. The paper puts the
# principal crossing (Rcross) at 1.2-1.4 Rd and the most recent (Rdown) at
# 1.2-1.3 Rd. With decay the live crossings tighten over time; this is what item 4
# (orbit re-tuning) checks -- if the later crossings fall below 1.2 Rd, widen
# target_Rcross in compute_M51.py so the decayed orbit lands back in range.
def find_crossings(times, planes, seps):
    out = []
    for i in range(1, len(times)):
        p0, p1 = planes[i - 1], planes[i]
        if p0 == 0.0 or p0 * p1 >= 0.0:
            continue
        f = p0 / (p0 - p1)                       # linear interp to plane = 0
        tc = times[i - 1] + f * (times[i] - times[i - 1])
        rc = seps[i - 1] + f * (seps[i] - seps[i - 1])
        out.append((tc, rc))
    return out

sim_times = [r["t"] for r in post]
sim_seps = [r[SEPK] for r in post]
if have_halo:
    sim_planes = [sum((r["h1" + ax] - r["h0" + ax]) * n_primary[k]
                      for k, ax in enumerate("xyz")) for r in post]
else:
    sim_planes = [sum((r["c1" + ax] - r["c0" + ax]) * n_primary[k]
                      for k, ax in enumerate("xyz")) for r in post]
sim_cross = find_crossings(sim_times, sim_planes, sim_seps)
ana_cross = find_crossings(ana_t, ana_plane, ana_sep)

print()
print("  Disc-plane crossings (paper: Rcross 1.2-1.4 Rd, Rdown 1.2-1.3 Rd):")
if not sim_cross:
    print("    none detected in the run window.")
else:
    print(f"    {'#':>2} {'t':>7} {'sim R/Rd':>9} {'kpc':>7} {'analytic R/Rd':>14}  in-range")
    for i, (tc, rc) in enumerate(sim_cross):
        ac = (min(ana_cross, key=lambda x: abs(x[0] - tc))[1] / Rd_primary
              if ana_cross else float("nan"))
        rrd = rc / Rd_primary
        flag = "yes" if 1.2 <= rrd <= 1.4 else ("LOW" if rrd < 1.2 else "high")
        print(f"    {i+1:>2} {tc:7.2f} {rrd:9.3f} {rc*du:7.2f} {ac:14.3f}  {flag}")
    # Rcross/Rdown are the crossing pair that BRACKETS an apocentre (paper Fig. 1:
    # "the apocentre is between the two disc crossings"), matching compute_M51.py.
    # This is robust to the orbit changing, unlike a hardcoded observation time.
    apo_times = [tt for tt, _, _ in apo_E]
    pair = None
    for i in range(len(sim_cross) - 1):
        ta, tb = sim_cross[i][0], sim_cross[i + 1][0]
        if any(ta < at < tb for at in apo_times):
            pair = (sim_cross[i], sim_cross[i + 1])
            break
    if pair:
        (t_rc, r_rc), (t_rd, r_rd) = pair
        rc_rd, rd_rd = r_rc / Rd_primary, r_rd / Rd_primary
        print(f"    Rcross (principal, before apocentre): {rc_rd:.3f} Rd at t={t_rc:.2f}"
              f"   [paper 1.2-1.4]")
        print(f"    Rdown  (most recent, after apocentre): {rd_rd:.3f} Rd at t={t_rd:.2f}"
              f"   [paper 1.2-1.3]")
        rc_ok = 1.2 <= rc_rd <= 1.4
        rd_ok = 1.2 <= rd_rd <= 1.3
        if rc_ok and rd_ok:
            print("    -> both crossings in the paper's ranges; retune converged.")
        else:
            print(f"    -> scale target_Rcross by ~{1.25/rd_rd:.3f} to centre Rdown at "
                  f"1.25 Rd (Rcross moves with it) and re-run.")
    else:
        print("    (no crossing pair brackets an apocentre in the run window)")

# ---- optional plot ----------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ts = [r["t"] for r in post]
    fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax[0].plot(ts, [r[SEPK] * du for r in post], label=f"sim ({label})")
    ax[0].plot(ana_t, [s * du for s in ana_sep], "--", label="analytic (conservative)")
    if have_halo:
        ax[0].plot(ts, [r["sep"] * du for r in post], ":", color="gray",
                   alpha=0.7, label="sim (barycentre, debris-shifted)")
    ax[0].set_ylabel("separation (kpc)")
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[0].set_title("M51b–M51a separation: live simulation vs conservative analytic orbit")

    ax[1].plot(ts, [r[EK] for r in post], color="C3")
    ax[1].axhline(E0, ls=":", color="gray", label="initial E_orb")
    ax[1].set_ylabel("specific orbital energy  (code)")
    ax[1].set_xlabel("t (code units, 58.7 Myr each)")
    ax[1].legend(); ax[1].grid(alpha=0.3)

    out = os.path.join(os.path.dirname(os.path.abspath(csv_path)), "orbit_diagnostic.png")
    fig.tight_layout(); fig.savefig(out, dpi=110)
    print(f"\n  plot written: {out}")
except ImportError:
    print("\n  (matplotlib not available -- text summary only)")
