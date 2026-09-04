"""
TEMPORARY analyzer for the bar/arm diagnostic (bar_diagnostic.csv).

The simulation, with `BarDiagnostic <N>` set, logs per system and per cylindrical
radial bin the m = 1..6 Fourier amplitudes and phases of the disc's mass
distribution:

    t, system, bin, r_lo, r_hi, n, A1..A6, ph1..ph6

  A_m   = |sum_j m_j e^{i m phi_j}| / sum_j m_j   in [0,1]
  ph_m  = atan2(sum m sin(m phi), sum m cos(m phi)) / m   -- order-m orientation

Interpretation:
  - smooth disc:   A_m ~ 1/sqrt(N_bin) shot-noise floor (small)
  - BAR:           strong A2 in the INNER disc, radially COHERENT (near-constant
                   ph2) -- a bar is a single fixed-orientation m=2 mode
  - ARMS:          in the OUTER disc, whichever A_m dominates gives the arm number
                   (m=2 two-arm, m=3 three-arm, m=4 four-arm, ...); a phase ph_m
                   that WINDS with radius marks a trailing/leading spiral

Reports, per system:
  - bar strength A2_bar(t) in the inner disc (peak A2 + phase coherence + verdict);
  - the dominant outer-disc arm order m and its amplitude (the "number of arms"),
    with whether its phase winds (spiral) vs stays constant (bar/oval).
If matplotlib is present, plots A2_bar(t), the dominant arm-m over time, and the
A_m(R) profiles (all m) at the last logged time.

Usage:
    python analyze_bar_diagnostic.py [bar_diagnostic.csv]

Default CSV: build/Release/bar_diagnostic.csv relative to the repo root, then cwd.
Radii are code units (x0.06 -> kpc). Slated for removal with the diagnostic.
"""

import csv
import math
import os
import sys
from collections import defaultdict

KPC_PER_CODE = 0.06
MMAX = 6
R_BAR_MAX_KPC = 7.0          # inner region searched for a bar (MW bar ~3-5 kpc)
R_ARM_MIN_KPC = 8.0          # outer "arm region" starts beyond a typical bar (~5 kpc)
                             # so the bar's own m=2 doesn't masquerade as a 2-arm count
R_ARM_MAX_KPC = 18.0
A2_BAR_STRONG = 0.20
A2_BAR_WEAK = 0.10
A_ARM_MIN = 0.05             # min amplitude to call an outer m an "arm"
PHASE_COH_DEG = 20.0         # phase spread below this = coherent (bar/oval, not winding)


def find_csv(argv):
    if len(argv) > 1:
        return argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    for cand in (os.path.join(root, "build", "Release", "bar_diagnostic.csv"),
                 "bar_diagnostic.csv"):
        if os.path.exists(cand):
            return cand
    return os.path.join(root, "build", "Release", "bar_diagnostic.csv")


def load(path):
    by_sys_t = defaultdict(list)
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                rec = {
                    "t": float(r["t"]), "r_lo": float(r["r_lo"]),
                    "r_hi": float(r["r_hi"]), "n": float(r["n"]),
                    "A": {m: float(r[f"A{m}"]) for m in range(1, MMAX + 1)},
                    "ph": {m: float(r[f"ph{m}"]) for m in range(1, MMAX + 1)},
                }
            except (ValueError, TypeError, KeyError):
                continue
            by_sys_t[(int(r["system"]), rec["t"])].append(rec)
    for v in by_sys_t.values():
        v.sort(key=lambda d: d["r_lo"])
    return by_sys_t


def rmid_kpc(b):
    return 0.5 * (b["r_lo"] + b["r_hi"]) * KPC_PER_CODE


def circ_spread_deg(angles, m):
    """Circular std (deg) of order-m orientation angles (mod pi/... ), via m*angle."""
    if len(angles) < 2:
        return 0.0
    cx = sum(math.cos(m * a) for a in angles) / len(angles)
    cy = sum(math.sin(m * a) for a in angles) / len(angles)
    R = math.hypot(cx, cy)
    if R <= 1e-12:
        return 90.0
    return math.degrees(math.sqrt(max(0.0, -2.0 * math.log(R))) / m)


def bar_metrics(bins):
    inner = [b for b in bins if rmid_kpc(b) <= R_BAR_MAX_KPC and b["n"] > 0]
    if not inner:
        return 0.0, 0.0, 90.0
    peak = max(inner, key=lambda b: b["A"][2])
    a2 = peak["A"][2]
    sig = [b["ph"][2] for b in inner if b["A"][2] >= 0.5 * a2]
    return a2, rmid_kpc(peak), circ_spread_deg(sig, 2)


def bar_verdict(a2, coh):
    if a2 >= A2_BAR_STRONG and coh <= PHASE_COH_DEG:
        return "BAR"
    if a2 >= A2_BAR_WEAK:
        return "weak bar/oval" if coh <= PHASE_COH_DEG else "transient m=2"
    return "no bar"


def arm_metrics(bins):
    """Dominant outer arm order: (m, peak A_m, winds?) over the arm region."""
    arm = [b for b in bins if R_ARM_MIN_KPC <= rmid_kpc(b) <= R_ARM_MAX_KPC and b["n"] > 0]
    if not arm:
        return 0, 0.0, False
    best_m, best_A = 0, 0.0
    for m in range(2, MMAX + 1):
        peakA = max(b["A"][m] for b in arm)
        if peakA > best_A:
            best_A, best_m = peakA, m
    if best_A < A_ARM_MIN or best_m == 0:
        return 0, best_A, False
    sig = [b["ph"][best_m] for b in arm if b["A"][best_m] >= 0.5 * best_A]
    winds = circ_spread_deg(sig, best_m) > PHASE_COH_DEG
    return best_m, best_A, winds


def summarise(by_sys_t):
    systems = sorted({s for (s, _) in by_sys_t})
    for sysid in systems:
        times = sorted(t for (s, t) in by_sys_t if s == sysid)
        print(f"\n=== System {sysid}: bar + arm structure vs time ===")
        print(f"{'t':>8}{'A2_bar':>9}{'bar':>16}{'arm_m':>7}{'arm_A':>8}  arm type")
        step = max(1, len(times) // 25)
        for i, t in enumerate(times):
            if i % step and t != times[-1]:
                continue
            a2, rp, coh = bar_metrics(by_sys_t[(sysid, t)])
            m, aA, winds = arm_metrics(by_sys_t[(sysid, t)])
            arm = "-" if m == 0 else f"{m}-arm {'spiral' if winds else '(coherent)'}"
            print(f"{t:>8.2f}{a2:>9.3f}{bar_verdict(a2, coh):>16}{m:>7}{aA:>8.3f}  {arm}")
        # Late-time conclusion: dominant arm m averaged over the last few snapshots.
        late = times[max(0, len(times) - 5):]
        counts = defaultdict(float)
        for t in late:
            m, aA, _ = arm_metrics(by_sys_t[(sysid, t)])
            if m:
                counts[m] += aA
        dom = max(counts, key=counts.get) if counts else 0
        a2, rp, coh = bar_metrics(by_sys_t[(sysid, times[-1])])
        print(f"  -> late-time: bar A2={a2:.3f} ({bar_verdict(a2, coh)}); "
              f"dominant arm order m={dom if dom else 'none'} "
              f"({'that many arms' if dom else 'no coherent arms'})")


def plot(by_sys_t):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib unavailable: {e}; skipping plot)")
        return
    systems = sorted({s for (s, _) in by_sys_t})
    fig, axes = plt.subplots(len(systems), 2, figsize=(13, 4.0 * len(systems)),
                             squeeze=False)
    for row, sysid in enumerate(systems):
        times = sorted(t for (s, t) in by_sys_t if s == sysid)
        axL = axes[row][0]
        axL.plot(times, [bar_metrics(by_sys_t[(sysid, t)])[0] for t in times],
                 label="A2_bar (inner)")
        axL.plot(times, [arm_metrics(by_sys_t[(sysid, t)])[1] for t in times],
                 label="dominant arm A_m (outer)", ls="--")
        axL.axhline(A2_BAR_STRONG, color="r", ls=":", lw=0.8)
        axL.set_title(f"System {sysid}: bar & arm amplitude vs time")
        axL.set_xlabel("t [code units]"); axL.set_ylabel("Fourier amplitude")
        axL.legend(fontsize=8)
        axR = axes[row][1]
        last = by_sys_t[(sysid, times[-1])]
        rr = [rmid_kpc(b) for b in last]
        for m in range(1, MMAX + 1):
            axR.plot(rr, [b["A"][m] for b in last], marker=".", ms=2, label=f"m={m}")
        axR.axvspan(0, R_BAR_MAX_KPC, color="grey", alpha=0.08)
        axR.set_title(f"System {sysid}: A_m(R) at t={times[-1]:.2f} (grey = bar region)")
        axR.set_xlabel("R [kpc]"); axR.set_ylabel("A_m")
        axR.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    out = "bar_diagnostic.png"
    fig.savefig(out, dpi=110)
    print(f"\nWrote {out}")


def main():
    path = find_csv(sys.argv)
    if not os.path.exists(path):
        print(f"CSV not found: {path}\nRun a sim with `BarDiagnostic <N>` first.")
        sys.exit(1)
    print(f"Reading {path}")
    by_sys_t = load(path)
    if not by_sys_t:
        print("No rows logged (or old single-m CSV format -- re-run to regenerate).")
        sys.exit(1)
    summarise(by_sys_t)
    plot(by_sys_t)


if __name__ == "__main__":
    main()
