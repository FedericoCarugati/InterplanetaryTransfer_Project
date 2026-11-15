"""Level 0 - Validation of the Python code against the thesis MATLAB results.

Run from the PythonCode folder:

    python scripts/run_validation.py                 # fast checks (~10 s)
    python scripts/run_validation.py --matlab-grid   # + faithful port of the MATLAB grid (minutes)

Steps
1. Building blocks: Itokawa elements from the Horizons state (car2par) vs MATLAB.
2. Lambert evaluated at the MATLAB optimum (t0, T) vs the MATLAB delta-V.
3. Lambert porkchop on the same 2025-2028 window + SQP refinement, per Earth revolution.
4. (optional) Faithful vectorised port of the three nested MATLAB loops.

Outputs: results/validation/summary.json and figures in results/validation/.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from itokawa_opt.astro.kepler import time_of_flight
from itokawa_opt.constants import AU, DAY
from itokawa_opt.ephemeris.thesis_kepler import thesis_bodies
from itokawa_opt.optimizers.grid_lambert import porkchop
from itokawa_opt.optimizers.grid_matlab import ThesisScenario, run_grid_search
from itokawa_opt.optimizers.local import refine
from itokawa_opt.problems.direct_transfer import DirectTransfer

ROOT = Path(__file__).resolve().parents[1]
REF = json.loads((ROOT / "data/reference/matlab_scenario2_advanced.json").read_text())
OUT = ROOT / "results/validation"
EPOCH0 = np.datetime64("2025-01-04")


def date(t_seconds: float) -> str:
    return str(EPOCH0 + np.timedelta64(int(round(t_seconds)), "s"))[:10]


def rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--matlab-grid",
        action="store_true",
        help="also run the faithful port of the MATLAB grid search",
    )
    parser.add_argument(
        "--step-days", type=float, default=1.0, help="porkchop grid step in days (default 1)"
    )
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    sc = ThesisScenario()
    earth, itokawa = thesis_bodies(sc)
    summary: dict = {}

    # ---------------------------------------------------------------- 1. building blocks
    el = itokawa.elements
    ref_a = REF["asteroid_from_horizons"]
    summary["car2par"] = {
        "a_rel_err": rel(el.a, ref_a["a"]),
        "e_rel_err": rel(el.e, ref_a["e"]),
        "theta0_abs_err": abs(el.theta - ref_a["theta0"]),
    }
    print("1) Itokawa elements from Horizons state")
    print(
        f"   a  = {el.a:.6e} km   (MATLAB {ref_a['a']:.6e}, rel. err "
        f"{summary['car2par']['a_rel_err']:.1e})"
    )
    print(f"   e  = {el.e:.12f}   (MATLAB {ref_a['e']:.12f})")
    print(f"   th0= {el.theta:.10f}     (MATLAB {ref_a['theta0']:.10f})")

    # ---------------------------------------------------------------- 2. Lambert @ MATLAB optimum
    opt = REF["optimum"]
    t_earth = earth.period
    # MATLAB label = anomaly one time step before the state actually used, 3rd revolution
    t_label = float(
        time_of_flight(earth.elements.a, earth.elements.e, 0.0, opt["theta_earth_label"], sc.mu)
    )
    t0_m = t_label + sc.dt + 2 * t_earth
    tof_m = opt["tof"]
    problem = DirectTransfer(
        earth,
        itokawa,
        t0_bounds=(0.0, 3 * t_earth),
        tof_bounds=(30 * DAY, itokawa.period),
        r_min_sun=sc.r_min_sun,
    )
    ev = problem.evaluate(np.array([t0_m]), np.array([tof_m]))
    dv_lam = float(ev["dv"][0])
    summary["lambert_at_matlab_optimum"] = {
        "t0_days": t0_m / DAY,
        "departure": date(t0_m),
        "tof_days": tof_m / DAY,
        "dv": dv_lam,
        "dv1": float(ev["dv1"][0]),
        "dv2": float(ev["dv2"][0]),
        "dv_matlab": opt["dv_total"],
        "rel_err": rel(dv_lam, opt["dv_total"]),
        "transfer_angle_deg": float(np.degrees(ev["angle"][0])),
    }
    print("\n2) Lambert evaluated at the MATLAB optimum (same t0 and T)")
    print(f"   departure {date(t0_m)}, TOF {tof_m / DAY:.2f} d")
    print(
        f"   dV = {dv_lam:.4f} km/s  (MATLAB {opt['dv_total']:.4f}, rel. diff "
        f"{100 * rel(dv_lam, opt['dv_total']):.3f} %)"
    )

    # ---------------------------------------------------------------- 3. Lambert porkchop
    step = args.step_days * DAY
    t0_grid = np.arange(0.0, 3 * t_earth, step)
    tof_grid = np.arange(30 * DAY, itokawa.period, step)
    pc = porkchop(problem, t0_grid, tof_grid)
    print(
        f"\n3) Lambert porkchop: {pc.n_evaluated:,} points in {pc.runtime_s:.2f} s "
        f"(MATLAB grid: ~{REF['runtime_matlab_s'] / 3600:.0f} h)"
    )
    per_rev = []
    for r in range(3):
        mask = (t0_grid >= r * t_earth) & (t0_grid < (r + 1) * t_earth)
        sub = np.where(mask[:, None], pc.dv, np.nan)
        i, j = np.unravel_index(np.nanargmin(sub), sub.shape)
        loc = refine(problem, [t0_grid[i], tof_grid[j]])
        ev_r = problem.evaluate(np.array([loc.x[0]]), np.array([loc.x[1]]))
        per_rev.append(
            {
                "revolution": r + 1,
                "dv_grid": float(sub[i, j]),
                "dv_refined": loc.dv,
                "departure": date(loc.x[0]),
                "t0_days": loc.x[0] / DAY,
                "tof_days": loc.x[1] / DAY,
                "dv1": float(ev_r["dv1"][0]),
                "dv2": float(ev_r["dv2"][0]),
                "r_min_au": float(ev_r["r_min"][0] / AU),
                "transfer_angle_deg": float(np.degrees(ev_r["angle"][0])),
                "sqp_success": loc.success,
                "sqp_nfev": loc.nfev,
                "dv_matlab": REF["dv_min_per_revolution"][r],
            }
        )
    summary["porkchop"] = {
        "n_points": pc.n_evaluated,
        "runtime_s": pc.runtime_s,
        "step_days": args.step_days,
        "per_revolution": per_rev,
    }
    print("   rev | dV MATLAB | dV Lambert (refined) | departure  | TOF [d] | angle [deg]")
    for p in per_rev:
        print(
            f"    {p['revolution']}  |  {p['dv_matlab']:.4f}   |       {p['dv_refined']:.4f}"
            f"         | {p['departure']} | {p['tof_days']:7.2f} "
            f"| {p['transfer_angle_deg']:6.1f}"
        )
    best = min(per_rev, key=lambda p: p["dv_refined"])
    summary["best"] = best

    # ---------------------------------------------------------------- 4. faithful MATLAB port
    if args.matlab_grid:
        print("\n4) Faithful port of the MATLAB grid search (this takes a few minutes)...")
        g = run_grid_search(sc, matlab_compat=True, verbose=True)
        summary["matlab_port"] = {
            "dv_min": g.dv_min,
            "dv_min_per_rev": g.dv_min_per_rev,
            "k_opt_per_rev": g.k_opt_per_rev,
            "k_opt": g.k_opt,
            "s_opt": g.s_opt,
            "j_opt": g.j_opt,
            "tof_days": g.tof / DAY,
            "departure": date(g.t_dep),
            "transfer": g.transfer,
            "n_feasible": g.n_feasible,
            "n_evaluated": g.n_evaluated,
            "runtime_s": g.runtime_s,
        }
        per_rev_txt = [round(x, 4) for x in g.dv_min_per_rev]
        print(
            f"   dV min = {g.dv_min:.4f} km/s (MATLAB {opt['dv_total']:.4f}); "
            f"per rev {per_rev_txt}; {g.n_feasible} feasible conics "
            f"out of {g.n_evaluated:.2e}; {g.runtime_s / 60:.1f} min"
        )

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nSummary written to {OUT / 'summary.json'}")

    if not args.no_plots:
        from itokawa_opt.astro.lambert import lambert
        from itokawa_opt.viz.plots import porkchop_plot, transfer_plot

        porkchop_plot(
            pc,
            markers={
                "MATLAB optimum": (t0_m, tof_m),
                "Lambert + SQP": (best["t0_days"] * DAY, best["tof_days"] * DAY),
            },
            path=str(OUT / "porkchop_2025_2028.png"),
        )
        t0b, tofb = best["t0_days"] * DAY, best["tof_days"] * DAY
        r1, _ = earth.rv(t0b)
        r2, _ = itokawa.rv(t0b + tofb)
        v1, _ = lambert(r1, r2, tofb, sc.mu)
        transfer_plot(earth, itokawa, r1, v1, tofb, path=str(OUT / "transfer_best.png"))
        print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
