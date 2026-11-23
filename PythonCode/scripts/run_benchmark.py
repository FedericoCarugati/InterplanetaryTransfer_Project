"""Level 1 - Run the equal-budget benchmark (DE vs grid vs other baselines).

Run from the PythonCode folder:

    python scripts/run_benchmark.py                    # full benchmark (a few minutes)
    python scripts/run_benchmark.py --quick            # 5 runs per method, smoke test (~30 s)
    python scripts/run_benchmark.py --sensitivity      # + DE parameter study (F x CR)
    python scripts/run_benchmark.py --window-revs 10   # wider launch window (10 Earth years)

Steps
1. Problem: thesis Keplerian model, t0 in [0, N Earth revolutions], T in [30 d, T_Itokawa].
2. Reference and landscape analysis: fine grid + SQP from every grid local minimum gives
   the best-known solution and the list of basins; cost of one evaluation.
3. Deterministic baseline: Lambert grid at several steps, with and without SQP.
4. Stochastic methods, n independent seeded runs each: random search, multi-start SQP,
   DE (three variants), and the hybrid DE + SQP at several budgets.
5. Statistics: success rate (Wilson 95 % CI), median error, ERT, Mann-Whitney U with
   Holm correction and Vargha-Delaney A12 at a common budget.
6. Figures.

Outputs: results/benchmark/summary.json, records_*.json and figures in results/benchmark/.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from itokawa_opt.analysis.benchmark import (
    method_de,
    method_grid,
    method_hybrid,
    method_multistart_sqp,
    method_random_search,
    reference_solution,
    reference_to_dict,
    run_repeated,
    save_records,
)
from itokawa_opt.analysis.stats import (
    at_resolution,
    evaluations_to_target,
    expected_running_time,
    fisher_success,
    holm_correction,
    mann_whitney,
    success_rate,
    summarize,
)
from itokawa_opt.constants import DAY
from itokawa_opt.ephemeris.thesis_kepler import thesis_bodies
from itokawa_opt.optimizers.de import DEConfig, differential_evolution
from itokawa_opt.optimizers.grid_matlab import ThesisScenario
from itokawa_opt.optimizers.pipeline import hybrid_de_sqp, verify_transfer
from itokawa_opt.problems.direct_transfer import DirectTransfer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/benchmark"
EPOCH0 = np.datetime64("2025-01-04")

GRID_STEPS_DAYS = (40.0, 20.0, 10.0, 5.0, 2.0, 1.0)
BENCHMARK_SEED = 2026  # master seed of the benchmark runs
TUNING_SEED = 7  # different master seed for the parameter study (no tuning on test seeds)


def date(t_seconds: float) -> str:
    return str(EPOCH0 + np.timedelta64(int(round(t_seconds)), "s"))[:10]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--runs", type=int, default=30, help="independent runs per method")
    p.add_argument(
        "--budgets",
        type=int,
        nargs="+",
        default=[200, 300, 500, 1000, 2000, 5000],
        help="evaluation budgets at which the methods are compared",
    )
    p.add_argument(
        "--compare-budget", type=int, default=1000, help="budget used for the statistical tests"
    )
    p.add_argument("--tol", type=float, default=1e-3, help="success threshold on dV [km/s]")
    p.add_argument(
        "--resolution",
        type=float,
        default=1e-6,
        help="errors below this are treated as 0 in the rank tests [km/s] (default 1 mm/s)",
    )
    p.add_argument("--pop-size", type=int, default=20, help="DE population size NP")
    p.add_argument("--window-revs", type=int, default=3, help="launch window [Earth years]")
    p.add_argument("--ref-step", type=float, default=0.5, help="reference grid step [days]")
    p.add_argument("--sensitivity", action="store_true", help="also run the F x CR study")
    p.add_argument("--quick", action="store_true", help="5 runs, 3 budgets (smoke test)")
    p.add_argument("--no-plots", action="store_true")
    a = p.parse_args()
    if a.quick:
        a.runs, a.budgets, a.compare_budget = 5, [300, 1000, 2000], 1000
    if a.compare_budget not in a.budgets:
        a.budgets = sorted({*a.budgets, a.compare_budget})
    return a


def main() -> None:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    summary: dict = {"settings": vars(args)}
    b_max = max(args.budgets)

    # ---------------------------------------------------------------- 1. problem
    sc = ThesisScenario()
    earth, itokawa = thesis_bodies(sc)
    problem = DirectTransfer(
        earth,
        itokawa,
        t0_bounds=(0.0, args.window_revs * earth.period),
        tof_bounds=(30 * DAY, itokawa.period),
        r_min_sun=sc.r_min_sun,
    )
    print(
        f"1) Problem: departure {date(0)} .. {date(problem.t0_bounds[1])}, "
        f"TOF 30 .. {itokawa.period / DAY:.0f} days, r_min >= 0.7 AU"
    )

    # ---------------------------------------------------------------- 2. reference
    print(f"\n2) Reference: grid {args.ref_step:g} d + SQP from every grid local minimum")
    ref, pc_ref = reference_solution(problem, step_days=args.ref_step)
    f_ref, tol = ref.f, args.tol
    target = f_ref + tol
    check = verify_transfer(problem, ref.x)
    print(
        f"   {ref.grid_points:,} grid points, {100 * ref.infeasible_fraction:.1f} % infeasible; "
        f"one evaluation costs {ref.eval_time_batch_us:.2f} us in a batch, "
        f"{ref.eval_time_scalar_us:.0f} us alone"
    )
    print("   basin | dV [km/s] | departure  | TOF [d] | on bound")
    for k, b in enumerate(ref.basins, start=1):
        print(f"   {k:5d} | {b.f:9.6f} | {date(b.x[0])} | {b.x[1] / DAY:7.2f} | {b.on_bound}")
    gap = ref.basins[1].f - ref.basins[0].f if len(ref.basins) > 1 else np.inf
    print(
        f"   best-known dV = {f_ref:.9f} km/s; gap to 2nd basin {1e3 * gap:.0f} m/s; "
        f"success threshold {1e3 * tol:g} m/s; independent check: miss "
        f"{check['position_miss_km']:.1e} km"
    )
    if gap <= tol:
        print("   WARNING: threshold larger than the gap between basins - success is ambiguous")
    summary["reference"] = reference_to_dict(ref) | {"verification": check, "gap_km_s": gap}

    # ---------------------------------------------------------------- 3. grid baseline
    print("\n3) Deterministic baseline: Lambert grid (+ SQP from the best node)")
    grid_records = []
    for step in GRID_STEPS_DAYS:
        for polish in (False, True):
            grid_records.append(method_grid(step, polish)(problem, 0, None))
    print("   method              |      nfev | error [m/s] | success")
    for r in grid_records:
        err = 1e3 * (r.f - f_ref)
        print(f"   {r.method:19s} | {r.nfev:9,d} | {err:11.3e} | {err <= 1e3 * tol}")
    save_records(grid_records, OUT / "records_grid.json")

    # ---------------------------------------------------------------- 4. stochastic methods
    np_ = args.pop_size
    anytime_methods = {
        "random search": method_random_search(),
        "multi-start SQP": method_multistart_sqp(),
        "DE/rand/1/bin": method_de(DEConfig(pop_size=np_, strategy="rand1bin"), "DE/rand/1/bin"),
        "DE/best/1/bin": method_de(DEConfig(pop_size=np_, strategy="best1bin"), "DE/best/1/bin"),
        "jDE (rand/1)": method_de(DEConfig(pop_size=np_, self_adaptive=True), "jDE (rand/1)"),
    }
    print(f"\n4) Stochastic methods: {args.runs} runs each, up to {b_max} evaluations")
    anytime = {}
    for name, m in anytime_methods.items():
        anytime[name] = run_repeated(m, problem, b_max, args.runs, BENCHMARK_SEED)
        save_records(anytime[name], OUT / f"records_{_slug(name)}.json")
        print(f"   {name:16s} done ({sum(r.runtime_s for r in anytime[name]):.1f} s)")

    hybrid_name = "DE/rand/1/bin + SQP"
    hybrid_m = method_hybrid(DEConfig(pop_size=np_, strategy="rand1bin"), name=hybrid_name)
    hybrid = {}
    for b in args.budgets:
        hybrid[b] = run_repeated(hybrid_m, problem, b, args.runs, BENCHMARK_SEED)
    save_records([r for b in args.budgets for r in hybrid[b]], OUT / "records_hybrid.json")
    print(f"   {hybrid_name} done at budgets {args.budgets}")

    # ---------------------------------------------------------------- 5. statistics
    print(f"\n5) Statistics (target = best-known + {1e3 * tol:g} m/s)")
    table: dict[str, list] = {}  # name -> [(nfev, rate, lo, hi)]
    errors_at: dict[str, dict[int, np.ndarray]] = {}  # name -> budget -> errors [km/s]
    for name, recs in anytime.items():
        errors_at[name] = {b: np.array([r.best_at(b) - f_ref for r in recs]) for b in args.budgets}
        table[name] = [(b, *success_rate(errors_at[name][b], tol)) for b in args.budgets]
    errors_at[hybrid_name] = {b: np.array([r.f - f_ref for r in hybrid[b]]) for b in args.budgets}
    table[hybrid_name] = [
        (float(np.mean([r.nfev for r in hybrid[b]])), *success_rate(errors_at[hybrid_name][b], tol))
        for b in args.budgets
    ]

    fixed_budget = {}
    print("   success rate [95 % CI] and median error [m/s] at each budget")
    print("   " + "method".ljust(22) + "".join(f"{b:>26d}" for b in args.budgets))
    for name in errors_at:
        cells = []
        fixed_budget[name] = {}
        for b, row in zip(args.budgets, table[name], strict=True):
            med = 1e3 * float(np.median(errors_at[name][b]))
            cells.append(f"{row[1]:4.2f} [{row[2]:.2f},{row[3]:.2f}] {max(med, 0.0):8.1e}")
            fixed_budget[name][b] = {
                "nfev_mean": row[0],
                "success_rate": row[1],
                "ci95": [row[2], row[3]],
                "error_m_s": summarize(1e3 * errors_at[name][b]).__dict__,
            }
        print("   " + name.ljust(22) + "".join(f"{c:>26s}" for c in cells))
    summary["fixed_budget"] = fixed_budget

    # Fixed target: evaluations-to-target and ERT (anytime methods)
    evals = {
        n: np.array([evaluations_to_target(r.history_nfev, r.history_f, target) for r in recs])
        for n, recs in anytime.items()
    }
    summary["fixed_target"] = {}
    print("\n   evaluations to reach the target (successful runs) and ERT")
    for n, e in evals.items():
        spent = [r.nfev for r in anytime[n]]
        ert = expected_running_time(e, spent)
        ok = e[np.isfinite(e)]
        med = float(np.median(ok)) if ok.size else np.inf
        summary["fixed_target"][n] = {
            "successes": int(ok.size),
            "median_evals": med,
            "ert": ert,
        }
        print(f"   {n:16s}: {ok.size:3d}/{e.size} runs, median {med:8.0f}, ERT {ert:9.0f}")
    grid_ok = [r for r in grid_records if r.f - f_ref <= tol]
    cheapest = min(grid_ok, key=lambda r: r.nfev) if grid_ok else None
    if cheapest is not None:
        print(f"   cheapest successful grid: {cheapest.method} with {cheapest.nfev:,} evaluations")
        summary["fixed_target"]["cheapest grid"] = {
            "method": cheapest.method,
            "nfev": cheapest.nfev,
        }

    # Hypothesis tests at the common budget
    bc = args.compare_budget
    res_ = args.resolution
    print(
        f"\n   Tests at {bc} evaluations, {hybrid_name} vs each (Holm-corrected p-values)\n"
        f"   - Mann-Whitney U on the errors (differences below {1e3 * res_:g} m/s ignored)\n"
        f"   - Fisher exact test on the success counts"
    )
    others = list(anytime)
    comps = [
        mann_whitney(
            at_resolution(errors_at[hybrid_name][bc], res_),
            at_resolution(errors_at[n][bc], res_),
            hybrid_name,
            n,
        )
        for n in others
    ]
    fishers = [
        fisher_success(errors_at[hybrid_name][bc], errors_at[n][bc], tol, hybrid_name, n)
        for n in others
    ]
    p_mw = holm_correction([c.p_value for c in comps])
    p_fi = holm_correction([c.p_value for c in fishers])
    summary["tests"] = []
    for c, fi, pm, pf in zip(comps, fishers, p_mw, p_fi, strict=True):
        verdict = c.verdict if pm < 0.05 else "no significant difference"
        print(
            f"   vs {c.name_b:16s}: A12 = {c.a12:.2f}, p_MW = {pm:.1e} ({verdict}); "
            f"successes {fi.successes_a}/{fi.n_a} vs {fi.successes_b}/{fi.n_b}, p_Fisher = {pf:.1e}"
        )
        summary["tests"].append(
            {"mann_whitney": c.__dict__ | {"p_holm": pm}, "fisher": fi.__dict__ | {"p_holm": pf}}
        )

    # Precision of the local phase
    h = hybrid[bc]
    de_only = 1e3 * (np.array([r.extra["dv_de_only"] for r in h]) - f_ref)
    after = 1e3 * (np.array([r.f for r in h]) - f_ref)
    nloc = [r.extra["nfev_local"] for r in h]
    print(
        f"\n   effect of SQP at {bc} evaluations: median error {np.median(de_only):.1e} -> "
        f"{np.median(np.clip(after, 0, None)):.1e} m/s for {np.mean(nloc):.0f} extra evaluations"
    )
    summary["local_phase"] = {
        "median_error_before_m_s": float(np.median(de_only)),
        "median_error_after_m_s": float(np.median(after)),
        "mean_local_nfev": float(np.mean(nloc)),
    }

    # ---------------------------------------------------------------- optional: F x CR
    if args.sensitivity:
        summary["sensitivity"] = sensitivity_study(problem, target, args)

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"\nSummary written to {OUT / 'summary.json'}")

    # ---------------------------------------------------------------- 6. figures
    if not args.no_plots:
        make_figures(
            problem,
            pc_ref,
            ref,
            grid_records,
            anytime,
            hybrid,
            hybrid_name,
            table,
            errors_at,
            evals,
            args,
            summary,
        )
        print(f"Figures written to {OUT}")


def sensitivity_study(problem, target, args, budget: int = 500, runs: int = 20) -> dict:
    """Success rate of DE/rand/1/bin (constant F, no dither) on an F x CR grid.

    Uses a different master seed from the benchmark: parameters are not tuned on the
    same runs used to compare the methods.
    """
    F_values, CR_values = [0.3, 0.5, 0.7, 0.9], [0.1, 0.3, 0.5, 0.7, 0.9]
    rates = np.zeros((len(F_values), len(CR_values)))
    print(f"\n   Sensitivity study: {runs} runs per (F, CR) at {budget} evaluations")
    for i, F in enumerate(F_values):
        for j, CR in enumerate(CR_values):
            cfg = DEConfig(
                pop_size=args.pop_size,
                F=F,
                CR=CR,
                max_nfev=budget,
                stall_generations=10**9,
                xtol=0.0,
            )
            seeds = np.random.SeedSequence(TUNING_SEED).spawn(runs)
            f = [
                differential_evolution(problem.fitness_batch, problem.get_bounds(), cfg, s).f
                for s in seeds
            ]
            rates[i, j] = np.mean(np.array(f) <= target)
        print("   F = " + f"{F:.1f}: " + "  ".join(f"{r:.2f}" for r in rates[i]))
    return {"F": F_values, "CR": CR_values, "success_rate": rates.tolist(), "budget": budget}


def make_figures(
    problem,
    pc_ref,
    ref,
    grid_records,
    anytime,
    hybrid,
    hybrid_name,
    table,
    errors_at,
    evals,
    args,
    summary,
) -> None:
    from itokawa_opt.viz import benchmark_plots as bp

    f_ref, tol = ref.f, args.tol
    bp.landscape_plot(pc_ref, ref, OUT / "landscape.png")

    # One DE run with population snapshots (exploration -> exploitation)
    demo = hybrid_de_sqp(
        problem,
        DEConfig(pop_size=args.pop_size, max_nfev=1000, stall_generations=10**9, xtol=0.0),
        seed=np.random.SeedSequence(BENCHMARK_SEED).spawn(1)[0],
        snapshot_generations=(0, 5, 15, 45),
    )
    bp.population_plot(pc_ref, demo.de.snapshots, OUT / "de_population.png")

    grids = [r for r in grid_records if "SQP" not in r.method]
    grids_sqp = [r for r in grid_records if "SQP" in r.method]
    fixed = {
        "grid": [(r.nfev, r.f) for r in grids],
        "grid + SQP": [(r.nfev, r.f) for r in grids_sqp],
        hybrid_name: [
            (np.mean([r.nfev for r in hybrid[b]]), [r.f for r in hybrid[b]]) for b in args.budgets
        ],
    }
    bp.convergence_plot(anytime, fixed, f_ref, tol, OUT / "convergence.png")

    succ = dict(table)
    succ["grid + SQP"] = [
        (r.nfev, float(r.f - f_ref <= tol), *(float(r.f - f_ref <= tol),) * 2) for r in grids_sqp
    ]
    succ["grid"] = [
        (r.nfev, float(r.f - f_ref <= tol), *(float(r.f - f_ref <= tol),) * 2) for r in grids
    ]
    bp.success_plot(succ, OUT / "success_rate.png")

    det = {r.method: r.nfev for r in grid_records if r.f - f_ref <= tol}
    det = dict(sorted(det.items(), key=lambda kv: kv[1])[:2])  # the two cheapest
    bp.ecdf_plot(evals, max(args.budgets), det, OUT / "ecdf_evaluations.png")

    bc = args.compare_budget
    bp.final_error_boxplot(
        {n: errors_at[n][bc] for n in errors_at}, tol, bc, OUT / "final_errors.png"
    )
    ends = {n: anytime[n] for n in ("random search", "multi-start SQP", "DE/best/1/bin")}
    ends[hybrid_name] = hybrid[min(args.budgets)]
    bp.final_points_plot(pc_ref, ends, OUT / "end_points.png")

    if "sensitivity" in summary:
        s = summary["sensitivity"]
        bp.sensitivity_heatmap(
            s["F"], s["CR"], np.array(s["success_rate"]), s["budget"], OUT / "sensitivity_F_CR.png"
        )


def _slug(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name.lower()).strip("_")


if __name__ == "__main__":
    main()