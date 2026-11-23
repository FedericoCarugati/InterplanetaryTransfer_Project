"""Level 1 - Equal-budget benchmark of optimisers over many independent seeded runs.

Every method is wrapped in a function with the same signature

    method(problem, budget, seed) -> RunRecord

where ``budget`` is the maximum number of objective evaluations (Lambert solutions) and
``seed`` a ``numpy.random.SeedSequence``. The *number of evaluations* is the currency of
the comparison: it does not depend on the computer, unlike the wall-clock time (which is
recorded too). Deterministic methods (grid search) ignore the seed and are run once.

A ``RunRecord`` stores the final result and the *anytime* convergence history
(best feasible value found after n evaluations), from which the statistics of
``analysis.stats`` are computed at any budget.

Independent random streams for the runs come from ``SeedSequence(master).spawn(n)``:
consecutive integer seeds (0, 1, 2, ...) are fine for NumPy's generator too, but spawned
sequences are guaranteed independent by construction and are the documented best practice.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import minimum_filter

from itokawa_opt.constants import DAY
from itokawa_opt.optimizers.de import DEConfig, differential_evolution
from itokawa_opt.optimizers.grid_lambert import PorkchopResult, porkchop
from itokawa_opt.optimizers.local import refine
from itokawa_opt.optimizers.pipeline import hybrid_de_sqp
from itokawa_opt.problems.direct_transfer import DirectTransfer

FEAS_TOL = 1e-9  # constraint tolerance used to call a point feasible


# --------------------------------------------------------------------------------------
# Data containers
# --------------------------------------------------------------------------------------
@dataclass
class RunRecord:
    method: str
    seed: int | None  # entropy of the SeedSequence (None for deterministic methods)
    x: list[float]  # final [t0, T] [s]
    f: float  # final total delta-V [km/s] (inf if infeasible)
    feasible: bool
    nfev: int
    runtime_s: float
    history_nfev: list[int] = field(default_factory=list)
    history_f: list[float] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    def best_at(self, budget: float) -> float:
        """Best feasible value found within ``budget`` evaluations (step function)."""
        h_n = np.asarray(self.history_nfev)
        h_f = np.asarray(self.history_f)
        ok = h_n <= budget
        return float(h_f[ok].min()) if ok.any() else np.inf


@dataclass
class Basin:
    """A local minimum of the landscape (found by SQP from a grid local minimum)."""

    x: list[float]
    f: float
    on_bound: bool  # True if the minimum lies on the box (constraint active)


@dataclass
class Reference:
    """Best-known solution and landscape description, computed independently."""

    x: list[float]
    f: float
    basins: list[Basin]
    grid_step_days: float
    grid_points: int
    infeasible_fraction: float
    eval_time_batch_us: float  # cost of one evaluation in a population call
    eval_time_scalar_us: float  # cost of one evaluation called one at a time


# --------------------------------------------------------------------------------------
# Reference solution and landscape analysis
# --------------------------------------------------------------------------------------
def grid_axes(problem: DirectTransfer, step_days: float) -> tuple[NDArray, NDArray]:
    lb, ub = problem.get_bounds()
    step = step_days * DAY
    return np.arange(lb[0], ub[0], step), np.arange(lb[1], ub[1], step)


def find_basins(
    problem: DirectTransfer, pc: PorkchopResult, n_max: int = 30, window: int = 7
) -> list[Basin]:
    """Local minima of the grid, refined with SQP and de-duplicated.

    A grid node is a candidate if it is the minimum of its ``window x window``
    neighbourhood. SQP from each candidate converges to the bottom of its basin; minima
    closer than 1e-3 (normalised box) are the same basin.
    """
    dv = np.where(np.isfinite(pc.dv), pc.dv, np.inf)
    is_min = (minimum_filter(dv, size=window, mode="nearest") == dv) & np.isfinite(dv)
    ii, jj = np.nonzero(is_min)
    order = np.argsort(dv[ii, jj])[:n_max]
    lb, ub = (np.asarray(b, dtype=float) for b in problem.get_bounds())
    span = ub - lb
    basins: list[Basin] = []
    for k in order:
        loc = refine(problem, [pc.t0[ii[k]], pc.tof[jj[k]]])
        if loc.constraint > FEAS_TOL:
            continue
        z = (loc.x - lb) / span
        if any(np.linalg.norm(z - (np.asarray(b.x) - lb) / span) < 1e-3 for b in basins):
            continue
        on_bound = bool(np.any((z < 1e-6) | (z > 1 - 1e-6)))
        basins.append(Basin(x=loc.x.tolist(), f=float(loc.dv), on_bound=on_bound))
    return sorted(basins, key=lambda b: b.f)


def time_evaluations(problem: DirectTransfer, n_batch: int = 20_000, n_scalar: int = 300):
    """Average cost of one evaluation [us], in a batch call and in single calls."""
    lb, ub = (np.asarray(b, dtype=float) for b in problem.get_bounds())
    x = lb + np.random.default_rng(0).random((n_batch, 2)) * (ub - lb)
    problem.fitness_batch(x[:10])  # make sure the Numba kernels are compiled
    t = time.perf_counter()
    problem.fitness_batch(x)
    batch = (time.perf_counter() - t) / n_batch * 1e6
    t = time.perf_counter()
    for row in x[:n_scalar]:
        problem.fitness(row)
    scalar = (time.perf_counter() - t) / n_scalar * 1e6
    return batch, scalar


def reference_solution(
    problem: DirectTransfer, step_days: float = 0.5
) -> tuple[Reference, PorkchopResult]:
    """Fine grid + SQP from every grid local minimum: the "ground truth" of the benchmark.

    It is computed with a *different* procedure from the methods under test, and much
    more expensively, so that it is not biased towards any of them.
    """
    t0, tof = grid_axes(problem, step_days)
    pc = porkchop(problem, t0, tof)
    basins = find_basins(problem, pc)
    batch, scalar = time_evaluations(problem)
    best = basins[0]
    ref = Reference(
        x=best.x,
        f=best.f,
        basins=basins,
        grid_step_days=step_days,
        grid_points=pc.n_evaluated,
        infeasible_fraction=float(np.mean(~np.isfinite(pc.dv))),
        eval_time_batch_us=batch,
        eval_time_scalar_us=scalar,
    )
    return ref, pc


# --------------------------------------------------------------------------------------
# Methods under test (same signature)
# --------------------------------------------------------------------------------------
Method = Callable[[DirectTransfer, int, np.random.SeedSequence | None], RunRecord]


def method_grid(step_days: float, polish: bool = False) -> Method:
    """Thesis-like baseline: porkchop grid with Lambert (+ optional SQP from the best node).

    Deterministic; ``budget`` is ignored (the cost is fixed by the step).
    """
    name = f"grid {step_days:g} d" + (" + SQP" if polish else "")

    def run(problem, budget, seed=None):
        t_start = time.perf_counter()
        pc = porkchop(problem, *grid_axes(problem, step_days))
        x, f, nfev = [pc.best_t0, pc.best_tof], pc.best_dv, pc.n_evaluated
        hist_n, hist_f = [nfev], [f]
        if polish:
            loc = refine(problem, x)
            nfev += loc.nfev
            if loc.constraint <= FEAS_TOL and loc.dv < f:
                x, f = loc.x.tolist(), loc.dv
            hist_n.append(nfev)
            hist_f.append(f)
        return RunRecord(
            method=name,
            seed=None,
            x=[float(v) for v in x],
            f=float(f),
            feasible=True,
            nfev=int(nfev),
            runtime_s=time.perf_counter() - t_start,
            history_nfev=hist_n,
            history_f=[float(v) for v in hist_f],
            extra={"step_days": step_days},
        )

    run.__name__ = name
    return run


def method_random_search(batch: int = 100) -> Method:
    """Uniform random sampling of the box: the sanity floor every optimiser must beat."""

    def run(problem, budget, seed=None):
        t_start = time.perf_counter()
        rng = np.random.default_rng(seed)
        lb, ub = (np.asarray(b, dtype=float) for b in problem.get_bounds())
        best_f, best_x, nfev = np.inf, lb.copy(), 0
        hist_n, hist_f = [], []
        while nfev < budget:
            m = min(batch, budget - nfev)
            x = lb + rng.random((m, lb.size)) * (ub - lb)
            f, g = problem.fitness_batch(x)
            f = np.where(g <= FEAS_TOL, f, np.inf)
            nfev += m
            k = int(np.argmin(f))
            if f[k] < best_f:
                best_f, best_x = float(f[k]), x[k]
            hist_n.append(nfev)
            hist_f.append(best_f)
        return _record("random search", seed, best_x, best_f, nfev, t_start, hist_n, hist_f)

    return run


def method_multistart_sqp(maxiter: int = 100) -> Method:
    """SQP from uniformly random points until the budget is spent (local-only baseline)."""

    def run(problem, budget, seed=None):
        t_start = time.perf_counter()
        rng = np.random.default_rng(seed)
        lb, ub = (np.asarray(b, dtype=float) for b in problem.get_bounds())
        best_f, best_x, nfev, starts = np.inf, lb.copy(), 0, 0
        hist_n, hist_f = [], []
        while nfev < budget:
            loc = refine(problem, lb + rng.random(lb.size) * (ub - lb), maxiter=maxiter)
            nfev += loc.nfev
            starts += 1
            if loc.constraint <= FEAS_TOL and loc.dv < best_f:
                best_f, best_x = float(loc.dv), loc.x
            hist_n.append(nfev)
            hist_f.append(best_f)
        rec = _record("multi-start SQP", seed, best_x, best_f, nfev, t_start, hist_n, hist_f)
        rec.extra["starts"] = starts
        return rec

    return run


def method_de(config: DEConfig, name: str | None = None) -> Method:
    """Differential Evolution alone, with the whole budget (anytime history)."""
    label = name or _de_label(config)

    def run(problem, budget, seed=None):
        cfg = _with_budget(config, budget)
        res = differential_evolution(problem.fitness_batch, problem.get_bounds(), cfg, seed=seed)
        f = res.f if res.feasible else np.inf
        rec = RunRecord(
            method=label,
            seed=_entropy(seed),
            x=res.x.tolist(),
            f=float(f),
            feasible=res.feasible,
            nfev=res.nfev,
            runtime_s=res.runtime_s,
            history_nfev=res.history_nfev.tolist(),
            history_f=res.history_f.tolist(),
            extra={"ngen": res.ngen, "message": res.message},
        )
        return rec

    return run


def method_hybrid(
    config: DEConfig, local_reserve: int = 60, n_candidates: int = 3, name: str | None = None
) -> Method:
    """DE + SQP. DE gets ``budget - local_reserve`` evaluations, SQP the rest (about).

    The local phase is not truncated: its true cost is added to ``nfev``, so the total
    may slightly exceed ``budget``. The statistics use the actual ``nfev``.
    """
    label = name or (_de_label(config) + " + SQP")

    def run(problem, budget, seed=None):
        de_budget = max(config.pop_size, budget - local_reserve)
        res = hybrid_de_sqp(
            problem, _with_budget(config, de_budget), seed=seed, n_candidates=n_candidates
        )
        f = res.dv if res.feasible else np.inf
        return RunRecord(
            method=label,
            seed=_entropy(seed),
            x=res.x.tolist(),
            f=float(f),
            feasible=res.feasible,
            nfev=res.nfev,
            runtime_s=res.runtime_s,
            history_nfev=res.history_nfev.tolist(),
            history_f=res.history_f.tolist(),
            extra={
                "nfev_global": res.nfev_global,
                "nfev_local": res.nfev_local,
                "dv_de_only": res.de.f,
                "budget": budget,
            },
        )

    return run


# --------------------------------------------------------------------------------------
# Repeated runs and persistence
# --------------------------------------------------------------------------------------
def run_repeated(
    method: Method,
    problem: DirectTransfer,
    budget: int,
    n_runs: int,
    master_seed: int = 2026,
    verbose: bool = False,
) -> list[RunRecord]:
    """``n_runs`` independent runs, seeds spawned from ``master_seed``."""
    seeds = np.random.SeedSequence(master_seed).spawn(n_runs)
    records = []
    for k, s in enumerate(seeds):
        records.append(method(problem, budget, s))
        if verbose:
            r = records[-1]
            print(f"    run {k + 1:3d}/{n_runs}: dV = {r.f:.6f} km/s, nfev = {r.nfev}")
    return records


def save_records(records: list[RunRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([asdict(r) for r in records], indent=1, default=_json_default))


def load_records(path: Path) -> list[RunRecord]:
    return [RunRecord(**d) for d in json.loads(path.read_text())]


def reference_to_dict(ref: Reference) -> dict:
    return asdict(ref)


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _with_budget(config: DEConfig, budget: int) -> DEConfig:
    """Same configuration with the given budget and no early stop (anytime runs)."""
    return replace(config, max_nfev=int(budget), stall_generations=10**9, xtol=0.0)


def _de_label(c: DEConfig) -> str:
    return ("jDE" if c.self_adaptive else f"DE/{c.strategy}") + f" NP={c.pop_size}"


def _entropy(seed) -> int | None:
    if seed is None:
        return None
    if isinstance(seed, np.random.SeedSequence):
        return int(seed.entropy) if seed.spawn_key == () else int(seed.spawn_key[-1])
    return int(seed)


def _record(method, seed, x, f, nfev, t_start, hist_n, hist_f) -> RunRecord:
    return RunRecord(
        method=method,
        seed=_entropy(seed),
        x=[float(v) for v in np.asarray(x)],
        f=float(f),
        feasible=bool(np.isfinite(f)),
        nfev=int(nfev),
        runtime_s=time.perf_counter() - t_start,
        history_nfev=[int(v) for v in hist_n],
        history_f=[float(v) for v in hist_f],
    )


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))