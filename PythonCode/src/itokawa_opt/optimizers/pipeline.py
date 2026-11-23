"""Level 1 - Global search -> clustering of candidates -> local SQP -> independent check.

Hybrid scheme ("memetic" optimisation):

1. **Global phase** - Differential Evolution explores the whole (t0, T) box and returns
   its final population, which at the end is gathered in one or a few basins.
2. **Candidate selection** - the final population is sorted by fitness and reduced to
   at most ``n_candidates`` points that are mutually farther than ``min_distance`` in the
   normalised [0, 1]^2 box: one representative per basin (greedy "leader" clustering).
   Refining all NP individuals would mostly repeat the same descent.
3. **Local phase** - SQP from every candidate; the best feasible result is kept.
4. **Independent check** - the final transfer is propagated with Kepler's equation from
   (r1, v1) for T and compared with the asteroid position: this tests the Lambert
   solution with code that does not share the Lambert implementation.

The cost reported (``nfev``) is the sum of the evaluations of both phases, so that the
hybrid can be compared with the other methods at equal budget.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from itokawa_opt.astro.elements import car2par, par2car
from itokawa_opt.astro.kepler import propagate_true_anomaly
from itokawa_opt.astro.lambert import lambert
from itokawa_opt.optimizers.de import DEConfig, DEResult, differential_evolution
from itokawa_opt.optimizers.local import LocalResult, refine
from itokawa_opt.problems.direct_transfer import DirectTransfer


@dataclass
class HybridResult:
    x: NDArray[np.float64]  # best [t0, T] [s]
    dv: float  # best total delta-V [km/s]
    constraint: float  # <= 0 when feasible
    nfev: int  # DE + SQP evaluations
    nfev_global: int
    nfev_local: int
    runtime_s: float
    de: DEResult = field(repr=False)
    candidates: NDArray[np.float64] = field(repr=False)  # starting points of SQP
    local: list[LocalResult] = field(repr=False, default_factory=list)

    @property
    def feasible(self) -> bool:
        return self.constraint <= 1e-9

    @property
    def history_nfev(self) -> NDArray[np.int64]:
        """DE history followed by one point at the end of the local phase."""
        return np.append(self.de.history_nfev, self.nfev)

    @property
    def history_f(self) -> NDArray[np.float64]:
        best = min(self.de.history_f[-1], self.dv if self.feasible else np.inf)
        return np.append(self.de.history_f, best)


def select_candidates(
    population: NDArray,
    f: NDArray,
    violation: NDArray,
    bounds: tuple,
    n_candidates: int = 3,
    min_distance: float = 0.02,
) -> NDArray[np.float64]:
    """Best individuals that are pairwise farther than ``min_distance`` (normalised box).

    Feasible individuals come first (sorted by f), then infeasible ones (sorted by
    violation), consistently with Deb's rules.
    """
    lb, ub = (np.asarray(b, dtype=float) for b in bounds)
    z = (population - lb) / (ub - lb)
    order = np.lexsort((f, violation, violation > 0.0))  # keys: last = primary
    chosen: list[int] = []
    for i in order:
        if all(np.linalg.norm(z[i] - z[j]) > min_distance for j in chosen):
            chosen.append(int(i))
        if len(chosen) == n_candidates:
            break
    return population[chosen]


def hybrid_de_sqp(
    problem: DirectTransfer,
    de_config: DEConfig | None = None,
    seed: int | np.random.SeedSequence | None = None,
    n_candidates: int = 3,
    min_distance: float = 0.02,
    snapshot_generations: tuple[int, ...] = (),
) -> HybridResult:
    """DE on the whole box, then SQP from the best distinct candidates."""
    t_start = time.perf_counter()
    bounds = problem.get_bounds()
    de = differential_evolution(
        problem.fitness_batch,
        bounds,
        de_config,
        seed=seed,
        snapshot_generations=snapshot_generations,
    )
    cand = select_candidates(
        de.population, de.population_f, de.population_violation, bounds, n_candidates, min_distance
    )
    local = [refine(problem, c) for c in cand]
    nfev_local = sum(r.nfev for r in local)

    # Best result: the DE point itself is kept if SQP did not improve it (e.g. SQP failed)
    best_x, best_f, best_g = de.x, de.f, (problem.fitness(de.x)[1] if de.feasible else np.inf)
    for r in local:
        if r.constraint <= 1e-9 and (r.dv < best_f or best_g > 1e-9):
            best_x, best_f, best_g = r.x, r.dv, r.constraint
    return HybridResult(
        x=np.asarray(best_x, dtype=float),
        dv=float(best_f),
        constraint=float(best_g),
        nfev=de.nfev + nfev_local,
        nfev_global=de.nfev,
        nfev_local=nfev_local,
        runtime_s=time.perf_counter() - t_start,
        de=de,
        candidates=cand,
        local=local,
    )


def verify_transfer(problem: DirectTransfer, x) -> dict[str, float]:
    """Independent check of a transfer [t0, T].

    Lambert gives v1 at departure; the arc is propagated for T with Kepler's equation
    (elements -> anomaly propagation -> state) and must end on the asteroid. Returns the
    position miss [km], the delta-V recomputed from the propagated arc and the minimum
    Sun distance constraint.
    """
    t0, tof = float(x[0]), float(x[1])
    mu = problem.mu
    r1, v_dep = problem.departure.rv(t0)
    r2, v_arr = problem.arrival.rv(t0 + tof)
    v1, _ = lambert(r1, r2, tof, mu, problem.prograde)
    el = car2par(r1, v1, mu)
    th_end = propagate_true_anomaly(el.a, el.e, el.theta, tof, mu)
    r_end, v_end = par2car(el.a, el.e, el.i, el.raan, el.argp, th_end, mu)
    return {
        "position_miss_km": float(np.linalg.norm(r_end - r2)),
        "dv_propagated": float(np.linalg.norm(v1 - v_dep) + np.linalg.norm(v_arr - v_end)),
        "dv_lambert": float(problem.fitness([t0, tof])[0]),
        "constraint": float(problem.fitness([t0, tof])[1]),
    }