"""Level 1 - Differential Evolution written from scratch, vectorised over the population.

Reference: R. Storn, K. Price, "Differential Evolution - A Simple and Efficient Heuristic
for Global Optimization over Continuous Spaces", J. Global Optimization 11 (1997) 341-359.
Self-adaptation (jDE): J. Brest et al., IEEE Trans. Evol. Comput. 10 (2006) 646-657.
Constraint handling: K. Deb, Comput. Methods Appl. Mech. Engrg. 186 (2000) 311-338.

One generation, for every target vector x_i of the population (NP vectors, D variables):

1. mutation      v_i = base + F * (x_r1 - x_r2)          (difference vector)
2. crossover     u_ij = v_ij with probability CR, else x_ij (at least one j from v_i)
3. selection     x_i <- u_i if u_i is not worse than x_i (Deb's feasibility rules)

Strategies (``base`` and number of differences):

* ``rand1bin``          base = x_r0                       (explorative, robust)
* ``best1bin``          base = x_best                     (fast, greedier)
* ``currenttobest1bin`` base = x_i + F * (x_best - x_i)   (compromise)

Design notes
------------
* The optimiser works in the unit hypercube z in [0, 1]^D and maps to the physical
  variables only to call the objective. Here t0 and T are O(1e7-1e8) s: in z every
  variable has the same scale, so F and CR mean the same thing for all of them.
* The update is *generation-synchronous* ("deferred updating"): the trial vectors of a
  whole generation are built from the same population and evaluated with ONE call to
  the batch objective. This is what makes the population a matrix operation (and later
  a GPU kernel). SciPy's default is "immediate updating" (each winner is used at once),
  which cannot be vectorised.
* Constraints are handled with Deb's rules instead of a penalty: no penalty weight to
  tune, and a feasible point always beats an infeasible one.
* Randomness comes only from ``numpy.random.Generator`` built from ``seed``: the same
  seed gives bit-identical runs (reproducibility of the benchmark).
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray

Strategy = Literal["rand1bin", "best1bin", "currenttobest1bin"]

# Batch objective: X (n, D) in physical units -> (f (n,), g (n,) or (n, m) or None).
# g <= 0 means feasible (pygmo convention).
BatchFitness = Callable[[NDArray[np.float64]], tuple[NDArray, NDArray | None]]
ArrayLikeF = NDArray[np.float64] | list[float] | tuple[float, ...]


@dataclass(frozen=True)
class DEConfig:
    """Control parameters of Differential Evolution.

    ``F`` is either a constant or a range (F_min, F_max): in the second case a new value
    is drawn uniformly at every generation ("dither", Storn 1996), which reduces the
    sensitivity to the choice of F. ``self_adaptive=True`` switches to jDE, where every
    individual carries its own F and CR (``F`` and ``CR`` are then only initial values).
    """

    pop_size: int = 20  # NP, absolute number of individuals (>= 4 for rand1)
    strategy: Strategy = "rand1bin"
    F: float | tuple[float, float] = (0.5, 1.0)
    CR: float = 0.9
    self_adaptive: bool = False  # jDE
    max_nfev: int = 10_000  # budget: objective evaluations (initial population included)
    stall_generations: int = 200  # stop if the best value has not improved for this many
    ftol: float = 1e-12  # minimum improvement of the best value that resets the stall count
    xtol: float = 1e-12  # stop when the population collapsed (max spread in z) below this
    init: Literal["lhs", "uniform"] = "lhs"
    # jDE parameters (Brest et al. 2006)
    tau_F: float = 0.1
    tau_CR: float = 0.1
    F_lower: float = 0.1
    F_upper: float = 0.9

    def __post_init__(self) -> None:
        if self.pop_size < 5:
            raise ValueError("pop_size must be >= 5 (target + up to 4 distinct indices).")
        if not 0.0 <= self.CR <= 1.0:
            raise ValueError("CR must be in [0, 1].")
        if self.max_nfev < self.pop_size:
            raise ValueError("max_nfev must be at least pop_size.")


@dataclass
class DEResult:
    x: NDArray[np.float64]  # best point, physical units
    f: float  # best objective value
    violation: float  # constraint violation of the best point (0 = feasible)
    nfev: int  # objective evaluations used
    ngen: int  # generations performed
    message: str  # stopping reason
    runtime_s: float
    # Convergence history: after the initial population and after every generation
    history_nfev: NDArray[np.int64] = field(repr=False)
    history_f: NDArray[np.float64] = field(repr=False)  # best feasible value so far (inf if none)
    # Final population (physical units) and its fitness
    population: NDArray[np.float64] = field(repr=False)
    population_f: NDArray[np.float64] = field(repr=False)
    population_violation: NDArray[np.float64] = field(repr=False)
    # Optional snapshots {generation: population (physical units)}, for plots
    snapshots: dict[int, NDArray[np.float64]] = field(repr=False, default_factory=dict)

    @property
    def feasible(self) -> bool:
        return self.violation <= 0.0


# --------------------------------------------------------------------------------------
# Building blocks (pure functions of arrays, tested individually)
# --------------------------------------------------------------------------------------
def total_violation(g: NDArray | None, n: int) -> NDArray[np.float64]:
    """Sum of the positive parts of the inequality constraints g <= 0, shape (n,)."""
    if g is None:
        return np.zeros(n)
    g = np.asarray(g, dtype=float).reshape(n, -1)
    g = np.where(np.isnan(g), np.inf, g)
    return np.sum(np.maximum(g, 0.0), axis=1)


def deb_better(f_new: NDArray, v_new: NDArray, f_old: NDArray, v_old: NDArray) -> NDArray[np.bool_]:
    """Deb's feasibility rules: is the new point at least as good as the old one?

    1. feasible vs feasible: lower objective wins;
    2. feasible vs infeasible: the feasible one wins;
    3. infeasible vs infeasible: lower total violation wins.

    Ties go to the new point ("<="): it lets the population drift along flat regions.
    NaN objectives are treated as +inf.
    """
    f_new = np.where(np.isnan(f_new), np.inf, f_new)
    f_old = np.where(np.isnan(f_old), np.inf, f_old)
    feas_new, feas_old = v_new <= 0.0, v_old <= 0.0
    both_feasible = feas_new & feas_old
    both_infeasible = ~feas_new & ~feas_old
    return (
        (both_feasible & (f_new <= f_old))
        | (feas_new & ~feas_old)
        | (both_infeasible & (v_new <= v_old))
    )


def best_index(f: NDArray, v: NDArray) -> int:
    """Index of the best individual according to Deb's rules."""
    feasible = v <= 0.0
    if feasible.any():
        f_masked = np.where(feasible & ~np.isnan(f), f, np.inf)
        return int(np.argmin(f_masked))
    return int(np.argmin(v))


def initial_population(
    rng: np.random.Generator, n: int, d: int, method: str = "lhs"
) -> NDArray[np.float64]:
    """Initial points in [0, 1]^d.

    ``lhs`` (Latin hypercube): each variable's range is split in n equal strata and every
    stratum receives exactly one point; the strata are paired at random across variables.
    Coverage is more even than with independent uniform samples.
    """
    if method == "uniform":
        return rng.random((n, d))
    if method != "lhs":
        raise ValueError(f"Unknown init method {method!r}.")
    strata = np.stack([rng.permutation(n) for _ in range(d)], axis=1)  # (n, d)
    return (strata + rng.random((n, d))) / n


def distinct_indices(rng: np.random.Generator, n: int, k: int) -> NDArray[np.int64]:
    """For every row i, k distinct random indices in [0, n) all different from i.

    Vectorised trick: draw an (n, n) matrix of uniform keys, set the diagonal to +inf so
    that i is never selected, and take the k smallest keys of each row. Each row is a
    uniformly random k-subset (in random order) of {0..n-1} without i.
    """
    keys = rng.random((n, n))
    np.fill_diagonal(keys, np.inf)
    return np.argpartition(keys, k, axis=1)[:, :k]


def mutate(
    z: NDArray, idx: NDArray, i_best: int, F: NDArray, strategy: Strategy
) -> NDArray[np.float64]:
    """Mutant vectors for the whole population. ``F`` has shape (n, 1)."""
    if strategy == "rand1bin":
        return z[idx[:, 0]] + F * (z[idx[:, 1]] - z[idx[:, 2]])
    if strategy == "best1bin":
        return z[i_best] + F * (z[idx[:, 0]] - z[idx[:, 1]])
    if strategy == "currenttobest1bin":
        return z + F * (z[i_best] - z) + F * (z[idx[:, 0]] - z[idx[:, 1]])
    raise ValueError(f"Unknown strategy {strategy!r}.")


def binomial_crossover(
    rng: np.random.Generator, z: NDArray, v: NDArray, CR: NDArray
) -> NDArray[np.float64]:
    """u_ij = v_ij if rand_ij <= CR_i or j == j_rand(i), else z_ij. ``CR`` shape (n, 1)."""
    n, d = z.shape
    take = rng.random((n, d)) <= CR
    take[np.arange(n), rng.integers(d, size=n)] = True  # at least one gene from v
    return np.where(take, v, z)


def bounce_back(u: NDArray, z: NDArray) -> NDArray[np.float64]:
    """Repair coordinates outside [0, 1]: midpoint between the parent and the bound.

    Clipping would pile up individuals exactly on the bounds; re-sampling at random
    throws away information. The midpoint keeps the move in the same direction and is
    the usual choice in Price, Storn & Lampinen (2005).
    """
    u = np.where(u < 0.0, 0.5 * z, u)
    return np.where(u > 1.0, 0.5 * (z + 1.0), u)


# --------------------------------------------------------------------------------------
# Main routine
# --------------------------------------------------------------------------------------
def differential_evolution(
    fitness_batch: BatchFitness,
    bounds: tuple[ArrayLikeF, ArrayLikeF],
    config: DEConfig | None = None,
    seed: int | np.random.SeedSequence | None = None,
    snapshot_generations: tuple[int, ...] = (),
    callback: Callable[[int, int, float], None] | None = None,
) -> DEResult:
    """Minimise ``fitness_batch`` over the box ``bounds = (lower, upper)``.

    ``snapshot_generations``: generations at which the population is stored (0 = initial).
    ``callback(generation, nfev, best_f)`` is called after every generation.
    """
    cfg = config or DEConfig()
    t_start = time.perf_counter()
    rng = np.random.default_rng(seed)
    lb, ub = (np.asarray(b, dtype=float) for b in bounds)
    span = ub - lb
    n, d = cfg.pop_size, lb.size

    def evaluate(z: NDArray) -> tuple[NDArray, NDArray]:
        f, g = fitness_batch(lb + z * span)
        return np.asarray(f, dtype=float).reshape(n), total_violation(g, n)

    def best_feasible(f: NDArray, v: NDArray) -> float:
        ok = (v <= 0.0) & ~np.isnan(f)
        return float(f[ok].min()) if ok.any() else np.inf

    # ---------------------------------------------------------------- initial population
    z = initial_population(rng, n, d, cfg.init)
    f, v = evaluate(z)
    nfev = n
    F_ind = np.full(n, cfg.F if np.isscalar(cfg.F) else 0.5, dtype=float)  # jDE only
    CR_ind = np.full(n, cfg.CR)

    hist_nfev, hist_f = [nfev], [best_feasible(f, v)]
    snapshots = {0: lb + z * span} if 0 in snapshot_generations else {}
    best_so_far, stall, gen = hist_f[0], 0, 0
    message = "maximum number of evaluations reached"

    while nfev + n <= cfg.max_nfev:
        gen += 1
        # ------------------------------------------------------------ control parameters
        if cfg.self_adaptive:
            new_F = np.where(
                rng.random(n) < cfg.tau_F, cfg.F_lower + cfg.F_upper * rng.random(n), F_ind
            )
            new_CR = np.where(rng.random(n) < cfg.tau_CR, rng.random(n), CR_ind)
        else:
            f_gen = cfg.F if np.isscalar(cfg.F) else rng.uniform(*cfg.F)  # type: ignore[misc]
            new_F, new_CR = np.full(n, f_gen), np.full(n, cfg.CR)

        # ------------------------------------------------------------ variation
        i_best = best_index(f, v)
        idx = distinct_indices(rng, n, 3)
        mutant = mutate(z, idx, i_best, new_F[:, None], cfg.strategy)
        trial = bounce_back(binomial_crossover(rng, z, mutant, new_CR[:, None]), z)

        # ------------------------------------------------------------ selection
        f_t, v_t = evaluate(trial)
        nfev += n
        win = deb_better(f_t, v_t, f, v)
        z = np.where(win[:, None], trial, z)
        f, v = np.where(win, f_t, f), np.where(win, v_t, v)
        if cfg.self_adaptive:  # jDE: successful parameters survive with their individual
            F_ind, CR_ind = np.where(win, new_F, F_ind), np.where(win, new_CR, CR_ind)

        # ------------------------------------------------------------ bookkeeping
        current = best_feasible(f, v)
        hist_nfev.append(nfev)
        hist_f.append(current)
        if gen in snapshot_generations:
            snapshots[gen] = lb + z * span
        if callback is not None:
            callback(gen, nfev, current)

        if best_so_far - current > cfg.ftol * max(1.0, abs(current)):
            best_so_far, stall = current, 0
        else:
            stall += 1
        if stall >= cfg.stall_generations:
            message = f"no improvement in {cfg.stall_generations} generations"
            break
        if np.max(np.ptp(z, axis=0)) < cfg.xtol:
            message = "population collapsed (spread < xtol)"
            break

    ib = best_index(f, v)
    return DEResult(
        x=lb + z[ib] * span,
        f=float(f[ib]),
        violation=float(v[ib]),
        nfev=nfev,
        ngen=gen,
        message=message,
        runtime_s=time.perf_counter() - t_start,
        history_nfev=np.asarray(hist_nfev, dtype=np.int64),
        history_f=np.asarray(hist_f, dtype=float),
        population=lb + z * span,
        population_f=f,
        population_violation=v,
        snapshots=snapshots,
    )
