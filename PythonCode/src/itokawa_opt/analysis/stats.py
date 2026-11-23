"""Level 1 - Success rate, ECDF of evaluations-to-target, Mann-Whitney U tests.

Conventions: values are *errors* e = f - f_ref >= 0 (lower is better); a run is a
*success* when e <= tol. Failed runs have e = inf where needed (e.g. infeasible).

Why non-parametric statistics: the final errors of a stochastic optimiser are not
Gaussian. Runs that find the global basin pile up near zero, runs that stop in a
secondary basin sit at the gap between basins (here ~0.9 km/s): the distribution is
skewed or bimodal, so mean and standard deviation describe it badly. Median, quartiles,
success rates and rank-based tests (Mann-Whitney) make no assumption on its shape.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import stats as sps


# --------------------------------------------------------------------------------------
# Descriptive statistics
# --------------------------------------------------------------------------------------
@dataclass
class Summary:
    n: int
    best: float
    q1: float
    median: float
    q3: float
    worst: float
    mean: float
    std: float


def summarize(values: ArrayLike) -> Summary:
    v = np.asarray(values, dtype=float)
    fin = v[np.isfinite(v)]
    q1, med, q3 = np.percentile(v, [25, 50, 75]) if v.size else (np.nan,) * 3
    return Summary(
        n=int(v.size),
        best=float(np.min(v)),
        q1=float(q1),
        median=float(med),
        q3=float(q3),
        worst=float(np.max(v)),
        mean=float(fin.mean()) if fin.size else np.inf,
        std=float(fin.std(ddof=1)) if fin.size > 1 else 0.0,
    )


# --------------------------------------------------------------------------------------
# Success rate with a confidence interval
# --------------------------------------------------------------------------------------
def wilson_interval(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n.

    The textbook interval p +- z sqrt(p(1-p)/n) collapses to [1, 1] when every run
    succeeds (p = 1) and can leave [0, 1]; Wilson's interval stays inside [0, 1] and
    still says "with 30/30 successes the true rate is only known to be > 0.89".
    """
    if n == 0:
        return 0.0, 1.0
    z = sps.norm.ppf(0.5 + confidence / 2)
    p = k / n
    den = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    lo = 0.0 if k == 0 else max(0.0, centre - half)
    hi = 1.0 if k == n else min(1.0, centre + half)
    return float(lo), float(hi)


def success_rate(errors: ArrayLike, tol: float) -> tuple[float, float, float]:
    """(rate, low, high) of the runs with error <= tol."""
    e = np.asarray(errors, dtype=float)
    k = int(np.sum(e <= tol))
    lo, hi = wilson_interval(k, e.size)
    return k / e.size if e.size else 0.0, lo, hi


# --------------------------------------------------------------------------------------
# Fixed-target analysis
# --------------------------------------------------------------------------------------
def evaluations_to_target(history_nfev: ArrayLike, history_f: ArrayLike, target: float) -> float:
    """First number of evaluations at which the best-so-far value is <= target (inf if never)."""
    h_n = np.asarray(history_nfev, dtype=float)
    hit = np.nonzero(np.asarray(history_f, dtype=float) <= target)[0]
    return float(h_n[hit[0]]) if hit.size else np.inf


def expected_running_time(evals: ArrayLike, budgets: ArrayLike) -> float:
    """ERT (Hansen et al., COCO platform): expected evaluations to reach the target.

    ``evals``: evaluations-to-target of each run (inf if it failed); ``budgets``:
    evaluations actually spent by each run. ERT = (sum of evaluations spent by *all*
    runs until success or end) / (number of successes). It is what a "restart until
    success" strategy would cost on average.
    """
    e = np.asarray(evals, dtype=float)
    b = np.asarray(budgets, dtype=float)
    ok = np.isfinite(e)
    if not ok.any():
        return np.inf
    return float((e[ok].sum() + b[~ok].sum()) / ok.sum())


def ecdf(values: ArrayLike) -> tuple[NDArray, NDArray]:
    """Empirical CDF: x sorted, y = fraction of runs with value <= x (inf values excluded
    from x, so the curve saturates at the success fraction)."""
    v = np.sort(np.asarray(values, dtype=float))
    y = np.arange(1, v.size + 1) / v.size
    fin = np.isfinite(v)
    return v[fin], y[fin]


# --------------------------------------------------------------------------------------
# Comparison of two methods
# --------------------------------------------------------------------------------------
def vargha_delaney_a12(a: ArrayLike, b: ArrayLike) -> float:
    """Probability that a run of A gives a LOWER value than a run of B (ties count 1/2).

    0.5 = no difference; 1 = A always better; 0 = B always better. Unlike the p-value it
    measures *how large* the difference is, not only whether it is detectable.
    """
    a = np.asarray(a, dtype=float)[:, None]
    b = np.asarray(b, dtype=float)[None, :]
    return float(np.mean((a < b) + 0.5 * (a == b)))


@dataclass
class Comparison:
    name_a: str
    name_b: str
    median_a: float
    median_b: float
    u_statistic: float
    p_value: float
    a12: float  # P(A better than B)

    @property
    def verdict(self) -> str:
        if self.p_value >= 0.05:
            return "no significant difference"
        return f"{self.name_a if self.a12 > 0.5 else self.name_b} is better"


def mann_whitney(a: ArrayLike, b: ArrayLike, name_a: str = "A", name_b: str = "B") -> Comparison:
    """Two-sided Mann-Whitney U test on the errors of two methods (lower is better).

    H0: a run of A and a run of B are equally likely to be the better one. A small p-value
    says the observed difference is unlikely under H0. Infinite values (failures) are
    allowed: they are simply the largest ranks.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    big = np.nanmax(np.where(np.isfinite(np.r_[a, b]), np.r_[a, b], -np.inf))
    big = 1.0 if not np.isfinite(big) else abs(big) + 1.0
    a_r, b_r = np.where(np.isfinite(a), a, 10 * big), np.where(np.isfinite(b), b, 10 * big)
    if np.all(a_r == a_r[0]) and np.all(b_r == a_r[0]):  # identical samples: no test
        u, p = 0.5 * a.size * b.size, 1.0
    else:
        u, p = sps.mannwhitneyu(a_r, b_r, alternative="two-sided")
    return Comparison(
        name_a=name_a,
        name_b=name_b,
        median_a=float(np.median(a)),
        median_b=float(np.median(b)),
        u_statistic=float(u),
        p_value=float(p),
        a12=vargha_delaney_a12(a, b),
    )


def at_resolution(errors: ArrayLike, resolution: float) -> NDArray[np.float64]:
    """Errors below ``resolution`` (and negative round-off) are set to 0.

    Without this, a rank test happily "detects" that 1e-13 km/s beats 1e-11 km/s:
    statistically real, physically meaningless. The resolution states which differences
    we care about (e.g. 1e-6 km/s = 1 mm/s).
    """
    e = np.asarray(errors, dtype=float)
    return np.where(e < resolution, 0.0, e)


@dataclass
class SuccessComparison:
    name_a: str
    name_b: str
    successes_a: int
    successes_b: int
    n_a: int
    n_b: int
    p_value: float


def fisher_success(
    errors_a: ArrayLike, errors_b: ArrayLike, tol: float, name_a: str = "A", name_b: str = "B"
) -> SuccessComparison:
    """Fisher's exact test on the 2x2 table (success / failure) x (method A / method B).

    Answers "is one method more *likely* to find the global minimum?", independently of
    how precisely it converges. Exact (no large-sample approximation), so valid with 30 runs.
    """
    a = np.asarray(errors_a, dtype=float) <= tol
    b = np.asarray(errors_b, dtype=float) <= tol
    table = [[int(a.sum()), int((~a).sum())], [int(b.sum()), int((~b).sum())]]
    _, p = sps.fisher_exact(table)
    return SuccessComparison(name_a, name_b, table[0][0], table[1][0], a.size, b.size, float(p))


def holm_correction(p_values: Sequence[float]) -> list[float]:
    """Holm-Bonferroni adjusted p-values for several comparisons made together.

    With 5 tests at the 5 % level, the chance of at least one false "significant" result
    is up to ~23 %; Holm's method keeps that family-wise error at 5 %.
    """
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    m = p.size
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj.tolist()