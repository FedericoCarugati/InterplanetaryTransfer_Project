"""Gradient-based local refinement (SQP) of a grid or global-search solution.

Variables are scaled to [0, 1] inside the optimiser: t0 and T are O(1e7-1e8) s,
and unscaled variables make finite-difference gradients and the SQP Hessian
approximation badly conditioned.

Fair accounting of the cost
---------------------------
SLSQP asks for the objective and for the constraint *separately*, at the same points
(including the finite-difference points of the gradients), and each request would solve
Lambert again. A small cache keyed on the point removes the duplicates, and a counter
reports the number of *distinct* problem evaluations (``nfev``). SciPy's own ``res.nfev``
counts only objective calls and would under-estimate the cost of the local phase in the
benchmark.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from itokawa_opt.problems.direct_transfer import DirectTransfer


@dataclass
class LocalResult:
    x: np.ndarray  # [t0, T] in seconds
    dv: float
    success: bool
    message: str
    nfev: int  # distinct problem evaluations (Lambert solutions)
    constraint: float  # <= 0 when feasible
    nit: int = 0  # SQP iterations


class _CachedFitness:
    """Memoise ``problem.fitness`` on the last few points and count true evaluations."""

    def __init__(self, problem: DirectTransfer, size: int = 32):
        self.problem = problem
        self.size = size
        self.cache: OrderedDict[bytes, list[float]] = OrderedDict()
        self.nfev = 0

    def __call__(self, x: np.ndarray) -> list[float]:
        key = np.asarray(x, dtype=float).tobytes()
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        out = self.problem.fitness(x)
        self.nfev += 1
        self.cache[key] = out
        if len(self.cache) > self.size:
            self.cache.popitem(last=False)
        return out


def refine(
    problem: DirectTransfer,
    x0,
    method: str = "SLSQP",
    tol: float = 1e-12,
    maxiter: int = 500,
) -> LocalResult:
    lb, ub = (np.asarray(b, dtype=float) for b in problem.get_bounds())
    span = ub - lb
    fit = _CachedFitness(problem)

    def to_x(z):
        return lb + np.asarray(z) * span

    def obj(z):
        return fit(to_x(z))[0]

    def ineq(z):  # SciPy convention: >= 0 when feasible
        return -fit(to_x(z))[1]

    z0 = (np.asarray(x0, dtype=float) - lb) / span
    res = minimize(
        obj,
        z0,
        method=method,
        bounds=[(0.0, 1.0)] * len(z0),
        constraints=[{"type": "ineq", "fun": ineq}],
        options={"ftol": tol, "maxiter": maxiter},
    )
    x = to_x(res.x)
    f, g = fit(x)
    return LocalResult(
        x=x,
        dv=f,
        success=bool(res.success),
        message=str(res.message),
        nfev=fit.nfev,
        constraint=g,
        nit=int(getattr(res, "nit", 0)),
    )
