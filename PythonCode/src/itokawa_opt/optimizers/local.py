"""Gradient-based local refinement (SQP) of a grid or global-search solution.

Variables are scaled to [0, 1] inside the optimiser: t0 and T are O(1e7-1e8) s,
and unscaled variables make finite-difference gradients and the SQP Hessian
approximation badly conditioned.
"""

from __future__ import annotations

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
    nfev: int
    constraint: float  # <= 0 when feasible


def refine(problem: DirectTransfer, x0, method: str = "SLSQP", tol: float = 1e-12) -> LocalResult:
    lb, ub = (np.asarray(b, dtype=float) for b in problem.get_bounds())
    span = ub - lb

    def to_x(z):
        return lb + np.asarray(z) * span

    def obj(z):
        return problem.fitness(to_x(z))[0]

    def ineq(z):  # SciPy convention: >= 0 when feasible
        return -problem.fitness(to_x(z))[1]

    z0 = (np.asarray(x0, dtype=float) - lb) / span
    res = minimize(
        obj,
        z0,
        method=method,
        bounds=[(0.0, 1.0)] * len(z0),
        constraints=[{"type": "ineq", "fun": ineq}],
        options={"ftol": tol, "maxiter": 500},
    )
    x = to_x(res.x)
    f, g = problem.fitness(x)
    return LocalResult(
        x=x,
        dv=f,
        success=bool(res.success),
        message=str(res.message),
        nfev=int(res.nfev),
        constraint=g,
    )
