"""Porkchop grid search on (t0, T) with Lambert's solver.

This replaces the three nested MATLAB loops (Earth position, asteroid position,
argument of pericentre): for each (t0, T) Lambert's solver returns the *unique*
zero-revolution conic that reaches the asteroid in exactly T, so the rendezvous
is exact and no tolerance window is needed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from itokawa_opt.problems.direct_transfer import DirectTransfer


@dataclass
class PorkchopResult:
    t0: NDArray[np.float64]  # (n0,) departure times [s]
    tof: NDArray[np.float64]  # (nT,) times of flight [s]
    dv: NDArray[np.float64]  # (n0, nT) total delta-V [km/s], NaN if infeasible
    dv_raw: NDArray[np.float64] = field(repr=False)  # without constraints
    r_min: NDArray[np.float64] = field(repr=False)
    angle: NDArray[np.float64] = field(repr=False)
    best_t0: float = np.nan
    best_tof: float = np.nan
    best_dv: float = np.nan
    n_evaluated: int = 0
    runtime_s: float = 0.0


def porkchop(
    problem: DirectTransfer,
    t0: NDArray[np.float64],
    tof: NDArray[np.float64],
    long_way_only: bool = False,
) -> PorkchopResult:
    """Evaluate the full (t0, T) grid and return the best feasible point.

    ``long_way_only=True`` keeps only transfer angles > 180 deg, i.e. the subset
    of conics that the thesis MATLAB code could generate (useful for a like-for-like
    comparison).
    """
    t_start = time.perf_counter()
    t0_g, tof_g = np.meshgrid(t0, tof, indexing="ij")
    res = problem.evaluate(t0_g, tof_g)
    feasible = np.isfinite(res["dv"]) & (res["r_min"] > problem.r_min_sun)
    if long_way_only:
        feasible &= res["angle"] > np.pi
    dv = np.where(feasible, res["dv"], np.nan)
    if np.all(np.isnan(dv)):
        raise RuntimeError("No feasible point on the grid.")
    i, j = np.unravel_index(np.nanargmin(dv), dv.shape)
    return PorkchopResult(
        t0=t0,
        tof=tof,
        dv=dv,
        dv_raw=res["dv"],
        r_min=res["r_min"],
        angle=res["angle"],
        best_t0=float(t0[i]),
        best_tof=float(tof[j]),
        best_dv=float(dv[i, j]),
        n_evaluated=dv.size,
        runtime_s=time.perf_counter() - t_start,
    )
