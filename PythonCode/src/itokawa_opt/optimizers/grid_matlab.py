"""Faithful, vectorised Python port of the thesis MATLAB grid search (Scenario 2, advanced).

The MATLAB script loops over three indices:

* ``k``  Earth departure time step (3 Earth revolutions, step dt = 1e4 s)
* ``s``  asteroid arrival time step (one asteroid period after departure)
* ``j``  argument of pericentre of the transfer conic (n = 100 values)

and keeps the conics whose time of flight matches the asteroid's travel time
within a tolerance window (rendezvous). Here the ``s`` and ``j`` loops are
replaced by NumPy broadcasting, so a single Python loop over ``k`` remains.

``matlab_compat=True`` reproduces two quirks of the original code so that the
numbers match exactly (see the README of the validation for details):

1. the position/velocity arrays are filled with the state at ``theta(j+1)``
   while the anomaly array is indexed with ``theta(j)`` (one-step shift);
2. the last column of the position/velocity arrays is never written and keeps
   the initial value ``[1, 1, 1]``.

Note on the geometry: the transfer plane normal is ``r_A x r_E`` (as in the
MATLAB code), i.e. only prograde transfers with a transfer angle larger than
180 deg (or retrograde ones) are generated. Lambert's solver in Level 1 will
explore both cases.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from itokawa_opt.astro.elements import car2par, par2car
from itokawa_opt.astro.kepler import (
    TWO_PI,
    period,
    propagate_true_anomaly,
    time_of_flight,
    wrap_2pi,
)
from itokawa_opt.constants import AU, M_ITOKAWA, M_SUN_THESIS, MU_SUN_THESIS


# --------------------------------------------------------------------------------------
# Scenario definition (all numbers copied from the MATLAB script)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ThesisScenario:
    """Inputs of the MATLAB script ``Scenario02_Advanced``."""

    mu: float = MU_SUN_THESIS
    # Earth heliocentric elements (fixed, osculating) [km, -, rad, rad, rad]
    a_earth: float = 1.4946e08
    e_earth: float = 0.016
    i_earth: float = 9.1920e-05
    raan_earth: float = 2.7847
    argp_earth: float = 5.2643
    theta0_earth: float = 0.0  # t = 0 at Earth perihelion (4 Jan 2025)
    # Itokawa heliocentric state at t = 0 from JPL Horizons [km, km/s]
    r0_ast: tuple[float, float, float] = (
        1.257283091380454e08,
        2.179025686550188e08,
        -1.109788941692904e06,
    )
    v0_ast: tuple[float, float, float] = (-17.36883104272374, 8.840929257352350, 0.5486875533296964)
    # Discretisation
    dt: float = 10_000.0  # time step [s]
    n_omega: int = 100  # number of transfer arguments of pericentre
    n_earth_revs: int = 3  # search window = 3 Earth revolutions
    # Constraints / tolerances
    r_min_sun: float = 0.7 * AU  # minimum distance from the Sun [km]
    window_scale: float = 100.0  # MATLAB: Deltat_rif = 100 * timeWindow(...)
    epoch0: str = "2025-01-04"  # calendar date of t = 0 (Earth perihelion)


@dataclass
class GridResult:
    """Output of the grid search."""

    dv_min: float
    k_opt: int  # 1-based, as in MATLAB
    s_opt: int  # 1-based
    j_opt: int  # 1-based
    rev_opt: int  # Earth revolution (0, 1, 2)
    t_dep: float  # departure time from t = 0 [s]
    tof: float  # transfer time of flight [s]
    dv1: float
    dv2: float
    transfer: dict[str, float]  # a, e, i, raan, argp, theta1, theta2
    theta_earth_label: float  # "thetai_T" as stored by MATLAB
    theta_ast_label: float  # "thetaf_A" as stored by MATLAB
    dv_min_per_rev: list[float]
    k_opt_per_rev: list[int]
    n_feasible: int
    n_evaluated: int
    runtime_s: float
    dv_map: list[NDArray[np.float64]] = field(repr=False, default_factory=list)
    theta_earth_grid: NDArray[np.float64] = field(repr=False, default_factory=lambda: np.empty(0))
    omega_grid: NDArray[np.float64] = field(repr=False, default_factory=lambda: np.empty(0))


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def time_window(a: float, e: float, theta: float, r_soi: float, mu: float) -> float:
    """
    Time needed to travel a distance r_soi at the orbital speed at ``theta`` [s].
    """
    p = a * (1.0 - e**2)
    vr = np.sqrt(mu / p) * e * np.sin(theta)
    vt = np.sqrt(mu / p) * (1.0 + e * np.cos(theta))
    return float(r_soi / np.hypot(vr, vt))


def _discretise_orbit(a, e, i, raan, argp, theta0, t_end, dt, mu, matlab_compat):
    """True anomaly samples and states every dt seconds (MATLAB ``t_vect``)."""
    t = np.arange(0.0, t_end + 0.5 * dt * 1e-9, dt)  # = 0:dt:t_end
    if t[-1] > t_end:
        t = t[:-1]
    theta = propagate_true_anomaly(a, e, theta0, t, mu)
    theta[0] = theta0
    r, v = par2car(a, e, i, raan, argp, theta, mu)
    if matlab_compat:
        # MATLAB stores state(theta(j+1)) in column j and never fills the last one
        r = np.vstack([r[1:], np.ones((1, 3))])
        v = np.vstack([v[1:], np.ones((1, 3))])
    return theta, r, v


def _wrap_index(idx: NDArray[np.int64] | int, n: int):
    """MATLAB pattern ``if idx > n: if idx > 2n: idx-2n else idx-n`` (1-based)."""
    idx = np.asarray(idx)
    return np.where(idx > 2 * n, idx - 2 * n, np.where(idx > n, idx - n, idx))


# --------------------------------------------------------------------------------------
# Main routine
# --------------------------------------------------------------------------------------
def run_grid_search(
    sc: ThesisScenario | None = None,
    matlab_compat: bool = True,
    k_range: tuple[int, int] | None = None,
    verbose: bool = True,
) -> GridResult:
    """Run the thesis grid search. ``k_range`` (1-based, inclusive) restricts the
    departure indices, useful for quick tests."""
    sc = sc or ThesisScenario()
    t_start = time.perf_counter()
    mu = sc.mu

    # Asteroid elements from the Horizons state (MATLAB: car2par)
    ast = car2par(sc.r0_ast, sc.v0_ast, mu)
    a_a, e_a = ast.a, ast.e
    p_a = a_a * (1.0 - e_a**2)

    t_ast = float(period(a_a, mu))
    t_earth = float(period(sc.a_earth, mu))

    th_e, r_e, v_e = _discretise_orbit(
        sc.a_earth,
        sc.e_earth,
        sc.i_earth,
        sc.raan_earth,
        sc.argp_earth,
        sc.theta0_earth,
        t_earth,
        sc.dt,
        mu,
        matlab_compat,
    )
    th_a, r_a, v_a = _discretise_orbit(
        a_a, e_a, ast.i, ast.raan, ast.argp, ast.theta, t_ast, sc.dt, mu, matlab_compat
    )
    n_e, n_a = th_e.size, th_a.size

    # Rendezvous tolerance window (MATLAB: Deltat)
    r_soi_1 = p_a / (1.0 + e_a) * (M_ITOKAWA / M_SUN_THESIS) ** 0.4
    r_soi_2 = p_a / (1.0 - e_a) * (M_ITOKAWA / M_SUN_THESIS) ** 0.4
    r_soi_med = 0.5 * (r_soi_1 + r_soi_2)
    dt_ref = sc.window_scale * time_window(a_a, e_a, np.pi, r_soi_med, mu)
    v_ref = np.sqrt(mu / p_a) * (1.0 - e_a)
    v_ast = np.sqrt(mu / p_a) * np.hypot(e_a * np.sin(th_a), 1.0 + e_a * np.cos(th_a))
    window = (v_ref / v_ast) * dt_ref

    omega = np.linspace(0.0, TWO_PI, sc.n_omega)[None, :]  # (1, n)
    cos_w, sin_w = np.cos(omega), np.sin(omega)

    s_idx = np.arange(1, n_a + 1)  # 1-based s
    k_first, k_last = (1, sc.n_earth_revs * n_e) if k_range is None else k_range

    best = dict(dv=np.inf)
    dv_min_rev = [np.inf] * sc.n_earth_revs
    k_opt_rev = [-1] * sc.n_earth_revs
    dv_map = [np.full((n_e, sc.n_omega), np.nan) for _ in range(sc.n_earth_revs)]
    n_feasible = 0

    # Degenerate geometries (e.g. the never-filled [1,1,1] column) give NaN: silenced
    with np.errstate(divide="ignore", invalid="ignore"):
        for k in range(k_first, k_last + 1):
            rev = (k - 1) // n_e  # 0, 1, 2
            ke = k - rev * n_e  # Earth index (1-based)
            r1v, v1e = r_e[ke - 1], v_e[ke - 1]
            th_a_i = th_a[int(_wrap_index(k, n_a)) - 1]

            f_idx = _wrap_index(s_idx + k, n_a) - 1  # 0-based, (N,)
            r2v, v2a = r_a[f_idx], v_a[f_idx]
            tol = window[f_idx]

            # Transfer plane: h ~ r_A x r_E (MATLAB: cross(r_vect_A_f, r_vect_T))
            h = np.cross(r2v, r1v)
            h /= np.linalg.norm(h, axis=1, keepdims=True)
            node = np.stack([-h[:, 1], h[:, 0], np.zeros(n_a)], axis=1)
            node /= np.linalg.norm(node, axis=1, keepdims=True)
            m_hat = np.cross(h, node)
            i_t = np.arccos(np.clip(h[:, 2], -1.0, 1.0))
            raan_t = np.arccos(np.clip(node[:, 0], -1.0, 1.0))
            raan_t = np.where(node[:, 1] >= 0.0, raan_t, TWO_PI - raan_t)

            # Argument of latitude of the two points, then true anomalies on the transfer
            u1 = np.arctan2(m_hat @ r1v, node @ r1v)[:, None]  # (N, 1)
            u2 = np.arctan2(np.einsum("ij,ij->i", m_hat, r2v), np.einsum("ij,ij->i", node, r2v))[
                :, None
            ]
            th1 = wrap_2pi(u1 - omega)  # (N, n)
            th2 = wrap_2pi(u2 - omega)
            rn1 = np.linalg.norm(r1v)
            rn2 = np.linalg.norm(r2v, axis=1)[:, None]

            with np.errstate(divide="ignore", invalid="ignore"):
                e_t = (rn1 - rn2) / (rn2 * np.cos(th2) - rn1 * np.cos(th1))
                ok = (e_t >= 0.0) & (e_t < 1.0)
                a_t = rn1 * (1.0 + e_t * np.cos(th1)) / (1.0 - e_t**2)
                p_t = a_t * (1.0 - e_t**2)
                # Minimum Sun distance along the arc (same logic as MATLAB)
                r_end1 = p_t / (1.0 + e_t * np.cos(th1))
                r_end2 = p_t / (1.0 + e_t * np.cos(th2))
                r_min = np.where(
                    th2 > th1, np.where((TWO_PI - th2) > th1, r_end1, r_end2), a_t * (1.0 - e_t)
                )
                ok &= r_min > sc.r_min_sun
                if not ok.any():
                    continue
                tof_t = np.where(
                    ok,
                    time_of_flight(np.where(ok, a_t, 1.0), np.where(ok, e_t, 0.0), th1, th2, mu),
                    np.nan,
                )
            tof_a = time_of_flight(a_a, e_a, th_a_i, th_a[f_idx], mu)[:, None]
            ok &= (tof_a + tol[:, None] > tof_t) & (tof_t > tof_a - tol[:, None])
            if not ok.any():
                continue

            # Delta-V only where feasible
            si, ji = np.nonzero(ok)
            e_s, p_s = e_t[si, ji], p_t[si, ji]
            cw, sw = cos_w[0, ji], sin_w[0, ji]
            pv = cw[:, None] * node[si] + sw[:, None] * m_hat[si]  # perifocal P
            qv = -sw[:, None] * node[si] + cw[:, None] * m_hat[si]  # perifocal Q
            kk = np.sqrt(mu / p_s)[:, None]
            t1, t2 = th1[si, ji][:, None], th2[si, ji][:, None]
            vt1 = kk * (-np.sin(t1) * pv + (e_s[:, None] + np.cos(t1)) * qv)
            vt2 = kk * (-np.sin(t2) * pv + (e_s[:, None] + np.cos(t2)) * qv)
            dv1 = np.linalg.norm(vt1 - v1e, axis=1)
            dv2 = np.linalg.norm(v2a[si] - vt2, axis=1)
            dv = dv1 + dv2
            n_feasible += dv.size

            # min over s for each omega (MATLAB map keeps the last s: we keep the minimum)
            row = dv_map[rev][ke - 1]
            np.fmin.at(row, ji, dv)

            b = int(np.argmin(dv))
            if dv[b] < dv_min_rev[rev]:
                dv_min_rev[rev], k_opt_rev[rev] = float(dv[b]), k
            if dv[b] < best["dv"]:
                s_b, j_b = int(si[b]), int(ji[b])
                best = dict(
                    dv=float(dv[b]),
                    dv1=float(dv1[b]),
                    dv2=float(dv2[b]),
                    k=k,
                    s=s_b + 1,
                    j=j_b + 1,
                    rev=rev,
                    ke=ke,
                    tof=float(tof_t[s_b, j_b]),
                    transfer=dict(
                        a=float(a_t[s_b, j_b]),
                        e=float(e_t[s_b, j_b]),
                        i=float(i_t[s_b]),
                        raan=float(raan_t[s_b]),
                        argp=float(omega[0, j_b]),
                        theta1=float(th1[s_b, j_b]),
                        theta2=float(th2[s_b, j_b]),
                    ),
                    theta_e=float(th_e[ke - 1]),
                    theta_a=float(th_a[f_idx[s_b]]),
                )
            if verbose and k % 1000 == 0:
                print(
                    f"  k = {k:5d}/{k_last}  best dV = {best['dv']:.6f} km/s  "
                    f"({time.perf_counter() - t_start:.0f} s)",
                    flush=True,
                )

    if not np.isfinite(best["dv"]):
        raise RuntimeError("No feasible transfer found in the requested range.")

    # Departure time from t = 0 (the Earth state used is at index ke, shifted in compat mode)
    shift = 1 if matlab_compat else 0
    t_dep = (best["ke"] - 1 + shift) * sc.dt + best["rev"] * t_earth

    return GridResult(
        dv_min=best["dv"],
        k_opt=best["k"],
        s_opt=best["s"],
        j_opt=best["j"],
        rev_opt=best["rev"],
        t_dep=t_dep,
        tof=best["tof"],
        dv1=best["dv1"],
        dv2=best["dv2"],
        transfer=best["transfer"],
        theta_earth_label=best["theta_e"],
        theta_ast_label=best["theta_a"],
        dv_min_per_rev=dv_min_rev,
        k_opt_per_rev=k_opt_rev,
        n_feasible=n_feasible,
        n_evaluated=(k_last - k_first + 1) * n_a * sc.n_omega,
        runtime_s=time.perf_counter() - t_start,
        dv_map=dv_map,
        theta_earth_grid=th_e,
        omega_grid=omega[0],
    )
