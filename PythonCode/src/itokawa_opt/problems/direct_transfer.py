"""Direct Earth -> Itokawa rendezvous transfer as a function of (t0, T).

Decision vector x = [t0, T]:
    t0  departure time from the reference epoch [s]
    T   time of flight [s]

Objective (as in the thesis: heliocentric impulses at departure and arrival)

    J = |v1 - v_E(t0)| + |v_A(t0 + T) - v2|,   (v1, v2) = Lambert(r_E(t0), r_A(t0 + T), T)

Constraint: minimum Sun distance along the transfer arc >= r_min_sun.

The numerical kernel :func:`_evaluate_kernel` works on batches and is the piece
that Level 2 ports to the GPU.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
from numpy.typing import NDArray

from itokawa_opt.astro.lambert import lambert_izzo
from itokawa_opt.constants import AU
from itokawa_opt.ephemeris.thesis_kepler import KeplerEphemeris


@njit(cache=True)
def _arc_min_radius(rx, ry, rz, vx, vy, vz, r2x, r2y, r2z, v2x, v2y, v2z, mu):
    """Minimum distance from the attractor along the arc r1 -> r2 (zero revolutions)."""
    r1 = math.sqrt(rx * rx + ry * ry + rz * rz)
    r2 = math.sqrt(r2x * r2x + r2y * r2y + r2z * r2z)
    hx, hy, hz = ry * vz - rz * vy, rz * vx - rx * vz, rx * vy - ry * vx
    h2 = hx * hx + hy * hy + hz * hz
    # eccentricity vector e = (v x h)/mu - r/|r|
    ex = (vy * hz - vz * hy) / mu - rx / r1
    ey = (vz * hx - vx * hz) / mu - ry / r1
    ez = (vx * hy - vy * hx) / mu - rz / r1
    e = math.sqrt(ex * ex + ey * ey + ez * ez)
    p = h2 / mu
    rp = p / (1.0 + e)
    if e < 1e-12:
        return min(r1, r2)
    # signed true anomalies in (-pi, pi]
    hn = math.sqrt(h2)
    # q = h x e / |h| (perpendicular to e in the orbit plane)
    qx = (hy * ez - hz * ey) / hn
    qy = (hz * ex - hx * ez) / hn
    qz = (hx * ey - hy * ex) / hn
    th1 = math.atan2(rx * qx + ry * qy + rz * qz, rx * ex + ry * ey + rz * ez)
    th2 = math.atan2(r2x * qx + r2y * qy + r2z * qz, r2x * ex + r2y * ey + r2z * ez)
    if e < 1.0:
        two_pi = 2.0 * math.pi
        a1 = th1 % two_pi
        a2 = th2 % two_pi
        passes_peri = a2 < a1
    else:
        passes_peri = th1 < 0.0 < th2
    if passes_peri:
        return rp
    return min(r1, r2)


@njit(parallel=True, cache=True)
def _evaluate_kernel(r1, v_dep_body, r2, v_arr_body, tof, mu, prograde):
    """Batch evaluation. Returns dv1, dv2, r_min, transfer_angle (each (n,))."""
    n = tof.shape[0]
    dv1 = np.empty(n)
    dv2 = np.empty(n)
    rmin = np.empty(n)
    theta = np.empty(n)
    for k in prange(n):
        out = lambert_izzo(
            r1[k, 0], r1[k, 1], r1[k, 2], r2[k, 0], r2[k, 1], r2[k, 2], tof[k], mu, prograde
        )
        if out[6] < 0:
            dv1[k] = math.nan
            dv2[k] = math.nan
            rmin[k] = math.nan
            theta[k] = math.nan
            continue
        a0 = out[0] - v_dep_body[k, 0]
        a1 = out[1] - v_dep_body[k, 1]
        a2 = out[2] - v_dep_body[k, 2]
        dv1[k] = math.sqrt(a0 * a0 + a1 * a1 + a2 * a2)
        b0 = v_arr_body[k, 0] - out[3]
        b1 = v_arr_body[k, 1] - out[4]
        b2 = v_arr_body[k, 2] - out[5]
        dv2[k] = math.sqrt(b0 * b0 + b1 * b1 + b2 * b2)
        rmin[k] = _arc_min_radius(
            r1[k, 0],
            r1[k, 1],
            r1[k, 2],
            out[0],
            out[1],
            out[2],
            r2[k, 0],
            r2[k, 1],
            r2[k, 2],
            out[3],
            out[4],
            out[5],
            mu,
        )
        # prograde transfer angle in [0, 2pi): > pi means "long way"
        cz = r1[k, 0] * r2[k, 1] - r1[k, 1] * r2[k, 0]
        dot = r1[k, 0] * r2[k, 0] + r1[k, 1] * r2[k, 1] + r1[k, 2] * r2[k, 2]
        cx = r1[k, 1] * r2[k, 2] - r1[k, 2] * r2[k, 1]
        cy = r1[k, 2] * r2[k, 0] - r1[k, 0] * r2[k, 2]
        ang = math.atan2(math.sqrt(cx * cx + cy * cy + cz * cz), dot)
        theta[k] = ang if cz >= 0.0 else 2.0 * math.pi - ang
    return dv1, dv2, rmin, theta


@dataclass
class DirectTransfer:
    """Direct rendezvous problem, pygmo-style interface (fitness / get_bounds)."""

    departure: KeplerEphemeris
    arrival: KeplerEphemeris
    t0_bounds: tuple[float, float]
    tof_bounds: tuple[float, float]
    r_min_sun: float = 0.7 * AU
    prograde: bool = True

    @property
    def mu(self) -> float:
        return self.departure.mu

    def get_bounds(self) -> tuple[list[float], list[float]]:
        return ([self.t0_bounds[0], self.tof_bounds[0]], [self.t0_bounds[1], self.tof_bounds[1]])

    def evaluate(self, t0: NDArray[np.float64], tof: NDArray[np.float64]) -> dict[str, NDArray]:
        """Vectorised evaluation on arrays of equal shape."""
        t0 = np.asarray(t0, dtype=float)
        tof = np.asarray(tof, dtype=float)
        shape = np.broadcast(t0, tof).shape
        t0f = np.broadcast_to(t0, shape).ravel()
        tff = np.broadcast_to(tof, shape).ravel()
        r1, v1b = self.departure.rv(t0f)
        r2, v2b = self.arrival.rv(t0f + tff)
        dv1, dv2, rmin, ang = _evaluate_kernel(
            np.ascontiguousarray(r1),
            np.ascontiguousarray(v1b),
            np.ascontiguousarray(r2),
            np.ascontiguousarray(v2b),
            tff,
            self.mu,
            self.prograde,
        )
        out = dict(dv1=dv1, dv2=dv2, dv=dv1 + dv2, r_min=rmin, angle=ang)
        return {k: v.reshape(shape) for k, v in out.items()}

    def fitness(self, x) -> list[float]:
        """pygmo convention: [objective, inequality constraint (<= 0)]."""
        res = self.evaluate(np.array([x[0]]), np.array([x[1]]))
        dv = float(res["dv"][0])
        g = (self.r_min_sun - float(res["r_min"][0])) / AU  # <= 0 when feasible
        return [dv if np.isfinite(dv) else 1e3, 
                g if np.isfinite(g) else 1e3]

    def fitness_batch(self, x: NDArray[np.float64]) -> tuple[NDArray, NDArray]:
        """Objective and constraint for a whole population in one call.

        ``x`` has shape (n, 2), one row per candidate [t0, T]. Returns ``f`` (n,) and
        ``g`` (n,), with the same conventions as :meth:`fitness`: ``g <= 0`` when the arc
        stays farther than ``r_min_sun`` from the Sun, and failed Lambert solutions
        (NaN) are replaced by 1e3 in both, so that they lose every comparison.

        A single call costs ~1 us per candidate, a call of :meth:`fitness` ~250 us
        (Python overhead): population-based optimisers must use this method.
        """
        x = np.atleast_2d(np.asarray(x, dtype=float))
        res = self.evaluate(x[:, 0], x[:, 1])
        f = res["dv"]
        g = (self.r_min_sun - res["r_min"]) / AU
        return np.where(np.isfinite(f), f, 1e3), np.where(np.isfinite(g), g, 1e3)
    
    def get_nic(self) -> int:
        return 1
