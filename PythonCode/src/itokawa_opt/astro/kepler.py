"""Two-body Keplerian relations for elliptic orbits (vectorised with NumPy).

Python counterparts of the thesis MATLAB functions:

* ``TOF(a, e, th1, th2, mu)``          -> :func:`time_of_flight`
* ``KepleroDiretto(a, e, th, dt, mu)`` -> :func:`propagate_true_anomaly`

All functions accept scalars or NumPy arrays (broadcasting rules apply).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

TWO_PI = 2.0 * np.pi


def wrap_2pi(angle: ArrayLike) -> NDArray[np.float64]:
    """Map an angle to the interval [0, 2*pi)."""
    return np.mod(np.asarray(angle, dtype=float), TWO_PI)


def period(a: ArrayLike, mu: float) -> NDArray[np.float64]:
    """Orbital period of an ellipse [s]."""
    return TWO_PI * np.sqrt(np.asarray(a, dtype=float) ** 3 / mu)


def mean_motion(a: ArrayLike, mu: float) -> NDArray[np.float64]:
    """Mean motion n = sqrt(mu / a^3) [rad/s]."""
    return np.sqrt(mu / np.asarray(a, dtype=float) ** 3)


def true_to_eccentric(theta: ArrayLike, e: ArrayLike) -> NDArray[np.float64]:
    """True anomaly -> eccentric anomaly, result in [0, 2*pi)."""
    theta = np.asarray(theta, dtype=float)
    e = np.asarray(e, dtype=float)
    ecc = 2.0 * np.arctan2(
        np.sqrt(1.0 - e) * np.sin(theta / 2.0), np.sqrt(1.0 + e) * np.cos(theta / 2.0)
    )
    return wrap_2pi(ecc)


def eccentric_to_true(ecc: ArrayLike, e: ArrayLike) -> NDArray[np.float64]:
    """Eccentric anomaly -> true anomaly, result in [0, 2*pi)."""
    ecc = np.asarray(ecc, dtype=float)
    e = np.asarray(e, dtype=float)
    theta = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(ecc / 2.0), np.sqrt(1.0 - e) * np.cos(ecc / 2.0)
    )
    return wrap_2pi(theta)


def eccentric_to_mean(ecc: ArrayLike, e: ArrayLike) -> NDArray[np.float64]:
    """Kepler's equation M = E - e sin E."""
    ecc = np.asarray(ecc, dtype=float)
    return ecc - np.asarray(e, dtype=float) * np.sin(ecc)


def solve_kepler(
    mean_anom: ArrayLike, e: ArrayLike, tol: float = 1e-14, max_iter: int = 50
) -> NDArray[np.float64]:
    """Solve M = E - e sin E for E with Newton's method (vectorised).

    A fixed iteration budget with a convergence mask keeps the routine
    branch-free.
    """
    m = wrap_2pi(mean_anom)
    e = np.broadcast_to(np.asarray(e, dtype=float), m.shape)
    ecc = np.where(e < 0.8, m, np.pi * np.ones_like(m))  # standard starting guess
    for _ in range(max_iter):
        f = ecc - e * np.sin(ecc) - m
        step = f / (1.0 - e * np.cos(ecc))
        ecc = ecc - step
        if np.all(np.abs(step) < tol):
            break
    return wrap_2pi(ecc)


def time_of_flight(
    a: ArrayLike, e: ArrayLike, theta1: ArrayLike, theta2: ArrayLike, mu: float
) -> NDArray[np.float64]:
    """Time to go from true anomaly theta1 to theta2 moving forward [s].

    The result is in [0, T): if theta2 precedes theta1 the pericentre passage is
    included (one orbital period is added), as in the thesis ``TOF.m``.
    """
    m1 = eccentric_to_mean(true_to_eccentric(theta1, e), e)
    m2 = eccentric_to_mean(true_to_eccentric(theta2, e), e)
    dt = (m2 - m1) / mean_motion(a, mu)
    return np.where(dt < 0.0, dt + period(a, mu), dt)


def propagate_true_anomaly(
    a: ArrayLike, e: ArrayLike, theta0: ArrayLike, dt: ArrayLike, mu: float
) -> NDArray[np.float64]:
    """True anomaly reached after a time dt starting from theta0 (``KepleroDiretto``)."""
    m0 = eccentric_to_mean(true_to_eccentric(theta0, e), e)
    m = m0 + mean_motion(a, mu) * np.asarray(dt, dtype=float)
    return eccentric_to_true(solve_kepler(m, e), e)
