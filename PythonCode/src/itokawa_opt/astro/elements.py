"""Conversions between Keplerian elements and Cartesian state.

Python counterparts of the thesis MATLAB functions ``par2car`` and ``car2par``.
Rotation convention (same as the MATLAB code): the matrix

    T = R3(omega) @ R1(i) @ R3(Omega)

maps inertial (heliocentric ecliptic) components to perifocal components; its
rows are therefore the perifocal unit vectors P, Q, W expressed in the inertial
frame.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from itokawa_opt.astro.kepler import TWO_PI


class Elements(NamedTuple):
    """Classical orbital elements (km, rad)."""

    a: float
    e: float
    i: float
    raan: float  # Omega, right ascension of the ascending node
    argp: float  # omega, argument of pericentre
    theta: float  # true anomaly


def inertial_to_perifocal(i: ArrayLike, raan: ArrayLike, argp: ArrayLike) -> NDArray[np.float64]:
    """Rotation matrix (..., 3, 3) from inertial to perifocal frame."""
    i, raan, argp = (np.asarray(x, dtype=float) for x in (i, raan, argp))
    co, so = np.cos(raan), np.sin(raan)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(argp), np.sin(argp)
    # Closed form of R3(w) @ R1(i) @ R3(O)
    t = np.empty(np.broadcast(i, raan, argp).shape + (3, 3))
    t[..., 0, 0] = cw * co - sw * ci * so
    t[..., 0, 1] = cw * so + sw * ci * co
    t[..., 0, 2] = sw * si
    t[..., 1, 0] = -sw * co - cw * ci * so
    t[..., 1, 1] = -sw * so + cw * ci * co
    t[..., 1, 2] = cw * si
    t[..., 2, 0] = si * so
    t[..., 2, 1] = -si * co
    t[..., 2, 2] = ci
    return t


def par2car(
    a: ArrayLike,
    e: ArrayLike,
    i: ArrayLike,
    raan: ArrayLike,
    argp: ArrayLike,
    theta: ArrayLike,
    mu: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Keplerian elements -> position and velocity, shape (..., 3) [km, km/s]."""
    a, e, theta = (np.asarray(x, dtype=float) for x in (a, e, theta))
    p = a * (1.0 - e**2)
    r = p / (1.0 + e * np.cos(theta))
    r_pf = np.stack([r * np.cos(theta), r * np.sin(theta), np.zeros_like(r)], axis=-1)
    k = np.sqrt(mu / p)
    v_pf = np.stack([-k * np.sin(theta), k * (e + np.cos(theta)), np.zeros_like(r)], axis=-1)
    t = inertial_to_perifocal(i, raan, argp)
    # perifocal -> inertial uses the transpose: x_in = T^T x_pf
    r_in = np.einsum("...ji,...j->...i", t, r_pf)
    v_in = np.einsum("...ji,...j->...i", t, v_pf)
    return r_in, v_in


def car2par(r: ArrayLike, v: ArrayLike, mu: float) -> Elements:
    """Cartesian state -> Keplerian elements (single state, elliptic orbit)."""
    r = np.asarray(r, dtype=float).reshape(3)
    v = np.asarray(v, dtype=float).reshape(3)
    rn, vn = np.linalg.norm(r), np.linalg.norm(v)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    i = np.arccos(h[2] / hn)
    e_vec = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(e_vec)
    a = 1.0 / (2.0 / rn - vn**2 / mu)
    k = np.array([0.0, 0.0, 1.0])
    n = np.cross(k, h)
    n = n / np.linalg.norm(n)
    raan = np.arccos(np.clip(n[0], -1.0, 1.0))
    if n[1] < 0.0:
        raan = TWO_PI - raan
    argp = np.arccos(np.clip(np.dot(n, e_vec) / e, -1.0, 1.0))
    if e_vec[2] < 0.0:
        argp = TWO_PI - argp
    theta = np.arccos(np.clip(np.dot(e_vec, r) / (e * rn), -1.0, 1.0))
    if np.dot(r, v) < 0.0:
        theta = TWO_PI - theta
    return Elements(float(a), float(e), float(i), float(raan), float(argp), float(theta))
