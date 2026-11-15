"""Lambert's problem solver (Izzo, 2015), zero-revolution branch, compiled with Numba.

Reference: D. Izzo, "Revisiting Lambert's problem", Celestial Mechanics and Dynamical
Astronomy 121 (2015) 1-15. The structure follows the reference implementation in
ESA's pykep (Householder iterations on the variable x, Battin series near x = 1).

Design notes
------------
* The scalar kernel :func:`lambert_izzo` is written with plain arithmetic only
  (no NumPy allocations inside the iteration), so the same code can later be
  compiled for the GPU with ``numba.cuda.jit`` (Level 2).
* The iteration count is bounded (``MAX_ITER``) to keep the cost predictable.
* Multi-revolution solutions (M >= 1) are left for Level 1 extensions; the API
  already exposes the ``prograde`` flag.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

MAX_ITER = 15
TOL = 1e-13


@njit(cache=True, fastmath=False)
def _hypergeometric_f(z: float, tol: float) -> float:
    sj, cj, err, j = 1.0, 1.0, 1.0, 0
    while err > tol and j < 100:
        cj1 = cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1.0)
        sj += cj1
        err = abs(cj1)
        cj = cj1
        j += 1
    return sj


@njit(cache=True)
def _x2tof(x: float, lam: float) -> float:
    """Non-dimensional time of flight T(x) for M = 0."""
    battin, lagrange = 0.01, 0.2
    dist = abs(x - 1.0)
    if battin <= dist < lagrange:
        # Lagrange's expression
        a = 1.0 / (1.0 - x * x)
        if a > 0.0:
            alfa = 2.0 * math.acos(x)
            beta = 2.0 * math.asin(math.sqrt(lam * lam / a))
            if lam < 0.0:
                beta = -beta
            return a * math.sqrt(a) * ((alfa - math.sin(alfa)) - (beta - math.sin(beta))) / 2.0
        alfa = 2.0 * math.acosh(x)
        beta = 2.0 * math.asinh(math.sqrt(-lam * lam / a))
        if lam < 0.0:
            beta = -beta
        return -a * math.sqrt(-a) * ((beta - math.sinh(beta)) - (alfa - math.sinh(alfa))) / 2.0

    k = lam * lam
    e = x * x - 1.0
    rho = abs(e)
    z = math.sqrt(1.0 + k * e)
    if dist < battin:
        # Battin's series
        eta = z - lam * x
        s1 = 0.5 * (1.0 - lam - x * eta)
        q = 4.0 / 3.0 * _hypergeometric_f(s1, 1e-11)
        return (eta**3 * q + 4.0 * lam * eta) / 2.0
    # Lancaster-Blanchard
    y = math.sqrt(rho)
    g = x * z - lam * e
    if e < 0.0:
        d = math.acos(g)
    else:
        f = y * (z - lam * x)
        d = math.log(f + g)
    return (x - lam * z - d / y) / e


@njit(cache=True)
def _dtdx(x: float, t: float, lam: float):
    l2 = lam * lam
    l3 = l2 * lam
    umx2 = 1.0 - x * x
    y = math.sqrt(1.0 - l2 * umx2)
    y2 = y * y
    y3 = y2 * y
    dt = 1.0 / umx2 * (3.0 * t * x - 2.0 + 2.0 * l3 * x / y)
    ddt = 1.0 / umx2 * (3.0 * t + 5.0 * x * dt + 2.0 * (1.0 - l2) * l3 / y3)
    dddt = 1.0 / umx2 * (7.0 * x * ddt + 8.0 * dt - 6.0 * (1.0 - l2) * l2 * l3 * x / y3 / y2)
    return dt, ddt, dddt


@njit(cache=True)
def lambert_izzo(r1x, r1y, r1z, r2x, r2y, r2z, tof, mu, prograde):
    """Solve Lambert's problem (zero revolutions).

    Returns ``(v1x, v1y, v1z, v2x, v2y, v2z, n_iter)``; ``n_iter < 0`` flags a
    failure (non-converged or invalid input), in which case velocities are NaN.
    """
    nan = math.nan
    if tof <= 0.0 or mu <= 0.0:
        return nan, nan, nan, nan, nan, nan, -1

    cx, cy, cz = r2x - r1x, r2y - r1y, r2z - r1z
    c = math.sqrt(cx * cx + cy * cy + cz * cz)
    r1 = math.sqrt(r1x * r1x + r1y * r1y + r1z * r1z)
    r2 = math.sqrt(r2x * r2x + r2y * r2y + r2z * r2z)
    s = 0.5 * (r1 + r2 + c)

    i1x, i1y, i1z = r1x / r1, r1y / r1, r1z / r1
    i2x, i2y, i2z = r2x / r2, r2y / r2, r2z / r2
    hx = i1y * i2z - i1z * i2y
    hy = i1z * i2x - i1x * i2z
    hz = i1x * i2y - i1y * i2x
    hn = math.sqrt(hx * hx + hy * hy + hz * hz)
    if hn < 1e-14:  # collinear position vectors: transfer plane undefined
        return nan, nan, nan, nan, nan, nan, -2
    hx, hy, hz = hx / hn, hy / hn, hz / hn

    lam = math.sqrt(max(0.0, 1.0 - c / s))
    if hz < 0.0:  # transfer angle > 180 deg for a prograde motion
        lam = -lam
        t1x, t1y, t1z = i1y * hz - i1z * hy, i1z * hx - i1x * hz, i1x * hy - i1y * hx
        t2x, t2y, t2z = i2y * hz - i2z * hy, i2z * hx - i2x * hz, i2x * hy - i2y * hx
    else:
        t1x, t1y, t1z = hy * i1z - hz * i1y, hz * i1x - hx * i1z, hx * i1y - hy * i1x
        t2x, t2y, t2z = hy * i2z - hz * i2y, hz * i2x - hx * i2z, hx * i2y - hy * i2x
    if not prograde:
        lam = -lam
        t1x, t1y, t1z = -t1x, -t1y, -t1z
        t2x, t2y, t2z = -t2x, -t2y, -t2z

    t_nd = math.sqrt(2.0 * mu / s**3) * tof

    # Initial guess (Izzo 2015 / pykep)
    l2 = lam * lam
    l3 = l2 * lam
    t00 = math.acos(lam) + lam * math.sqrt(1.0 - l2)
    t1 = 2.0 / 3.0 * (1.0 - l3)
    if t_nd >= t00:
        x = -(t_nd - t00) / (t_nd - t00 + 4.0)
    elif t_nd <= t1:
        x = t1 * (t1 - t_nd) / (2.0 / 5.0 * (1.0 - l2 * l3) * t_nd) + 1.0
    else:
        x = (t_nd / t00) ** (0.69314718055994529 / math.log(t1 / t00)) - 1.0

    # Householder iterations
    n_iter = -3
    for it in range(MAX_ITER):
        tof_x = _x2tof(x, lam)
        dt, ddt, dddt = _dtdx(x, tof_x, lam)
        delta = tof_x - t_nd
        dt2 = dt * dt
        x_new = x - delta * (dt2 - delta * ddt / 2.0) / (
            dt * (dt2 - delta * ddt) + dddt * delta * delta / 6.0
        )
        err = abs(x - x_new)
        x = x_new
        if err < TOL:
            n_iter = it + 1
            break
    if n_iter < 0:
        return nan, nan, nan, nan, nan, nan, n_iter

    # Velocities
    gamma = math.sqrt(mu * s / 2.0)
    rho = (r1 - r2) / c
    sigma = math.sqrt(max(0.0, 1.0 - rho * rho))
    y = math.sqrt(1.0 - l2 + l2 * x * x)
    vr1 = gamma * ((lam * y - x) - rho * (lam * y + x)) / r1
    vr2 = -gamma * ((lam * y - x) + rho * (lam * y + x)) / r2
    vt = gamma * sigma * (y + lam * x)
    vt1, vt2 = vt / r1, vt / r2
    return (
        vr1 * i1x + vt1 * t1x,
        vr1 * i1y + vt1 * t1y,
        vr1 * i1z + vt1 * t1z,
        vr2 * i2x + vt2 * t2x,
        vr2 * i2y + vt2 * t2y,
        vr2 * i2z + vt2 * t2z,
        n_iter,
    )


def lambert(r1, r2, tof: float, mu: float, prograde: bool = True):
    """Convenience wrapper: NumPy vectors in, ``(v1, v2)`` NumPy vectors out."""
    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float)
    out = lambert_izzo(
        r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], float(tof), float(mu), bool(prograde)
    )
    if out[6] < 0:
        raise ValueError(f"Lambert solver failed (code {out[6]}).")
    return np.array(out[0:3]), np.array(out[3:6])


@njit(parallel=True, cache=True)
def lambert_batch(r1, r2, tof, mu, prograde):
    """Solve many independent Lambert problems in parallel on the CPU.

    ``r1``, ``r2``: (n, 3) arrays; ``tof``: (n,) array. Returns ``v1``, ``v2``
    (n, 3) and the iteration counts (n,). This is the CPU twin of the GPU kernel.
    """
    n = tof.shape[0]
    v1 = np.empty((n, 3))
    v2 = np.empty((n, 3))
    iters = np.empty(n, dtype=np.int64)
    for k in prange(n):
        out = lambert_izzo(
            r1[k, 0], r1[k, 1], r1[k, 2], r2[k, 0], r2[k, 1], r2[k, 2], tof[k], mu, prograde
        )
        v1[k, 0], v1[k, 1], v1[k, 2] = out[0], out[1], out[2]
        v2[k, 0], v2[k, 1], v2[k, 2] = out[3], out[4], out[5]
        iters[k] = out[6]
    return v1, v2, iters
