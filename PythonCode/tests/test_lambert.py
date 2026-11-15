import numpy as np
import pytest

from itokawa_opt.astro.elements import car2par, par2car
from itokawa_opt.astro.kepler import propagate_true_anomaly
from itokawa_opt.astro.lambert import lambert, lambert_batch

MU = 1.327e11
AU = 1.495978707e8


def _random_cases(n, seed=0):
    rng = np.random.default_rng(seed)
    r1 = rng.normal(size=(n, 3))
    r1 *= (AU * rng.uniform(0.7, 2.0, n) / np.linalg.norm(r1, axis=1))[:, None]
    r2 = rng.normal(size=(n, 3))
    r2 *= (AU * rng.uniform(0.7, 2.0, n) / np.linalg.norm(r2, axis=1))[:, None]
    tof = rng.uniform(30, 1000, n) * 86400.0
    return r1, r2, tof


def test_lambert_hits_target_by_propagation():
    """Propagating (r1, v1) for T with Kepler's equation must land on r2."""
    r1, r2, tof = _random_cases(200)
    for k in range(len(tof)):
        v1, _ = lambert(r1[k], r2[k], tof[k], MU)
        el = car2par(r1[k], v1, MU)
        if el.e >= 1.0:  # hyperbolic arcs are checked against pykep below
            continue
        th = propagate_true_anomaly(el.a, el.e, el.theta, tof[k], MU)
        r_end, _ = par2car(el.a, el.e, el.i, el.raan, el.argp, th, MU)
        assert np.linalg.norm(r_end - r2[k]) < 1e-6 * np.linalg.norm(r2[k])


def test_lambert_matches_pykep():
    pk = pytest.importorskip("pykep")
    r1, r2, tof = _random_cases(500, seed=1)
    for k in range(len(tof)):
        for prograde in (True, False):
            v1, v2 = lambert(r1[k], r2[k], tof[k], MU, prograde=prograde)
            ref = pk.lambert_problem(
                r0=r1[k].tolist(),
                r1=r2[k].tolist(),
                tof=tof[k],
                mu=MU,
                cw=not prograde,
                multi_revs=0,
            )
            assert np.allclose(v1, ref.v0[0], rtol=1e-10, atol=1e-10)
            assert np.allclose(v2, ref.v1[0], rtol=1e-10, atol=1e-10)


def test_batch_equals_scalar():
    r1, r2, tof = _random_cases(1000, seed=2)
    v1b, v2b, it = lambert_batch(r1, r2, tof, MU, True)
    assert np.all(it > 0)
    for k in range(0, 1000, 97):
        v1, v2 = lambert(r1[k], r2[k], tof[k], MU)
        assert np.allclose(v1, v1b[k], rtol=0, atol=1e-12)
        assert np.allclose(v2, v2b[k], rtol=0, atol=1e-12)


def test_hohmann_like_transfer():
    """180 deg-ish transfer between circular coplanar orbits: classic Hohmann numbers."""
    r_a, r_b = AU, 1.524 * AU
    a_t = 0.5 * (r_a + r_b)
    tof = np.pi * np.sqrt(a_t**3 / MU)
    ang = np.radians(179.999)  # exactly 180 deg is degenerate (plane undefined)
    r1 = np.array([r_a, 0.0, 0.0])
    r2 = r_b * np.array([np.cos(ang), np.sin(ang), 0.0])
    v1, _ = lambert(r1, r2, tof, MU)
    v_peri = np.sqrt(MU * (2 / r_a - 1 / a_t))
    assert np.linalg.norm(v1) == pytest.approx(v_peri, rel=1e-3)
