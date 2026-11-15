import numpy as np
import pytest

from itokawa_opt.astro.elements import car2par, par2car
from itokawa_opt.astro.kepler import (
    eccentric_to_true,
    period,
    propagate_true_anomaly,
    solve_kepler,
    time_of_flight,
    true_to_eccentric,
)

MU = 1.327e11


@pytest.mark.parametrize("e", [0.0, 0.016, 0.28, 0.7, 0.95])
def test_anomaly_round_trip(e):
    th = np.linspace(0, 2 * np.pi, 361, endpoint=False)
    back = eccentric_to_true(true_to_eccentric(th, e), e)
    assert np.allclose(np.mod(back - th + np.pi, 2 * np.pi) - np.pi, 0.0, atol=1e-12)


@pytest.mark.parametrize("e", [0.0, 0.3, 0.9, 0.99])
def test_solve_kepler_residual(e):
    m = np.linspace(0, 2 * np.pi, 1000)
    ecc = solve_kepler(m, e)
    res = np.mod(ecc - e * np.sin(ecc) - m + np.pi, 2 * np.pi) - np.pi
    assert np.max(np.abs(res)) < 1e-12


def test_tof_full_period_and_wrap():
    a, e = 1.97e8, 0.28
    t = period(a, MU)
    # half orbit pericentre -> apocentre is T/2
    assert time_of_flight(a, e, 0.0, np.pi, MU) == pytest.approx(t / 2, rel=1e-12)
    # going "backwards" in anomaly wraps through pericentre
    assert time_of_flight(a, e, 3.0, 1.0, MU) == pytest.approx(
        t - time_of_flight(a, e, 1.0, 3.0, MU), rel=1e-12
    )


def test_propagation_consistent_with_tof():
    a, e = 1.97e8, 0.28
    th0, dt = 1.2, 1.7e7
    th1 = propagate_true_anomaly(a, e, th0, dt, MU)
    assert time_of_flight(a, e, th0, th1, MU) == pytest.approx(dt, rel=1e-10)


def test_par2car_car2par_round_trip():
    el = (1.966e8, 0.2838, 0.0283, 1.204, 2.851, 3.275)
    r, v = par2car(*el, mu=MU)
    back = car2par(r, v, MU)
    assert np.allclose(back, el, rtol=1e-10, atol=1e-12)


def test_energy_and_momentum_conserved():
    el = (1.966e8, 0.2838, 0.0283, 1.204, 2.851, 0.0)
    th = np.linspace(0, 2 * np.pi, 50)
    r, v = par2car(*el[:5], th, MU)
    energy = 0.5 * np.sum(v**2, axis=1) - MU / np.linalg.norm(r, axis=1)
    h = np.linalg.norm(np.cross(r, v), axis=1)
    assert np.ptp(energy) / abs(energy.mean()) < 1e-12
    assert np.ptp(h) / h.mean() < 1e-12
