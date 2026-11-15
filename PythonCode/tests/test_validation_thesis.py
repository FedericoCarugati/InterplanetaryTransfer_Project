"""Regression tests against the MATLAB results of the thesis (Scenario 2, advanced)."""

import json
from pathlib import Path

import numpy as np
import pytest

from itokawa_opt.astro.kepler import time_of_flight
from itokawa_opt.constants import DAY
from itokawa_opt.ephemeris.thesis_kepler import thesis_bodies
from itokawa_opt.optimizers.grid_lambert import porkchop
from itokawa_opt.optimizers.grid_matlab import ThesisScenario, run_grid_search
from itokawa_opt.optimizers.local import refine
from itokawa_opt.problems.direct_transfer import DirectTransfer

REF = json.loads(
    (
        Path(__file__).resolve().parents[1] / "data/reference/matlab_scenario2_advanced.json"
    ).read_text()
)
SC = ThesisScenario()


@pytest.fixture(scope="module")
def bodies():
    return thesis_bodies(SC)


@pytest.fixture(scope="module")
def problem(bodies):
    earth, itokawa = bodies
    return DirectTransfer(
        earth, itokawa, (0.0, 3 * earth.period), (30 * DAY, itokawa.period), r_min_sun=SC.r_min_sun
    )


def test_itokawa_elements_match_matlab(bodies):
    el = bodies[1].elements
    ref = REF["asteroid_from_horizons"]
    assert el.a == pytest.approx(ref["a"], rel=1e-12)
    assert el.e == pytest.approx(ref["e"], rel=1e-12)
    assert el.theta == pytest.approx(ref["theta0"], abs=1e-9)


def test_lambert_at_matlab_optimum(bodies, problem):
    """Same (t0, T) as the MATLAB optimum -> same delta-V within 0.05 %."""
    earth, _ = bodies
    opt = REF["optimum"]
    t_label = float(
        time_of_flight(earth.elements.a, earth.elements.e, 0.0, opt["theta_earth_label"], SC.mu)
    )
    t0 = t_label + SC.dt + 2 * earth.period
    dv = problem.evaluate(np.array([t0]), np.array([opt["tof"]]))["dv"][0]
    assert dv == pytest.approx(opt["dv_total"], rel=5e-4)


def test_lambert_beats_or_matches_grid(problem, bodies):
    """Continuous optimum in the 3rd revolution must not be worse than MATLAB's."""
    earth, itokawa = bodies
    t0 = np.arange(2 * earth.period, 3 * earth.period, 2 * DAY)
    tof = np.arange(30 * DAY, itokawa.period, 2 * DAY)
    pc = porkchop(problem, t0, tof)
    loc = refine(problem, [pc.best_t0, pc.best_tof])
    assert loc.success
    assert loc.constraint <= 1e-9
    assert loc.dv <= REF["optimum"]["dv_total"] + 1e-6


def test_matlab_port_smoke():
    """The faithful port runs and finds feasible conics on a short departure range."""
    res = run_grid_search(SC, k_range=(7500, 7700), verbose=False)
    assert res.n_feasible > 0
    assert 5.0 < res.dv_min < 9.0
