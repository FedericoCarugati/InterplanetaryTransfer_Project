"""Statistics of the benchmark: checked against values computed by hand."""

import numpy as np
import pytest

from itokawa_opt.analysis.stats import (
    at_resolution,
    ecdf,
    evaluations_to_target,
    expected_running_time,
    fisher_success,
    holm_correction,
    mann_whitney,
    success_rate,
    vargha_delaney_a12,
    wilson_interval,
)


def test_wilson_interval_known_values():
    # 30/30 successes: the true rate is only known to be above ~0.886
    lo, hi = wilson_interval(30, 30)
    assert lo == pytest.approx(0.8865, abs=1e-3)
    assert hi == 1.0
    # 15/30: symmetric around 0.5
    lo, hi = wilson_interval(15, 30)
    assert 0.5 - lo == pytest.approx(hi - 0.5, abs=1e-12)
    assert lo == pytest.approx(0.3315, abs=1e-3)


def test_success_rate():
    rate, lo, hi = success_rate([0.0, 1e-4, 0.5, np.inf], tol=1e-3)
    assert rate == 0.5
    assert lo < 0.5 < hi


def test_evaluations_to_target():
    n, f = [20, 40, 60, 80], [9.0, 5.0, 3.0, 1.0]
    assert evaluations_to_target(n, f, 3.0) == 60
    assert evaluations_to_target(n, f, 0.5) == np.inf


def test_ert():
    # two successes after 100 and 300 evaluations, one failure after a budget of 1000
    assert expected_running_time([100, 300, np.inf], [100, 300, 1000]) == pytest.approx(700)
    assert expected_running_time([np.inf], [1000]) == np.inf


def test_ecdf_saturates_at_success_fraction():
    x, y = ecdf([3.0, 1.0, np.inf, 2.0])
    assert x.tolist() == [1.0, 2.0, 3.0]
    assert y.tolist() == [0.25, 0.5, 0.75]


def test_a12():
    assert vargha_delaney_a12([1, 2, 3], [4, 5, 6]) == 1.0
    assert vargha_delaney_a12([4, 5, 6], [1, 2, 3]) == 0.0
    assert vargha_delaney_a12([1, 2], [1, 2]) == 0.5


def test_mann_whitney_detects_shift_and_handles_ties_and_inf():
    rng = np.random.default_rng(0)
    a, b = rng.normal(0, 1, 30), rng.normal(2, 1, 30)
    c = mann_whitney(a, b, "a", "b")
    assert c.p_value < 1e-6
    assert c.verdict == "a is better"
    same = mann_whitney(np.zeros(10), np.zeros(10))
    assert same.p_value == 1.0
    with_inf = mann_whitney([0.0] * 10, [np.inf] * 10)
    assert with_inf.p_value < 1e-3 and with_inf.a12 == 1.0


def test_resolution_and_fisher():
    e = at_resolution([-1e-12, 1e-8, 2e-3], 1e-6)
    assert e.tolist() == [0.0, 0.0, 2e-3]
    fi = fisher_success([0.0] * 10, [1.0] * 10, tol=1e-3)
    assert (fi.successes_a, fi.successes_b) == (10, 0)
    assert fi.p_value < 1e-4


def test_holm():
    adj = holm_correction([0.01, 0.04, 0.03])
    assert adj == pytest.approx([0.03, 0.06, 0.06])