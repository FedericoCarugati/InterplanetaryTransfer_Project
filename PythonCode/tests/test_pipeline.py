"""Hybrid DE + SQP on the thesis transfer problem and the benchmark machinery."""

import json
from pathlib import Path

import numpy as np
import pytest

from itokawa_opt.analysis.benchmark import (
    load_records,
    method_de,
    method_grid,
    method_hybrid,
    run_repeated,
    save_records,
)
from itokawa_opt.constants import DAY
from itokawa_opt.ephemeris.thesis_kepler import thesis_bodies
from itokawa_opt.optimizers.de import DEConfig
from itokawa_opt.optimizers.grid_matlab import ThesisScenario
from itokawa_opt.optimizers.local import refine
from itokawa_opt.optimizers.pipeline import hybrid_de_sqp, select_candidates, verify_transfer
from itokawa_opt.problems.direct_transfer import DirectTransfer

REF = json.loads(
    (
        Path(__file__).resolve().parents[1] / "data/reference/matlab_scenario2_advanced.json"
    ).read_text()
)
SC = ThesisScenario()
BEST_KNOWN = 5.5776549779  # 2025-2028 window, from run_validation / run_benchmark


@pytest.fixture(scope="module")
def problem():
    earth, itokawa = thesis_bodies(SC)
    return DirectTransfer(
        earth, itokawa, (0.0, 3 * earth.period), (30 * DAY, itokawa.period), r_min_sun=SC.r_min_sun
    )


def test_fitness_batch_equals_scalar_fitness(problem):
    lb, ub = (np.asarray(b) for b in problem.get_bounds())
    x = lb + np.random.default_rng(0).random((50, 2)) * (ub - lb)
    f, g = problem.fitness_batch(x)
    for k in range(0, 50, 7):
        assert [f[k], g[k]] == pytest.approx(problem.fitness(x[k]), rel=1e-14)


def test_refine_counts_distinct_evaluations(problem):
    calls = []
    original = problem.fitness

    class Spy:
        def __getattr__(self, name):
            return getattr(problem, name)

        def fitness(self, x):
            calls.append(tuple(x))
            return original(x)

    loc = refine(Spy(), [865 * DAY, 405 * DAY])
    assert loc.nfev == len(calls) == len(set(calls))


def test_hybrid_finds_global_minimum_and_beats_matlab(problem):
    res = hybrid_de_sqp(problem, DEConfig(max_nfev=1500), seed=0)
    assert res.feasible
    assert res.dv == pytest.approx(BEST_KNOWN, abs=1e-6)
    assert res.dv < REF["optimum"]["dv_total"]
    assert res.nfev == res.nfev_global + res.nfev_local
    assert res.history_f[-1] == pytest.approx(res.dv)


def test_independent_check_of_the_solution(problem):
    res = hybrid_de_sqp(problem, DEConfig(max_nfev=1500), seed=1)
    chk = verify_transfer(problem, res.x)
    assert chk["position_miss_km"] < 1.0  # vs ~1.5e8 km heliocentric distance
    assert chk["dv_propagated"] == pytest.approx(chk["dv_lambert"], rel=1e-9)
    assert chk["constraint"] <= 0.0


def test_candidates_are_distinct_and_best_first():
    pop = np.array([[0.0, 0.0], [0.001, 0.0], [1.0, 1.0], [0.5, 0.5]])
    f = np.array([1.0, 0.5, 2.0, 9.0])
    v = np.zeros(4)
    c = select_candidates(pop, f, v, ([0, 0], [1, 1]), n_candidates=3, min_distance=0.02)
    assert c.tolist() == [[0.001, 0.0], [1.0, 1.0], [0.5, 0.5]]


def test_benchmark_records_round_trip(problem, tmp_path):
    recs = run_repeated(method_de(DEConfig()), problem, 400, n_runs=3, master_seed=1)
    recs += run_repeated(method_hybrid(DEConfig()), problem, 400, n_runs=2, master_seed=1)
    recs.append(method_grid(40.0, polish=True)(problem, 0, None))
    assert len({r.seed for r in recs[:3]}) == 3  # independent spawned seeds
    assert all(r.nfev <= 400 for r in recs[:3])
    assert recs[0].best_at(20) >= recs[0].best_at(400)
    path = tmp_path / "records.json"
    save_records(recs, path)
    back = load_records(path)
    assert [r.f for r in back] == [r.f for r in recs]