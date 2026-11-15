# itokawa-opt

Global and local optimisation of a direct Earth → Itokawa rendezvous transfer,
extending the interplanetary phase of my BSc thesis (Politecnico di Milano, 2025).

The thesis solved the transfer in MATLAB with a brute-force grid over three geometric
parameters (Earth departure position, asteroid arrival position, argument of
pericentre of the transfer conic), imposing the rendezvous through a time-tolerance
window. This package reformulates the problem as a **Lambert problem in (t0, T)**,
validates it against the MATLAB results, and then (Levels 1–2) adds global
optimisers, gradient-based refinement and GPU acceleration.

| Level | Content | Status |
| --- | --- | --- |
| 0 – Validation | Kepler/elements, Izzo's Lambert solver, faithful port of the MATLAB grid, Lambert porkchop + SQP | done |
| 1 – Global + local | Differential Evolution from scratch, pygmo algorithms, statistical benchmark vs grid | planned |
| 2 – GPU | Lambert and DE population on the GPU (Numba CUDA, CuPy), speed-up study | planned |

## Level 0 – Validation results (thesis Keplerian model, 2025–2028)

`python scripts/run_validation.py` (about 10 s):

| Check | Python | MATLAB | Difference |
| --- | --- | --- | --- |
| Itokawa elements from the Horizons state (`car2par`) | a = 1.966476e8 km, e = 0.2838488 | same | 1e-16 (rel.) |
| Lambert evaluated at the MATLAB optimum (same t0, T) | 5.7206 km/s | 5.7204 km/s | 0.003 % |
| Best transfer, 3rd Earth revolution (Lambert grid 1 day + SQP) | **5.5777 km/s**, dep. 2027-05-19, TOF 410.0 d | 5.7204 km/s, dep. ≈ 2027-06-01, TOF 393.4 d | −0.143 km/s (−2.5 %) |
| Best transfer, 1st revolution | 6.9396 km/s | 7.9893 km/s | −1.05 km/s |
| Best transfer, 2nd revolution | 7.1840 km/s | 7.9034 km/s | −0.72 km/s |
| Faithful Python port of the MATLAB loops (`--matlab-grid`) | 5.6660 km/s, dep. 2027-05-08, TOF 424.9 d; 350 feasible conics out of 4.5·10⁹ | 5.7204 km/s | see note below |
| Cost of the search | Lambert grid: 570 495 solutions in < 1 s; port of the loops: ≈ 15 min (NumPy, 1 core) | ≈ 4.5·10⁹ conics, ≈ 6 h | |

Why Lambert finds lower values with far fewer evaluations:

* For given departure and arrival positions, the MATLAB loop over the argument of
  pericentre (100 values) is a coarse, sampled way of finding the conic with the
  right time of flight. Lambert's solver finds that conic exactly.
* With a ±12 s tolerance and a 1e4 s time step, only 350 of the 4.5·10⁹ grid
  points satisfy the rendezvous condition, so the grid optimum depends on which
  samples happen to fall inside the window. The Python port uses the same loops
  but propagates the anomaly analytically from t = 0, while MATLAB chains
  `KepleroDiretto` step by step (the most likely source of the small differences,
  of order 1e-4 rad, observed in the anomalies): that is enough to change which
  conics pass the window. That is why the port finds 5.666 km/s and MATLAB
  5.720 km/s, while Lambert evaluated at the MATLAB (t0, T) reproduces 5.7204 km/s.
* SQP then moves continuously in (t0, T) instead of being restricted to grid nodes.

Notes on the MATLAB code found during the port (documented in
`optimizers/grid_matlab.py`):

1. The transfer plane normal is `r_A × r_E`: only prograde arcs with a transfer
   angle above 180° are generated. In this window the optimum is indeed a long-way
   arc (272°), so this did not change the result, but Lambert explores both cases.
2. The state arrays are filled with the state at `theta(j+1)` while the anomaly
   array uses `theta(j)`: the stored anomalies are one time step behind the
   positions actually used.
3. Earth and asteroid indices wrap with `N·dt`, which is not exactly the orbital
   period. After the wraps the Earth departure time and the asteroid "clock" differ
   by about 8.5·10³ s at the optimum, much more than the ±12 s rendezvous window.
4. `deltaV_map` stores the last feasible value over `s`, not the minimum.

## Installation

```bash
cd PythonCode
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install
```

## Usage

```bash
python scripts/run_validation.py                 # validation, ~10 s
python scripts/run_validation.py --matlab-grid   # + faithful port of the MATLAB loops
pytest                                           # test suite
```

Figures and `summary.json` are written to `results/validation/`.

## Structure

```text
src/itokawa_opt/
  astro/        kepler.py, elements.py (par2car, car2par), lambert.py (Izzo 2015, Numba)
  ephemeris/    thesis_kepler.py (thesis model), spice.py (Level 1)
  problems/     direct_transfer.py (objective + Sun-distance constraint, batch evaluation)
  optimizers/   grid_matlab.py (MATLAB port), grid_lambert.py (porkchop), local.py (SQP),
                de.py / pygmo_adapter.py / pipeline.py (Level 1)
  gpu/          Level 2
  analysis/     benchmarks and statistics (Level 1)
  viz/          porkchop and orbit plots
tests/          pytest suite (Lambert checked against ESA pykep to 1e-10)
data/reference/ MATLAB reference values
```

## Conventions

Units: km, s, km/s, rad. `t = 0` is the Earth perihelion of 4 January 2025, the
reference instant of the thesis; Itokawa's state at that instant comes from JPL
Horizons. The Level 0 model uses the thesis values (μ_Sun = 1.327e11 km³/s²,
fixed osculating elements); Level 1 switches to SPICE ephemerides.

## References

* D. Izzo, "Revisiting Lambert's problem", *Celest. Mech. Dyn. Astr.* 121 (2015).
* W. Martens et al., "SALTO: an expert-informed global trajectory design and
  optimisation toolkit", ESA GNC-ICATT 2023.
