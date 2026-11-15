# Interplanetary Transfer: Earth → Itokawa

Mission design to the near-Earth asteroid 25143 Itokawa, from my BSc thesis in
Aerospace Engineering (Politecnico di Milano, 2025), and its extension in Python with a
Lambert-based formulation of the interplanetary transfer and modern optimisation methods.

| Folder | Content |
| --- | --- |
| [`MatlabCode/`](MatlabCode) | Original thesis code (MATLAB) |
| [`PythonCode/`](PythonCode) | Python re-implementation and extension of the interplanetary phase |

## 1. The thesis project (MATLAB)

The mission is split into three phases, each a separate scenario:

| Scenario | Phase | Reference frame | Method |
| --- | --- | --- | --- |
| 1 | Transfer from the initial Earth orbit to a parking orbit | Geocentric | Comparison of impulsive strategies (standard and alternative, e.g. bi-elliptic) |
| 2 | Interplanetary transfer Earth → Itokawa | Heliocentric | Grid search over the transfer geometry |
| 3 | Escape from the parking orbit and capture at Itokawa | Patched conics | Hyperbolic escape and capture manoeuvres |

Scenario 2 was solved with a brute-force search: three nested loops over the Earth
departure position, the asteroid arrival position and the argument of pericentre of the
transfer conic. Two versions were developed:

* **Standard**: the energetically optimal conic between the two orbits, without imposing
  that Itokawa is at the arrival point at the right time: 4.3008 km/s, TOF 222.8 days
  (no rendezvous).
* **Advanced**: the rendezvous is imposed through a time-tolerance window on the arrival
  epoch: **5.7204 km/s**, departure 30 May 2027, TOF 393.4 days. Running time about 6 hours.

## 2. The extension (Python): from grid search to Lambert's problem

For given departure and arrival positions, the loop over the argument of pericentre
scans the one-parameter family of conics through the two points, looking for the one
whose time of flight matches the asteroid's. This is exactly **Lambert's problem**,
which can be solved directly: given the departure epoch t0 and the time of flight T,
Lambert's solver returns the unique (zero-revolution) conic that reaches the asteroid in
exactly T. The transfer becomes a smooth function of two variables,

    ΔV(t0, T) = |v1 − v_Earth(t0)| + |v_Itokawa(t0 + T) − v2|,   (v1, v2) = Lambert(r_Earth(t0), r_Itokawa(t0 + T), T)

with an exact rendezvous and a constraint on the minimum distance from the Sun (≥ 0.7 AU,
as in the thesis).

Lambert's problem is solved with Izzo's algorithm (2015), compiled with Numba. The
optimisation follows the standard two-step scheme: a global search over (t0, T) finds
the best basin, then a gradient-based method (SQP) refines it to the exact minimum.

### Validation against MATLAB (thesis Keplerian model, 2025–2028)

| Check | Python | MATLAB |
| --- | --- | --- |
| Itokawa elements from the JPL Horizons state | identical (rel. error 1e-16) | — |
| Lambert evaluated at the MATLAB optimum (same t0, T) | 5.7206 km/s | 5.7204 km/s (0.003 %) |
| Best transfer, 3rd Earth revolution (Lambert grid + SQP) | **5.5777 km/s**, dep. 19 May 2027, TOF 410.0 d | 5.7204 km/s |
| Best transfer, 1st / 2nd revolution | 6.94 / 7.18 km/s | 7.99 / 7.90 km/s |
| Faithful Python port of the three MATLAB loops | 5.6660 km/s | 5.7204 km/s |
| Cost of the search | 570 495 Lambert solutions in < 1 s | 4.5·10⁹ conics in ≈ 6 h |

The second row shows that both codes compute the same cost function. The third shows
that the Lambert formulation finds a transfer 0.14 km/s cheaper in the same model.

### Why Lambert is better, and not by chance

* **The rendezvous is solved, not sampled.** It is an equality constraint: in the
  three-dimensional grid its solutions lie on a surface, which grid points almost never
  hit. With a ±12 s window and a 10⁴ s step, only 350 of the 4.5·10⁹ MATLAB conics pass
  the test (about one in 13 million). Lambert eliminates the third variable and every
  evaluation is a valid rendezvous.
* **The search is denser where it matters.** 570 thousand valid transfers on a regular
  1-day grid in (t0, T), against 350 valid transfers scattered where the tolerance window
  happened to catch them.
* **The minimum is continuous.** SQP moves freely in (t0, T) and converges to the bottom
  of the basin, instead of being restricted to grid nodes.

The difference between the MATLAB result (5.720) and the Python port of the same loops
(5.666) is instead due to chance: with so few valid conics, tiny numerical differences in
the propagation change which ones fall inside the window.

### Notes on the original MATLAB code

Found while porting it (details in `PythonCode/src/itokawa_opt/optimizers/grid_matlab.py`):

1. The transfer plane normal is `r_A × r_E`, so only prograde arcs longer than 180° are
   generated. Here the optimum is a 272° arc, so the result is unaffected; Lambert explores
   both cases.
2. The stored anomalies are one time step behind the positions actually used.
3. The time indices wrap with `N·dt`, not exactly the orbital period: after the wraps the
   Earth and asteroid clocks differ by about 8.5·10³ s, much more than the ±12 s window.
4. `deltaV_map` stores the last feasible value, not the minimum.

## 3. Next steps

| Step | Content | Status |
| --- | --- | --- |
| Validation | Kepler and element conversions, Lambert solver, MATLAB port, Lambert grid + SQP | done |
| Global optimisation | Differential Evolution written from scratch, comparison with pygmo algorithms, statistical benchmark against the grid | in progress |
| GPU | Lambert solver and DE population on the GPU (Numba CUDA), speed-up study | planned |

## Running the Python code

```bash
cd PythonCode
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run_validation.py     # validation, ~10 s (add --matlab-grid for the MATLAB port, ~15 min)
python -m pytest                     # automated tests
```

Results and figures are written to `PythonCode/results/validation/`.

Code structure (`PythonCode/src/itokawa_opt/`):

| Module | Role |
| --- | --- |
| `astro/` | Two-body mechanics: Kepler's equation, element conversions, Lambert solver |
| `ephemeris/` | Positions of Earth and Itokawa at any time |
| `problems/` | The optimisation problem: ΔV and constraint as a function of (t0, T) |
| `optimizers/` | MATLAB grid port, Lambert grid, SQP refinement |
| `viz/` | Porkchop plot and orbit plot |

Conventions: units km, s, km/s, rad; t = 0 is the Earth perihelion of 4 January 2025;
the validation uses the thesis values (μ_Sun = 1.327·10¹¹ km³/s², fixed Keplerian elements).

## References

* D. Izzo, "Revisiting Lambert's problem", *Celestial Mechanics and Dynamical Astronomy* 121 (2015).
* W. Martens et al., "SALTO: an expert-informed global trajectory design and optimisation
  toolkit", ESA GNC-ICATT 2023.