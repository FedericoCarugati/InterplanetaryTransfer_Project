"""Keplerian ephemerides of the thesis (fixed osculating elements, two-body motion).

``t`` is measured in seconds from ``t = 0`` = Earth perihelion of 4 Jan 2025, the
reference instant used in the MATLAB script (Itokawa's state at that instant was
taken from JPL Horizons).

Level 1 will add a SPICE-based ephemeris with the same ``rv(t)`` interface.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from itokawa_opt.astro.elements import Elements, car2par, par2car
from itokawa_opt.astro.kepler import period, propagate_true_anomaly
from itokawa_opt.optimizers.grid_matlab import ThesisScenario


@dataclass(frozen=True)
class KeplerEphemeris:
    """Body on a fixed Keplerian orbit; ``elements.theta`` is the anomaly at t = 0."""

    elements: Elements
    mu: float
    name: str = ""

    @property
    def period(self) -> float:
        return float(period(self.elements.a, self.mu))

    def true_anomaly(self, t: ArrayLike) -> NDArray[np.float64]:
        el = self.elements
        return propagate_true_anomaly(el.a, el.e, el.theta, t, self.mu)

    def rv(self, t: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Position and velocity at time(s) t, shape (..., 3) [km, km/s]."""
        el = self.elements
        return par2car(el.a, el.e, el.i, el.raan, el.argp, self.true_anomaly(t), self.mu)


def thesis_bodies(sc: ThesisScenario | None = None) -> tuple[KeplerEphemeris, KeplerEphemeris]:
    """Earth and Itokawa exactly as defined in the MATLAB script."""
    sc = sc or ThesisScenario()
    earth = KeplerEphemeris( Elements(sc.a_earth, sc.e_earth, sc.i_earth, sc.raan_earth, 
                                      sc.argp_earth, sc.theta0_earth), sc.mu, "Earth")
    itokawa = KeplerEphemeris(car2par(sc.r0_ast, sc.v0_ast, sc.mu), sc.mu, "Itokawa")
    return earth, itokawa
