"""Figures for the validation: porkchop plot and 3D transfer orbit."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from itokawa_opt.constants import AU, DAY
from itokawa_opt.ephemeris.thesis_kepler import KeplerEphemeris

EPOCH0 = np.datetime64("2025-01-04")


def _to_dates(t_seconds):
    return EPOCH0 + np.round(np.asarray(t_seconds, dtype=float)).astype("timedelta64[s]")


def porkchop_plot(
    pc,
    markers: dict[str, tuple[float, float]] | None = None,
    dv_max: float = 10.0,
    path: str | None = None,
):
    """Contour of total delta-V over departure date and time of flight.

    ``markers``: name -> (t0 [s], T [s]) points to overlay (e.g. MATLAB optimum).
    """
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    dates = _to_dates(pc.t0)
    dv = np.clip(pc.dv, None, dv_max)
    levels = np.linspace(np.nanmin(dv), dv_max, 25)
    cs = ax.contourf(dates, pc.tof / DAY, dv.T, levels=levels, cmap="viridis_r")
    fig.colorbar(cs, ax=ax, label=r"$\Delta V$ total [km/s]")
    for (name, (t0, tof)), mk in zip((markers or {}).items(), ["o", "s", "^", "D"], strict=False):
        ax.plot(_to_dates(t0), tof / DAY, mk, ms=9, mfc="none", mec="red", mew=2, label=name)
    ax.set_xlabel("Departure date")
    ax.set_ylabel("Time of flight [days]")
    ax.set_title("Earth -> Itokawa rendezvous, thesis Keplerian model (Lambert, 0 rev)")
    if markers:
        ax.legend(loc="upper right")
    if path:
        fig.savefig(path, dpi=150)
    return fig


def transfer_plot(
    earth: KeplerEphemeris, asteroid: KeplerEphemeris, r1, v1, tof: float, path: str | None = None
):
    """Orbits of Earth and Itokawa plus the transfer arc (ecliptic projection)."""
    from itokawa_opt.astro.elements import car2par, par2car
    from itokawa_opt.astro.kepler import propagate_true_anomaly

    fig, ax = plt.subplots(figsize=(6.5, 6.5), constrained_layout=True)
    for body, color in ((earth, "tab:blue"), (asteroid, "tab:gray")):
        th = np.linspace(0, 2 * np.pi, 400)
        el = body.elements
        r, _ = par2car(el.a, el.e, el.i, el.raan, el.argp, th, body.mu)
        ax.plot(r[:, 0] / AU, r[:, 1] / AU, color=color, label=f"{body.name} orbit")
    tr = car2par(r1, v1, earth.mu)
    t = np.linspace(0.0, tof, 400)
    th = propagate_true_anomaly(tr.a, tr.e, tr.theta, t, earth.mu)
    r, _ = par2car(tr.a, tr.e, tr.i, tr.raan, tr.argp, th, earth.mu)
    ax.plot(r[:, 0] / AU, r[:, 1] / AU, color="tab:red", lw=2, label="Transfer")
    ax.plot(r[0, 0] / AU, r[0, 1] / AU, "o", color="tab:green", label="Departure")
    ax.plot(r[-1, 0] / AU, r[-1, 1] / AU, "o", color="tab:red", label="Arrival")
    ax.plot(0, 0, "*", color="gold", ms=14, mec="k", label="Sun")
    ax.set_aspect("equal")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)
    if path:
        fig.savefig(path, dpi=150)
    return fig
