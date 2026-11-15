"""Time of flight of the grid-search conics as a function of the argument of pericentre.

For the departure and arrival points of the selected grid-search transfer, the MATLAB loop over
omega_t generates 100 conics through the two points. This script plots their time of flight
against the asteroid travel time (the rendezvous condition), and marks the exact solution
found by Lambert's solver.

Run from the PythonCode folder:

    python scripts/plot_tof_vs_omega.py

Output: results/validation/tof_vs_omega.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from itokawa_opt.astro.elements import car2par
from itokawa_opt.astro.kepler import TWO_PI, time_of_flight, wrap_2pi
from itokawa_opt.astro.lambert import lambert
from itokawa_opt.constants import DAY
from itokawa_opt.ephemeris.thesis_kepler import thesis_bodies
from itokawa_opt.optimizers.grid_matlab import ThesisScenario

ROOT = Path(__file__).resolve().parents[1]
REF = json.loads((ROOT / "data/reference/matlab_scenario2_advanced.json").read_text())
OUT = ROOT / "results/validation"


def conic_tof(omega, r1v, r2v, mu):
    """Time of flight of the conics through r1, r2 with argument of pericentre omega.

    Same geometry as the MATLAB grid search: plane normal along r2 x r1, conic through the
    two points (e and a from the polar equation). Returns NaN for non-elliptic conics.
    """
    h = np.cross(r2v, r1v)
    h /= np.linalg.norm(h)
    node = np.array([-h[1], h[0], 0.0])
    node /= np.linalg.norm(node)
    m_hat = np.cross(h, node)
    u1 = np.arctan2(m_hat @ r1v, node @ r1v)  # argument of latitude of the two points
    u2 = np.arctan2(m_hat @ r2v, node @ r2v)
    r1, r2 = np.linalg.norm(r1v), np.linalg.norm(r2v)
    th1 = wrap_2pi(u1 - omega)
    th2 = wrap_2pi(u2 - omega)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = (r1 - r2) / (r2 * np.cos(th2) - r1 * np.cos(th1))
        a = r1 * (1.0 + e * np.cos(th1)) / (1.0 - e**2)
    ok = (e >= 0.0) & (e < 1.0)
    tof = time_of_flight(np.where(ok, a, 1.0), np.where(ok, e, 0.0), th1, th2, mu)
    return np.where(ok, tof, np.nan)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sc = ThesisScenario()
    earth, itokawa = thesis_bodies(sc)

    # Departure epoch and time of flight of the MATLAB optimum (same reconstruction as
    # run_validation.py: anomaly label + one time step, third Earth revolution)
    opt = REF["optimum"]
    t_label = float(
        time_of_flight(earth.elements.a, earth.elements.e, 0.0, opt["theta_earth_label"], sc.mu)
    )
    t0 = t_label + sc.dt + 2 * earth.period
    tof_target = opt["tof"]  # asteroid travel time [s]

    r1v, _ = earth.rv(t0)
    r2v, _ = itokawa.rv(t0 + tof_target)

    # The 100 values of the MATLAB loop and a fine curve
    omega_grid = np.linspace(0.0, TWO_PI, sc.n_omega)
    tof_grid = conic_tof(omega_grid, r1v, r2v, sc.mu)
    omega_fine = np.linspace(0.0, TWO_PI, 20001)
    tof_fine = conic_tof(omega_fine, r1v, r2v, sc.mu)

    # Exact rendezvous: root of TOF(omega) - T, bracketed by the two grid samples around it
    valid = np.isfinite(tof_grid)
    idx = np.where(
        valid[:-1]
        & valid[1:]
        & (np.sign(tof_grid[:-1] - tof_target) != np.sign(tof_grid[1:] - tof_target))
    )[0]
    k = int(idx[0])
    omega_star = brentq(
        lambda w: float(conic_tof(np.array(w), r1v, r2v, sc.mu)) - tof_target,
        omega_grid[k],
        omega_grid[k + 1],
        xtol=1e-14,
    )

    # Cross-check with Lambert: argument of pericentre of the Lambert transfer orbit
    v1, _ = lambert(r1v, r2v, tof_target, sc.mu)
    omega_lambert = car2par(r1v, v1, sc.mu).argp

    # Numbers quoted in the report
    slope = np.gradient(tof_fine, omega_fine)[np.nanargmin(np.abs(tof_fine - tof_target))]
    d_omega = TWO_PI / (sc.n_omega - 1)
    print(f"elliptic conics in the loop: {int(valid.sum())} of {sc.n_omega}")
    print(f"omega* (root of TOF = T): {omega_star:.6f} rad")
    print(f"omega of the Lambert transfer orbit: {omega_lambert:.6f} rad")
    err = np.abs(tof_grid[[k, k + 1]] - tof_target)
    print(
        f"samples around the root: omega = {omega_grid[k]:.4f}, {omega_grid[k + 1]:.4f} rad; "
        f"TOF errors = {err[0]:.0f} s, {err[1]:.0f} s"
    )
    print(
        f"|dTOF/domega| = {abs(slope):.2e} s/rad; TOF spacing between samples = "
        f"{abs(slope) * d_omega / DAY:.1f} days; hit probability for +-12 s = "
        f"{24.0 / (abs(slope) * d_omega):.1e}"
    )

    # Figure
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    ax.plot(omega_fine, tof_fine / DAY, color="0.6", lw=1, label="TOF of the conic family")
    ax.plot(
        omega_grid,
        tof_grid / DAY,
        "o",
        ms=4,
        color="tab:blue",
        label=rf"Grid samples ($n = {sc.n_omega}$, {int(valid.sum())} elliptic)",
    )
    ax.axhline(
        tof_target / DAY,
        color="tab:red",
        ls="--",
        lw=1.2,
        label=f"Asteroid travel time ({tof_target / DAY:.1f} days)",
    )
    ax.plot(
        omega_star,
        tof_target / DAY,
        "*",
        ms=14,
        color="tab:red",
        mec="k",
        label=rf"Exact rendezvous (Lambert), $\omega_t^* = {omega_star:.3f}$ rad",
    )
    ax.set_yscale("log")
    ax.set_ylim(50, 3e4)
    ax.set_xlabel(r"Argument of pericentre of the transfer orbit $\omega_t$ [rad]")
    ax.set_ylabel("Time of flight [days]")
    ax.set_xlim(
        omega_fine[np.isfinite(tof_fine)].min() - 0.1, omega_fine[np.isfinite(tof_fine)].max() + 0.1
    )
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=9)

    # Zoom on the intersection
    ins = ax.inset_axes([0.17, 0.47, 0.34, 0.40])
    win = (omega_fine > omega_grid[k] - 0.05) & (omega_fine < omega_grid[k + 1] + 0.05)
    ins.plot(omega_fine[win], tof_fine[win] / DAY, color="0.6", lw=1)
    sel = (omega_grid > omega_grid[k] - 0.05) & (omega_grid < omega_grid[k + 1] + 0.05)
    ins.axhline(tof_target / DAY, color="tab:red", ls="--", lw=1)
    ins.plot(omega_star, tof_target / DAY, "*", ms=16, color="tab:red", mec="k", zorder=3)
    ins.plot(
        omega_grid[sel], tof_grid[sel] / DAY, "o", ms=5, color="tab:blue", mec="white", zorder=4
    )
    ins.set_title("Zoom: two consecutive samples are 16 days apart in TOF", fontsize=8)
    ins.tick_params(labelsize=7)
    ins.grid(alpha=0.3)

    path = OUT / "tof_vs_omega.png"
    fig.savefig(path, dpi=200)
    print(f"Figure written to {path}")


if __name__ == "__main__":
    main()
