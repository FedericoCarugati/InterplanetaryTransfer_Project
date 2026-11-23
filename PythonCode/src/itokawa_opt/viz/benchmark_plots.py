"""Figures of the optimiser benchmark (Level 1).

Errors are shown in m/s, e = (dV - dV_ref) * 1000, on a log scale; errors below
``FLOOR`` (numerical precision) are drawn at the floor.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from itokawa_opt.constants import DAY
from itokawa_opt.viz.plots import _to_dates

FLOOR = 1e-9  # m/s
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def _err_ms(f, f_ref):
    e = (np.asarray(f, dtype=float) - f_ref) * 1e3
    return np.clip(np.where(np.isfinite(e), e, np.nan), FLOOR, None)


def _porkchop(ax, pc, dv_max=10.0):
    dv = np.clip(pc.dv, None, dv_max)
    levels = np.linspace(np.nanmin(dv), dv_max, 25)
    cs = ax.contourf(_to_dates(pc.t0), pc.tof / DAY, dv.T, levels=levels, cmap="viridis_r")
    ax.set_xlabel("Departure date")
    ax.set_ylabel("Time of flight [days]")
    return cs


def _save(fig, path):
    if path:
        fig.savefig(path, dpi=150)
    return fig


def landscape_plot(pc, ref, path=None):
    """Porkchop of the reference grid with the local minima (basins) numbered."""
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    cs = _porkchop(ax, pc)
    fig.colorbar(cs, ax=ax, label=r"$\Delta V$ total [km/s]")
    for k, b in enumerate(ref.basins, start=1):
        ax.plot(_to_dates(b.x[0]), b.x[1] / DAY, "o", ms=9, mfc="none", mec="red", mew=2)
        ax.annotate(
            f"{k}: {b.f:.3f}",
            (_to_dates(b.x[0]), b.x[1] / DAY),
            xytext=(6, 6),
            textcoords="offset points",
            color="red",
            fontsize=9,
            fontweight="bold",
        )
    ax.set_title(
        f"Cost landscape: {len(ref.basins)} local minima (numbered by $\\Delta V$ [km/s]); "
        f"{100 * ref.infeasible_fraction:.0f} % of the box infeasible (white)"
    )
    return _save(fig, path)


def population_plot(pc, snapshots: dict, path=None):
    """DE population at a few generations, over the porkchop (exploration -> exploitation)."""
    gens = sorted(snapshots)
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True, sharex=True)
    for ax, g in zip(axes.flat, gens, strict=False):
        cs = _porkchop(ax, pc)
        pop = snapshots[g]
        ax.plot(_to_dates(pop[:, 0]), pop[:, 1] / DAY, "o", ms=5, mfc="white", mec="red")
        ax.set_title(f"Generation {g}")
    fig.colorbar(cs, ax=axes, label=r"$\Delta V$ total [km/s]", shrink=0.8)
    return _save(fig, path)


def convergence_plot(anytime: dict, fixed: dict, f_ref: float, tol: float, path=None):
    """Median (line) and interquartile range (band) of the error vs evaluations.

    ``anytime``: name -> list of RunRecord with a full history (DE, random, multi-start).
    ``fixed``: name -> list of (nfev, error array or scalar) points (grid, hybrid runs at
    several budgets).
    """
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    for name, recs in anytime.items():
        n_max = max(r.nfev for r in recs)
        grid = np.unique(np.geomspace(min(r.history_nfev[0] for r in recs), n_max, 200).astype(int))
        curves = np.array([[r.best_at(b) for b in grid] for r in recs])
        e = _err_ms(curves, f_ref)
        med = np.nanmedian(e, axis=0)
        q1, q3 = np.nanpercentile(e, 25, axis=0), np.nanpercentile(e, 75, axis=0)
        (line,) = ax.plot(grid, med, lw=2, label=f"{name} (median, IQR)")
        ax.fill_between(grid, q1, q3, color=line.get_color(), alpha=0.2)
    for (name, pts), mk in zip(fixed.items(), MARKERS, strict=False):
        nfev = np.array([p[0] for p in pts], dtype=float)
        vals = [np.atleast_1d(p[1]) for p in pts]
        med = np.array([np.nanmedian(_err_ms(v, f_ref)) for v in vals])
        lo = np.array([np.nanpercentile(_err_ms(v, f_ref), 25) for v in vals])
        hi = np.array([np.nanpercentile(_err_ms(v, f_ref), 75) for v in vals])
        ax.errorbar(
            nfev, med, yerr=[med - lo, hi - med], fmt=mk + "--", ms=7, capsize=3, label=name
        )
    ax.axhline(tol * 1e3, color="k", ls=":", lw=1, label=f"success threshold ({tol * 1e3:g} m/s)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Objective evaluations (Lambert solutions)")
    ax.set_ylabel(r"$\Delta V - \Delta V_{ref}$ [m/s]")
    ax.set_title("Convergence: best value found vs cost")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower left")
    return _save(fig, path)


def success_plot(table: dict, path=None):
    """Success rate vs evaluations with 95 % Wilson intervals.

    ``table``: name -> list of (nfev, rate, low, high).
    """
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for (name, rows), mk in zip(table.items(), MARKERS, strict=False):
        rows = sorted(rows)
        n = np.array([r[0] for r in rows], dtype=float)
        p, lo, hi = (np.array([r[k] for r in rows]) for k in (1, 2, 3))
        ax.errorbar(n, p, yerr=[p - lo, hi - p], fmt=mk + "-", ms=6, capsize=3, label=name)
    ax.set_xscale("log")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Objective evaluations (Lambert solutions)")
    ax.set_ylabel("Success rate")
    ax.set_title("Probability of reaching the global minimum (95 % Wilson intervals)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower right")
    return _save(fig, path)


def ecdf_plot(evals: dict, budget: float, deterministic: dict | None = None, path=None):
    """ECDF of evaluations-to-target: fraction of runs solved within n evaluations.

    ``evals``: name -> array of evaluations-to-target (inf = never). ``deterministic``:
    name -> nfev of a deterministic method that reaches the target (vertical line).
    """
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for name, e in evals.items():
        e = np.sort(np.asarray(e, dtype=float))
        y = np.arange(1, e.size + 1) / e.size
        fin = np.isfinite(e)
        x = np.r_[e[fin], budget]
        yy = np.r_[y[fin], y[fin][-1] if fin.any() else 0.0]
        ax.step(
            np.r_[x[0] if fin.any() else budget, x], np.r_[0.0, yy], where="post", lw=2, label=name
        )
    for name, n in (deterministic or {}).items():
        ax.axvline(n, ls="--", lw=1.2, color="0.3")
        ax.text(n, 0.02, " " + name, rotation=90, va="bottom", fontsize=8, color="0.3")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Objective evaluations")
    ax.set_ylabel("Fraction of runs that reached the target")
    ax.set_title("Fixed-target view: ECDF of the evaluations needed")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="upper left")
    return _save(fig, path)


def final_error_boxplot(errors: dict, tol: float, budget: int, path=None):
    """Distribution of the final errors at the same budget (box = IQR, dots = runs)."""
    names = list(errors)
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    data = [_err_ms(errors[n], 0.0) for n in names]
    ax.boxplot([d[~np.isnan(d)] for d in data], whis=(0, 100), widths=0.5)
    rng = np.random.default_rng(0)
    for k, d in enumerate(data, start=1):
        ax.plot(k + rng.uniform(-0.12, 0.12, d.size), d, ".", color="tab:red", alpha=0.6)
    ax.axhline(tol * 1e3, color="k", ls=":", lw=1, label=f"success threshold ({tol * 1e3:g} m/s)")
    ax.set_xticks(range(1, len(names) + 1), names, rotation=20, ha="right", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\Delta V - \Delta V_{ref}$ [m/s]")
    ax.set_title(f"Final error of every run at a budget of ~{budget} evaluations")
    ax.grid(alpha=0.3, axis="y", which="both")
    ax.legend(fontsize=8)
    return _save(fig, path)


def final_points_plot(pc, records: dict, path=None):
    """Where the runs of each method ended, over the porkchop (which basins they found)."""
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    cs = _porkchop(ax, pc)
    fig.colorbar(cs, ax=ax, label=r"$\Delta V$ total [km/s]")
    for (name, recs), mk in zip(records.items(), MARKERS, strict=False):
        x = np.array([r.x for r in recs])
        ax.plot(_to_dates(x[:, 0]), x[:, 1] / DAY, mk, ms=8, mfc="none", mew=1.5, label=name)
    ax.set_title("End points of the runs")
    ax.legend(fontsize=8, loc="upper right")
    return _save(fig, path)


def sensitivity_heatmap(F_values, CR_values, rates, budget, path=None):
    """Success rate of DE for each (F, CR) pair."""
    fig, ax = plt.subplots(figsize=(6.5, 5), constrained_layout=True)
    im = ax.imshow(rates, origin="lower", cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(CR_values)), [f"{c:g}" for c in CR_values])
    ax.set_yticks(range(len(F_values)), [f"{f:g}" for f in F_values])
    for i in range(len(F_values)):
        for j in range(len(CR_values)):
            ax.text(j, i, f"{rates[i, j]:.2f}", ha="center", va="center", fontsize=9)
    ax.set_xlabel("CR (crossover probability)")
    ax.set_ylabel("F (scale factor)")
    ax.set_title(f"DE success rate at {budget} evaluations")
    fig.colorbar(im, ax=ax, label="Success rate")
    return _save(fig, path)