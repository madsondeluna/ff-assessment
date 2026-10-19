#!/usr/bin/env python3
# =============================================================================
# 04_plot_comparison.py
# Comparative figures and summary statistics for the four-way comparison
# (CHARMM36m, AMBER, GROMOS 53A6, GROMOS 54A7), from the outputs of
# 03_analysis.sh (analysis/<ff>/rep<N>/). Any subset of force fields and any
# number of replicas is accepted; difference panels use one force field as
# reference (--reference, default charmm36m).
#
# Two levels of comparison, as in a multi-replica benchmark:
#   between replicas   fig7_replicas: each replica of each force field on
#                      its own, to judge convergence and run-to-run spread
#   between force fields   all other figures: replica means with SD bands,
#                      distributions pooled over replicas, and "ff minus
#                      reference" panels
#
# Figures (figures/):
#   fig1_timeseries          RMSD, Rg, H-bonds (intra, water), SASA: replica
#                            mean (running mean) with min-max band
#   fig2_distributions       distributions over the equilibrated window
#   fig3_rmsf                C-alpha RMSF per residue + difference vs reference
#   fig4_contact_maps        mean residue distance maps + differences vs reference
#   fig5_sasa_per_residue    heavy-atom SASA per residue + difference vs reference
#   fig6_quality_checks      EM, NVT, NPT and production thermodynamics
#   fig7_replicas            per-replica time series and RMSF for each force field
# Tables:
#   summary_statistics.csv   mean, SD and block SEM per observable, replica
#                            and force field (replica "all" = mean +/- SD of
#                            the replica means)
#   summary_table.csv        one row per force field: mean +/- SD over
#                            replicas for every observable (the table to
#                            report), plus the difference vs the reference
#
# Usage:
#   python3 scripts/04_plot_comparison.py
#   python3 scripts/04_plot_comparison.py --reference amber --smooth-ns 2
# Requires: numpy, matplotlib (scipy optional, for Welch t-test p-values)
# =============================================================================

import argparse
import csv
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless backend: works on servers without a display
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

try:
    from scipy import stats as _scipy_stats
except ImportError:  # p-values are then reported as nan
    _scipy_stats = None

# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------
# Fixed categorical order and colors (Okabe-Ito, colorblind-safe): a force
# field always keeps its color, whatever subset of data is plotted. Labels not
# listed here (e.g. a fifth force field) get the remaining palette entries.
FF_ORDER = ["charmm36m", "amber", "gromos53a6", "gromos54a7"]
FF_COLOR = {
    "charmm36m": "#0072b2",   # blue
    "amber": "#009e73",       # green
    "gromos53a6": "#e69f00",  # orange
    "gromos54a7": "#d55e00",  # vermillion
}
EXTRA_COLORS = ["#cc79a7", "#56b4e9", "#8c8c8c", "#000000"]
# Fallback names when meta.txt has no display_name (older analyses)
FF_NAME_DEFAULT = {"charmm36m": "CHARMM36m", "charmm36": "CHARMM36m", "amber": "AMBER",
                   "gromos53a6": "GROMOS 53A6", "gromos54a7": "GROMOS 54A7"}
FF_NAME = dict(FF_NAME_DEFAULT)  # filled from meta.txt in discover()

REP_STYLE = ["-", "--", ":", "-.", (0, (5, 1, 1, 1))]  # one line style per replica

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"

# Sequential (one hue, dark = close contact) and diverging (blue - gray - red)
SEQ_CMAP = LinearSegmentedColormap.from_list(
    "seq_blue", ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#fcfcfb"])
DIV_CMAP = LinearSegmentedColormap.from_list(
    "div_blue_red", ["#184f95", "#6da7ec", "#f0efec", "#ef8a7e", "#b8322f"])

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "axes.edgecolor": INK_2,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.5,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def color_of(ff):
    if ff not in FF_COLOR:
        FF_COLOR[ff] = EXTRA_COLORS[len(FF_COLOR) % len(EXTRA_COLORS)]
    return FF_COLOR[ff]


def name_of(ff):
    return FF_NAME.get(ff, ff)


# -----------------------------------------------------------------------------
# Readers
# -----------------------------------------------------------------------------
def read_xvg(path):
    """Read a GROMACS .xvg file.

    Returns (data, meta). data is a 2D array (rows = frames). meta holds the
    axis labels and legends. Time columns are converted to ns when the x-axis
    label says ps, so all time series share the same unit.
    """
    meta = {"xlabel": "", "ylabel": "", "legends": []}
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("@"):
                m = re.search(r'xaxis\s+label\s+"(.*)"', line)
                if m:
                    meta["xlabel"] = m.group(1)
                m = re.search(r'yaxis\s+label\s+"(.*)"', line)
                if m:
                    meta["ylabel"] = m.group(1)
                m = re.search(r's\d+\s+legend\s+"(.*)"', line)
                if m:
                    meta["legends"].append(m.group(1))
                continue
            parts = line.split()
            if parts:
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    continue  # tolerate stray text lines
    if not rows:
        raise ValueError(f"no numeric data in {path}")
    width = min(len(r) for r in rows)
    data = np.array([r[:width] for r in rows])
    if "(ps)" in meta["xlabel"]:
        data[:, 0] /= 1000.0
        meta["xlabel"] = meta["xlabel"].replace("(ps)", "(ns)")
    return data, meta


def read_xpm(path):
    """Read a GROMACS .xpm matrix (e.g. gmx mdmat -mean).

    Each color line carries its numeric value in a comment:
        "AB c #FFFFFF " /* "0.5" */,
    Symbols can be 1 or 2 characters and may include non-ASCII bytes, so the
    file is read as latin-1 and parsed line by line.
    Returns (matrix, x_axis, y_axis, title) with matrix[i, j] at y_axis[i], x_axis[j].
    """
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    text = "\n".join(lines)
    title = re.search(r'title:\s*"(.*)"', text)
    title = title.group(1) if title else ""

    quoted = [ln for ln in lines if ln.startswith('"')]
    ncol, nrow, ncolors, cpp = (int(v) for v in quoted[0].strip('",').split()[:4])

    symbol_value = {}
    for ln in quoted[1:1 + ncolors]:
        sym = ln[1:1 + cpp]
        val = re.search(r'/\*\s*"([^"]*)"\s*\*/', ln)
        symbol_value[sym] = float(val.group(1)) if val and val.group(1) else np.nan

    def axis_values(name):
        vals = []
        for chunk in re.findall(rf"/\*\s*{name}:\s*(.*?)\*/", text):
            vals.extend(float(v) for v in chunk.split())
        return np.array(vals)

    x_axis, y_axis = axis_values("x-axis"), axis_values("y-axis")

    mat = np.empty((nrow, ncol))
    for i, ln in enumerate(quoted[1 + ncolors:1 + ncolors + nrow]):
        row = ln[1:ln.rindex('"')]
        for j in range(ncol):
            mat[i, j] = symbol_value[row[j * cpp:(j + 1) * cpp]]
    mat = mat[::-1, :]  # xpm stores the top row (largest y) first
    # some versions write one axis value per bin edge; keep one per column
    x_axis = x_axis[:ncol] if x_axis.size >= ncol else np.arange(1, ncol + 1)
    y_axis = y_axis[:nrow] if y_axis.size >= nrow else np.arange(1, nrow + 1)
    return mat, x_axis, y_axis, title


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------
def running_mean(t, y, window_ns):
    """Centred running mean over a time window; NaN-padded at the edges."""
    if len(t) < 3:
        return y.copy()
    dt = np.median(np.diff(t))
    n = max(1, int(round(window_ns / dt))) if dt > 0 else 1
    if n < 2 or n >= len(y):
        return y.copy()
    kernel = np.ones(n) / n
    out = np.full_like(y, np.nan, dtype=float)
    valid = np.convolve(y, kernel, mode="valid")
    start = (n - 1) // 2
    out[start:start + len(valid)] = valid
    return out


def block_sem(y, nblocks=5):
    """Standard error from block averaging.

    Consecutive MD frames are correlated, so SD/sqrt(N) underestimates the
    uncertainty. Means of contiguous blocks are closer to independent. With
    a single short run this is still an approximation; independent replicas
    are the robust way to estimate uncertainty.
    """
    y = np.asarray(y)
    if len(y) < nblocks * 2:
        return np.nan
    blocks = np.array_split(y, nblocks)
    means = np.array([b.mean() for b in blocks])
    return means.std(ddof=1) / np.sqrt(nblocks)


def welch(a, b):
    """Welch's t-test between two sets of replica means.

    Returns (difference of means, SE of the difference, t, p). With three
    replicas per group the test has about 2-4 degrees of freedom, so only
    large, consistent differences reach significance; the point of the
    table is the difference and its uncertainty, not the p-value. p is nan
    when scipy is not installed or a group has fewer than two replicas.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    diff = a.mean() - b.mean()
    if len(a) < 2 or len(b) < 2:
        return diff, np.nan, np.nan, np.nan
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(va + vb)
    if se == 0:
        return diff, 0.0, np.nan, np.nan
    t = diff / se
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    p = 2 * _scipy_stats.t.sf(abs(t), df) if _scipy_stats else np.nan
    return diff, se, t, p


# -----------------------------------------------------------------------------
# Data discovery
# -----------------------------------------------------------------------------
def read_meta(rep_dir):
    meta = {}
    for line in (rep_dir / "meta.txt").read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            meta[k] = v
    return meta


def discover(ana_dir, only=None):
    """Return {ff: [replica_dir, ...]} for every analysed replica.

    Force fields are listed in FF_ORDER first, then any other label found
    under the analysis directory. Display names are taken from meta.txt.
    """
    labels = [p.name for p in sorted(ana_dir.iterdir()) if p.is_dir()] if ana_dir.exists() else []
    ordered = [f for f in FF_ORDER if f in labels] + [f for f in labels if f not in FF_ORDER]
    if only:
        ordered = [f for f in ordered if f in only]
    found = {}
    for ff in ordered:
        reps = sorted((ana_dir / ff).glob("rep*"), key=lambda p: int(re.sub(r"\D", "", p.name) or 0))
        reps = [r for r in reps if (r / "meta.txt").exists()]
        if reps:
            found[ff] = reps
            name = read_meta(reps[0]).get("display_name")
            if name:
                FF_NAME[ff] = name
            color_of(ff)
    return found


# Short panel titles for the compact distribution figure
SHORT_TITLE = {"rmsd": "Backbone RMSD", "rg": "Radius of gyration", "hb_intra": "H-bonds, intra",
               "hb_water": "H-bonds, protein-water", "sasa_heavy": "SASA, heavy atoms"}

# Time-series observables: key -> (file, column, axis label, panel title)
OBSERVABLES = [
    ("rmsd", "rmsd_backbone.xvg", 1, "RMSD (nm)", "Backbone RMSD vs minimized structure"),
    ("rg", "gyrate.xvg", 1, "Rg (nm)", "Radius of gyration"),
    ("hb_intra", "hbnum_intra.xvg", 1, "H-bonds", "Intramolecular H-bonds"),
    ("hb_water", "hbnum_water.xvg", 1, "H-bonds", "Protein-water H-bonds"),
    ("sasa_heavy", "sasa_heavy.xvg", 1, "SASA (nm$^2$)", "SASA, heavy atoms (comparable)"),
    ("sasa_native", "sasa_native.xvg", 1, "SASA (nm$^2$)", "SASA, all atoms of each model"),
]
# Observables reported in summary_table.csv (time series + mean RMSF)
TABLE_OBSERVABLES = ["rmsd", "rg", "rmsf_mean", "hb_intra", "hb_water", "sasa_heavy"]
TABLE_UNITS = {"rmsd": "nm", "rg": "nm", "rmsf_mean": "nm", "hb_intra": "count",
               "hb_water": "count", "sasa_heavy": "nm^2"}


def load_series(rep_dir, fname, col):
    path = rep_dir / fname
    if not path.exists():
        return None
    data, _ = read_xvg(path)
    if data.shape[1] <= col:
        return None
    return data[:, 0], data[:, col]


def load_by_legend(rep_dir, fname, legend):
    """Time series of one term in a multi-column xvg, located by its legend.

    gmx energy writes terms in energy-file order, not in the order requested,
    so columns must be identified by name.
    """
    path = rep_dir / fname
    if not path.exists():
        return None
    data, meta = read_xvg(path)
    names = [n.strip() for n in meta["legends"]]
    if legend not in names:
        return None
    col = names.index(legend) + 1
    return (data[:, 0], data[:, col]) if data.shape[1] > col else None


def window(t, y, skip_ns):
    """Equilibrated part of a series (whole series if it is too short)."""
    return y[t >= skip_ns] if np.any(t >= skip_ns) else y


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
def save(fig, fig_dir, name, formats):
    for ext in formats:
        fig.savefig(fig_dir / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote {fig_dir / name}.{{{','.join(formats)}}}")


def shade_skip(ax, skip_ns):
    if skip_ns > 0:
        ax.axvspan(0, skip_ns, color="#f0efec", zorder=0, lw=0)


def ff_legend(ax, ffs, loc="best", **kw):
    handles = [plt.Line2D([], [], color=color_of(f), lw=2, label=name_of(f)) for f in ffs]
    ax.legend(handles=handles, loc=loc, **kw)


def fig_timeseries(found, skip_ns, smooth_ns, fig_dir, formats):
    """Replica mean (running mean) per force field with a min-max band."""
    fig, axes = plt.subplots(3, 2, figsize=(10, 9.2), sharex=True)
    for ax, (key, fname, col, ylabel, title) in zip(axes.flat, OBSERVABLES):
        shade_skip(ax, skip_ns)
        plotted = False
        for ff, reps in found.items():
            curves = []
            for rep in reps:
                s = load_series(rep, fname, col)
                if s is not None:
                    curves.append((s[0], running_mean(s[0], s[1], smooth_ns)))
            if not curves:
                continue
            # replicas may differ in length (partial runs): align on the shortest
            n = min(len(c[0]) for c in curves)
            t = curves[0][0][:n]
            arr = np.vstack([c[1][:n] for c in curves])
            # the running mean is NaN-padded at the edges: those columns are
            # all-NaN and are simply not drawn
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                lo, hi, mean = np.nanmin(arr, axis=0), np.nanmax(arr, axis=0), np.nanmean(arr, axis=0)
            if len(curves) > 1:
                ax.fill_between(t, lo, hi, color=color_of(ff), alpha=0.18, lw=0)
            ax.plot(t, mean, color=color_of(ff), lw=1.5)
            plotted = True
        ax.set_title(title, loc="left", color=INK)
        ax.set_ylabel(ylabel)
        if not plotted:
            ax.text(0.5, 0.5, f"{fname} not found", ha="center", va="center",
                    transform=ax.transAxes, color=INK_2)
    for ax in axes[-1]:
        ax.set_xlabel("Time (ns)")
    ff_legend(axes[0, 0], list(found), loc="lower right")
    nrep = max(len(r) for r in found.values())
    fig.suptitle(f"Time series: {smooth_ns:g} ns running mean, averaged over replicas "
                 f"(band: min-max across replicas; gray: first {skip_ns:g} ns, excluded from averages)"
                 if nrep > 1 else
                 f"Time series: {smooth_ns:g} ns running mean (gray band: first {skip_ns:g} ns, excluded from averages)",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig1_timeseries", formats)


def fig_distributions(found, skip_ns, fig_dir, formats, stats_rows):
    keys = [o for o in OBSERVABLES if o[0] != "sasa_native"]
    fig, axes = plt.subplots(1, len(keys), figsize=(2.3 * len(keys), 2.8))
    for ax, (key, fname, col, xlabel, title) in zip(axes, keys):
        pooled = {}
        for ff, reps in found.items():
            vals = []
            for i, rep in enumerate(reps, start=1):
                s = load_series(rep, fname, col)
                if s is None:
                    continue
                t, y = s
                y = window(t, y, skip_ns)
                vals.append(y)
                stats_rows.append({
                    "observable": key, "force_field": ff, "replica": i,
                    "mean": y.mean(), "sd": y.std(ddof=1) if len(y) > 1 else np.nan,
                    "block_sem": block_sem(y), "n_frames": len(y),
                })
            if vals:
                pooled[ff] = np.concatenate(vals)
                rep_means = np.array([v.mean() for v in vals])
                stats_rows.append({
                    "observable": key, "force_field": ff, "replica": "all",
                    "mean": rep_means.mean(),
                    "sd": rep_means.std(ddof=1) if len(rep_means) > 1 else np.nan,
                    "block_sem": np.nan, "n_frames": len(pooled[ff]),
                })
        if pooled:
            lo = min(v.min() for v in pooled.values())
            hi = max(v.max() for v in pooled.values())
            # integer bins for counts (H-bonds), 40 bins otherwise
            if key.startswith("hb"):
                bins = np.arange(np.floor(lo) - 0.5, np.ceil(hi) + 1.5, 1.0)
            else:
                bins = np.linspace(lo, hi, 40) if hi > lo else 10
            for ff, v in pooled.items():
                ax.hist(v, bins=bins, density=True, histtype="step", lw=1.6, color=color_of(ff))
                ax.axvline(v.mean(), color=color_of(ff), lw=1, ls="--")
        ax.set_title(SHORT_TITLE.get(key, title), loc="left", color=INK, fontsize=8.5)
        ax.set_xlabel(xlabel)
        ax.set_yticks([])
        ax.grid(axis="y", visible=False)
    axes[0].set_ylabel("Probability density")
    handles = [plt.Line2D([], [], color=color_of(f), lw=2, label=name_of(f)) for f in found]
    fig.legend(handles=handles, loc="upper right", ncol=len(handles), bbox_to_anchor=(0.995, 1.02))
    fig.suptitle(f"Distributions after {skip_ns:g} ns, pooled over replicas (dashed: mean)",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig2_distributions", formats)


def per_residue_curves(found, fname, col):
    """{ff: (residues, [curve per replica])} for a per-residue xvg."""
    out = {}
    for ff, reps in found.items():
        curves = []
        for rep in reps:
            p = rep / fname
            if p.exists():
                d, _ = read_xvg(p)
                curves.append((d[:, 0], d[:, col]))
        if not curves:
            continue
        res = curves[0][0]
        out[ff] = (res, [c[1] for c in curves if len(c[0]) == len(res)])
    return out


def per_residue_mean(found, fname, col):
    """{ff: (residues, mean over replicas, sd over replicas or None)}"""
    out = {}
    for ff, (res, curves) in per_residue_curves(found, fname, col).items():
        arr = np.vstack(curves)
        out[ff] = (res, arr.mean(axis=0), arr.std(axis=0, ddof=1) if len(arr) > 1 else None)
    return out


def fig_per_residue(found, reference, fname, col, ylabel, title, diff_label, name, fig_dir, formats):
    data = per_residue_mean(found, fname, col)
    if not data:
        print(f"  skipped {name}: {fname} not found")
        return
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(10, 5.6), sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1]})
    for ff, (res, m, sd) in data.items():
        if sd is not None:
            ax.fill_between(res, m - sd, m + sd, color=color_of(ff), alpha=0.18, lw=0)
        ax.plot(res, m, color=color_of(ff), lw=1.5)
    ax.set_ylabel(ylabel)
    ax.set_title(title + (" (replica mean, band: SD over replicas)" if any(
        d[2] is not None for d in data.values()) else ""), loc="left", color=INK)
    ff_legend(ax, list(data), loc="upper right", ncol=min(len(data), 4))

    others = [f for f in data if f != reference]
    if reference in data and others:
        res_ref, m_ref, _ = data[reference]
        axd.axhline(0, color=INK_2, lw=0.6)
        drawn = 0
        for ff in others:
            res, m, _ = data[ff]
            if len(res) != len(res_ref):
                continue
            axd.plot(res, m - m_ref, color=color_of(ff), lw=1.3)
            axd.fill_between(res, 0, m - m_ref, color=color_of(ff), alpha=0.12, lw=0)
            drawn += 1
        axd.set_ylabel(diff_label)
        axd.set_title(f"Difference vs {name_of(reference)}", loc="left", color=INK, fontsize=9)
        if drawn == 0:
            axd.text(0.5, 0.5, "residue numbering differs between force fields",
                     ha="center", va="center", transform=axd.transAxes, color=INK_2)
    else:
        axd.text(0.5, 0.5, f"difference panel requires the reference force field ({reference})",
                 ha="center", va="center", transform=axd.transAxes, color=INK_2)
    axd.set_xlabel("Residue")
    fig.tight_layout()
    save(fig, fig_dir, name, formats)


def fig_contact_maps(found, reference, cutoff, fig_dir, formats):
    mats = {}
    for ff, reps in found.items():
        arrs = []
        for rep in reps:
            p = rep / "mdmat_mean.xpm"
            if p.exists():
                m, x, y, _ = read_xpm(p)
                arrs.append((m, x, y))
        if arrs:
            shape = arrs[0][0].shape
            stack = np.stack([a[0] for a in arrs if a[0].shape == shape])
            mats[ff] = (stack.mean(axis=0), arrs[0][1], arrs[0][2])
    if not mats:
        print("  skipped fig4_contact_maps: mdmat_mean.xpm not found")
        return

    others = [f for f in mats if f != reference]
    diffs = {}
    if reference in mats:
        m_ref = mats[reference][0]
        for ff in others:
            if mats[ff][0].shape == m_ref.shape:
                diffs[ff] = mats[ff][0] - m_ref
    ncols = len(mats)
    nrows = 2 if diffs else 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 4.2 * nrows), squeeze=False)

    for ax, (ff, (m, x, y)) in zip(axes[0], mats.items()):
        ext = [x[0] - 0.5, x[-1] + 0.5, y[0] - 0.5, y[-1] + 0.5]
        im = ax.imshow(m, origin="lower", extent=ext, cmap=SEQ_CMAP, vmin=0, vmax=cutoff,
                       interpolation="nearest", aspect="equal")
        ax.set_title(f"{name_of(ff)}: mean min. distance", loc="left", color=INK)
        ax.set_xlabel("Residue")
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("Distance (nm)")
    axes[0][0].set_ylabel("Residue")

    if diffs:
        lim = max(max(np.nanmax(np.abs(d)) for d in diffs.values()), 1e-3)
        x, y = mats[reference][1], mats[reference][2]
        ext = [x[0] - 0.5, x[-1] + 0.5, y[0] - 0.5, y[-1] + 0.5]
        for ax, (ff, d) in zip(axes[1], diffs.items()):
            im = ax.imshow(d, origin="lower", extent=ext, cmap=DIV_CMAP,
                           norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim),
                           interpolation="nearest", aspect="equal")
            ax.set_title(f"{name_of(ff)} minus {name_of(reference)}", loc="left", color=INK)
            ax.set_xlabel("Residue")
            ax.grid(False)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cb.set_label(f"$\\Delta$ distance (nm); <0: closer in {name_of(ff)}")
        axes[1][0].set_ylabel("Residue")
        for ax in axes[1][len(diffs):]:
            ax.axis("off")
    fig.suptitle(f"Residue contact maps, heavy atoms, replica mean (gmx mdmat, cutoff {cutoff:g} nm)",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig4_contact_maps", formats)


def fig_quality(found, ana_dir, fig_dir, formats):
    panels = [
        ("em", "Energy minimization", "Potential (kJ/mol)"),
        ("nvt_T", "NVT: temperature", "T (K)"),
        ("npt_P", "NPT: pressure", "P (bar)"),
        ("npt_rho", "NPT: density", "Density (kg/m$^3$)"),
        ("md_T", "Production: temperature", "T (K)"),
        ("md_rho", "Production: density", "Density (kg/m$^3$)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6))
    for ax, (key, title, ylabel) in zip(axes.flat, panels):
        for ff, reps in found.items():
            if key == "em":
                p = ana_dir / ff / "qc_em_potential.xvg"
                if p.exists():
                    d, _ = read_xvg(p)
                    # EM "time" is the step number; undo the ps -> ns conversion
                    ax.plot(d[:, 0] * 1000, d[:, 1], color=color_of(ff), lw=1.6)
                    ax.set_xlabel("EM step")
                continue
            for i, rep in enumerate(reps):
                fname, term = {
                    "nvt_T": ("qc_nvt_temperature.xvg", "Temperature"),
                    "npt_P": ("qc_npt_pressure.xvg", "Pressure"),
                    "npt_rho": ("qc_npt_density.xvg", "Density"),
                    "md_T": ("qc_md_thermo.xvg", "Temperature"),
                    "md_rho": ("qc_md_thermo.xvg", "Density"),
                }[key]
                s = load_by_legend(rep, fname, term)
                if s is None:
                    continue
                t, y = s
                if key.startswith(("nvt", "npt")):
                    t = t * 1000  # short runs: show ps
                ax.plot(t, y, color=color_of(ff), lw=0.8, alpha=0.9,
                        ls=REP_STYLE[i % len(REP_STYLE)])
                ax.set_xlabel("Time (ps)" if key.startswith(("nvt", "npt")) else "Time (ns)")
        ax.set_title(title, loc="left", color=INK)
        ax.set_ylabel(ylabel)
    ff_legend(axes[0, 1], list(found), loc="lower right")
    fig.suptitle("Equilibration and production quality checks (line style: replica)",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig6_quality_checks", formats)


def fig_replicas(found, skip_ns, smooth_ns, fig_dir, formats):
    """One row per force field, each replica as its own line: convergence
    and run-to-run spread within a force field, before any comparison
    between force fields."""
    cols = [o for o in OBSERVABLES if o[0] in ("rmsd", "rg", "hb_intra", "sasa_heavy")]
    nrows = len(found)
    fig, axes = plt.subplots(nrows, len(cols) + 1, figsize=(3.0 * (len(cols) + 1), 2.3 * nrows),
                             squeeze=False)
    rmsf = per_residue_curves(found, "rmsf_ca.xvg", 1)
    for i, (ff, reps) in enumerate(found.items()):
        for ax, (key, fname, col, ylabel, title) in zip(axes[i], cols):
            shade_skip(ax, skip_ns)
            for r, rep in enumerate(reps):
                s = load_series(rep, fname, col)
                if s is None:
                    continue
                t, y = s
                ax.plot(t, y, color=color_of(ff), lw=0.4, alpha=0.2)
                ax.plot(t, running_mean(t, y, smooth_ns), color=color_of(ff), lw=1.2,
                        ls=REP_STYLE[r % len(REP_STYLE)], label=f"rep {r + 1}")
            ax.set_ylabel(ylabel)
            if i == 0:
                ax.set_title(SHORT_TITLE.get(key, title), loc="left", color=INK)
            if i == nrows - 1:
                ax.set_xlabel("Time (ns)")
        axr = axes[i][-1]
        if ff in rmsf:
            res, curves = rmsf[ff]
            for r, c in enumerate(curves):
                axr.plot(res, c, color=color_of(ff), lw=1.0, ls=REP_STYLE[r % len(REP_STYLE)],
                         label=f"rep {r + 1}")
        axr.set_ylabel("RMSF (nm)")
        if i == 0:
            axr.set_title("C-alpha RMSF per residue", loc="left", color=INK)
        if i == nrows - 1:
            axr.set_xlabel("Residue")
        axr.legend(loc="upper right", fontsize=7.5, ncol=2)
        axes[i][0].text(-0.28, 0.5, name_of(ff), transform=axes[i][0].transAxes, rotation=90,
                        ha="center", va="center", color=color_of(ff), fontsize=10, fontweight="bold")
    fig.suptitle(f"Replicas within each force field ({smooth_ns:g} ns running mean; "
                 f"gray band: first {skip_ns:g} ns)", x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig7_replicas", formats)


# -----------------------------------------------------------------------------
# Tables
# -----------------------------------------------------------------------------
def rmsf_rows(found, stats_rows):
    """Mean C-alpha RMSF over residues per replica (the per-residue curves are
    plotted in fig3 and fig7); appended to stats_rows as observable rmsf_mean."""
    for ff, (res, curves) in per_residue_curves(found, "rmsf_ca.xvg", 1).items():
        means = []
        for i, c in enumerate(curves, start=1):
            means.append(c.mean())
            stats_rows.append({"observable": "rmsf_mean", "force_field": ff, "replica": i,
                               "mean": c.mean(), "sd": c.std(ddof=1), "block_sem": np.nan,
                               "n_frames": len(c)})
        means = np.array(means)
        stats_rows.append({"observable": "rmsf_mean", "force_field": ff, "replica": "all",
                           "mean": means.mean(), "sd": means.std(ddof=1) if len(means) > 1 else np.nan,
                           "block_sem": np.nan, "n_frames": len(res)})


def write_summary(rows, path):
    fields = ["observable", "force_field", "replica", "mean", "sd", "block_sem", "n_frames"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.4f}" if isinstance(r[k], float) else r[k]) for k in fields})
    print(f"  wrote {path}")


def write_table(rows, found, reference, path):
    """One row per force field and observable: mean +/- SD over replicas,
    plus the difference vs the reference with Welch's t-test."""
    ffs = list(found)
    rep_means = {}
    for key in TABLE_OBSERVABLES:
        for ff in ffs:
            vals = [r["mean"] for r in rows
                    if r["observable"] == key and r["force_field"] == ff and r["replica"] != "all"]
            if vals:
                rep_means[(key, ff)] = np.array(vals)

    fields = ["observable", "unit", "force_field", "n_replicas", "mean", "sd_replicas",
              "diff_vs_reference", "se_diff", "welch_t", "welch_p"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for key in TABLE_OBSERVABLES:
            for ff in ffs:
                if (key, ff) not in rep_means:
                    continue
                a = rep_means[(key, ff)]
                row = {"observable": key, "unit": TABLE_UNITS[key], "force_field": ff,
                       "n_replicas": len(a), "mean": f"{a.mean():.4f}",
                       "sd_replicas": f"{a.std(ddof=1):.4f}" if len(a) > 1 else "",
                       "diff_vs_reference": "", "se_diff": "", "welch_t": "", "welch_p": ""}
                if ff != reference and (key, reference) in rep_means:
                    d, se, t, p = welch(a, rep_means[(key, reference)])
                    row.update({"diff_vs_reference": f"{d:.4f}",
                                "se_diff": "" if np.isnan(se) else f"{se:.4f}",
                                "welch_t": "" if np.isnan(t) else f"{t:.2f}",
                                "welch_p": "" if np.isnan(p) else f"{p:.3f}"})
                w.writerow(row)
    print(f"  wrote {path}")

    # Console overview: mean +/- SD of replica means (block SEM for n = 1)
    colw = 24
    print("\n  observable      " + "".join(f"{name_of(f):>{colw}s}" for f in ffs))
    for key in TABLE_OBSERVABLES:
        line = f"  {key:<15s} "
        for ff in ffs:
            a = rep_means.get((key, ff))
            if a is None:
                line += f"{'-':>{colw}s}"
            elif len(a) == 1:
                sem = [r["block_sem"] for r in rows
                       if r["observable"] == key and r["force_field"] == ff and r["replica"] != "all"][0]
                line += f"{a[0]:>12.3f} +/- {sem:<7.3f}"
            else:
                line += f"{a.mean():>12.3f} +/- {a.std(ddof=1):<7.3f}"
        print(line)
    nrep = sorted({len(a) for a in rep_means.values()})
    print(f"  (n = {', '.join(map(str, nrep))} replicas; +/- SD of replica means, or block SEM when n = 1)")
    if reference in ffs:
        print(f"  differences vs {name_of(reference)} with Welch's t-test: {path.name}"
              + ("" if _scipy_stats else " (install scipy for p-values)"))


# -----------------------------------------------------------------------------
def main():
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description="Force field comparison plots (CHARMM36m, AMBER, GROMOS 53A6, GROMOS 54A7)")
    ap.add_argument("--analysis-dir", type=Path, default=root / "analysis")
    ap.add_argument("--fig-dir", type=Path, default=root / "figures")
    ap.add_argument("--reference", default="charmm36m",
                    help="force field label used for the difference panels and the table (default: charmm36m)")
    ap.add_argument("--only", nargs="+", default=None, help="restrict to these force field labels")
    ap.add_argument("--skip-ns", type=float, default=None,
                    help="exclude this initial time from averages (default: value used by 03_analysis.sh)")
    ap.add_argument("--smooth-ns", type=float, default=1.0, help="running-mean window (ns)")
    ap.add_argument("--mdmat-cutoff", type=float, default=None,
                    help="color scale limit for contact maps (default: cutoff used by 03_analysis.sh)")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    found = discover(args.analysis_dir, args.only)
    if not found:
        sys.exit(f"no analysed replicas under {args.analysis_dir}; run scripts/03_analysis.sh first")
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    first = read_meta(next(iter(found.values()))[0])
    if args.skip_ns is None:
        args.skip_ns = float(first.get("skip_ps", 0)) / 1000.0
    if args.mdmat_cutoff is None:
        args.mdmat_cutoff = float(first.get("mdmat_cutoff_nm", 1.5))

    print("replicas found: " + ", ".join(f"{name_of(f)} x{len(r)}" for f, r in found.items()))
    if len(found) < 2:
        print("  note: only one force field found; difference panels are skipped")
    if args.reference not in found:
        print(f"  note: reference '{args.reference}' not found; difference panels are skipped "
              f"(use --reference with one of: {', '.join(found)})")
    nrep = {len(r) for r in found.values()}
    if nrep != {1} and min(nrep) < 3:
        print("  note: fewer than 3 replicas for some force fields; SD over replicas is unreliable")

    stats_rows = []
    fig_timeseries(found, args.skip_ns, args.smooth_ns, args.fig_dir, args.formats)
    fig_distributions(found, args.skip_ns, args.fig_dir, args.formats, stats_rows)
    fig_per_residue(found, args.reference, "rmsf_ca.xvg", 1, "RMSF (nm)", "C-alpha RMSF per residue",
                    "$\\Delta$RMSF (nm)", "fig3_rmsf", args.fig_dir, args.formats)
    fig_contact_maps(found, args.reference, args.mdmat_cutoff, args.fig_dir, args.formats)
    fig_per_residue(found, args.reference, "sasa_heavy_res.xvg", 1, "SASA (nm$^2$)",
                    "Heavy-atom SASA per residue", "$\\Delta$SASA (nm$^2$)",
                    "fig5_sasa_per_residue", args.fig_dir, args.formats)
    fig_quality(found, args.analysis_dir, args.fig_dir, args.formats)
    fig_replicas(found, args.skip_ns, args.smooth_ns, args.fig_dir, args.formats)
    rmsf_rows(found, stats_rows)
    write_summary(stats_rows, args.fig_dir / "summary_statistics.csv")
    write_table(stats_rows, found, args.reference, args.fig_dir / "summary_table.csv")


if __name__ == "__main__":
    main()
