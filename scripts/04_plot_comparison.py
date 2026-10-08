#!/usr/bin/env python3
# =============================================================================
# 04_plot_comparison.py
# Comparative figures and summary statistics, CHARMM36(m) vs GROMOS 53A6,
# from the outputs of 03_analysis.sh (analysis/<ff>/rep<N>/).
#
# Figures (figures/):
#   fig1_timeseries          RMSD, Rg, H-bonds (intra, water), SASA, contacts
#   fig2_distributions       distributions over the equilibrated window
#   fig3_rmsf                C-alpha RMSF per residue + difference
#   fig4_contact_maps        mean residue distance maps + difference
#   fig5_sasa_per_residue    heavy-atom SASA per residue + difference
#   fig6_quality_checks      EM, NVT, NPT and production thermodynamics
# Table:
#   summary_statistics.csv   mean, SD and block-averaged SEM per observable
#
# Usage:
#   python3 scripts/04_plot_comparison.py
#   python3 scripts/04_plot_comparison.py --smooth-ns 2 --formats png pdf
# Requires: numpy, matplotlib
# =============================================================================

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless backend: works on servers without a display
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------
# Fixed categorical order (validated colorblind-safe pair): the force field
# always keeps its color, whatever subset of data is plotted.
FF_ORDER = ["charmm36", "gromos53a6"]
FF_NAME = {"charmm36": "CHARMM36m", "gromos53a6": "GROMOS 53A6"}
FF_COLOR = {"charmm36": "#2a78d6", "gromos53a6": "#eb6834"}

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


# -----------------------------------------------------------------------------
# Data discovery
# -----------------------------------------------------------------------------
def discover(ana_dir):
    """Return {ff: [replica_dir, ...]} for every analysed replica."""
    found = {}
    for ff in FF_ORDER:
        reps = sorted((ana_dir / ff).glob("rep*"), key=lambda p: int(re.sub(r"\D", "", p.name) or 0))
        reps = [r for r in reps if (r / "meta.txt").exists()]
        if reps:
            found[ff] = reps
    return found


def read_meta(rep_dir):
    meta = {}
    for line in (rep_dir / "meta.txt").read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            meta[k] = v
    return meta


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


def ff_legend(ax, ffs, loc="best"):
    handles = [plt.Line2D([], [], color=FF_COLOR[f], lw=2, label=FF_NAME[f]) for f in ffs]
    ax.legend(handles=handles, loc=loc)


def fig_timeseries(found, skip_ns, smooth_ns, fig_dir, formats):
    fig, axes = plt.subplots(3, 2, figsize=(10, 9.2), sharex=True)
    for ax, (key, fname, col, ylabel, title) in zip(axes.flat, OBSERVABLES):
        shade_skip(ax, skip_ns)
        plotted = False
        for ff, reps in found.items():
            for rep in reps:
                s = load_series(rep, fname, col)
                if s is None:
                    continue
                t, y = s
                ax.plot(t, y, color=FF_COLOR[ff], lw=0.5, alpha=0.25)
                ax.plot(t, running_mean(t, y, smooth_ns), color=FF_COLOR[ff], lw=1.6)
                plotted = True
        ax.set_title(title, loc="left", color=INK)
        ax.set_ylabel(ylabel)
        if not plotted:
            ax.text(0.5, 0.5, f"{fname} not found", ha="center", va="center",
                    transform=ax.transAxes, color=INK_2)
    for ax in axes[-1]:
        ax.set_xlabel("Time (ns)")
    ff_legend(axes[0, 0], list(found), loc="lower right")
    fig.suptitle(f"Time series (thin: raw; thick: {smooth_ns:g} ns running mean; "
                 f"gray band: first {skip_ns:g} ns, excluded from averages)",
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
                y = y[t >= skip_ns] if np.any(t >= skip_ns) else y
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
                ax.hist(v, bins=bins, density=True, histtype="step", lw=1.6, color=FF_COLOR[ff])
                ax.axvline(v.mean(), color=FF_COLOR[ff], lw=1, ls="--")
        ax.set_title(SHORT_TITLE.get(key, title), loc="left", color=INK, fontsize=8.5)
        ax.set_xlabel(xlabel)
        ax.set_yticks([])
        ax.grid(axis="y", visible=False)
    axes[0].set_ylabel("Probability density")
    handles = [plt.Line2D([], [], color=FF_COLOR[f], lw=2, label=FF_NAME[f]) for f in found]
    fig.legend(handles=handles, loc="upper right", ncol=len(handles), bbox_to_anchor=(0.995, 1.02))
    fig.suptitle(f"Distributions after {skip_ns:g} ns (dashed: mean)",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig2_distributions", formats)


def per_residue_mean(found, fname, col):
    """{ff: (residues, mean over replicas, sd over replicas or None)}"""
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
        same = [c[1] for c in curves if len(c[0]) == len(res)]
        arr = np.vstack(same)
        out[ff] = (res, arr.mean(axis=0), arr.std(axis=0, ddof=1) if len(arr) > 1 else None)
    return out


def fig_per_residue(found, fname, col, ylabel, title, diff_label, name, fig_dir, formats):
    data = per_residue_mean(found, fname, col)
    if not data:
        print(f"  skipped {name}: {fname} not found")
        return
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(10, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1]})
    for ff, (res, m, sd) in data.items():
        if sd is not None:
            ax.fill_between(res, m - sd, m + sd, color=FF_COLOR[ff], alpha=0.18, lw=0)
        ax.plot(res, m, color=FF_COLOR[ff], lw=1.6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", color=INK)
    ff_legend(ax, list(data), loc="upper right")

    if all(f in data for f in FF_ORDER) and len(data[FF_ORDER[0]][0]) == len(data[FF_ORDER[1]][0]):
        res = data[FF_ORDER[0]][0]
        diff = data[FF_ORDER[0]][1] - data[FF_ORDER[1]][1]
        axd.bar(res, diff, width=0.85, color=np.where(diff >= 0, FF_COLOR["charmm36"], FF_COLOR["gromos53a6"]))
        axd.axhline(0, color=INK_2, lw=0.6)
        axd.set_ylabel(diff_label)
    else:
        axd.text(0.5, 0.5, "difference requires both force fields with identical residue numbering",
                 ha="center", va="center", transform=axd.transAxes, color=INK_2)
    axd.set_xlabel("Residue")
    fig.tight_layout()
    save(fig, fig_dir, name, formats)


def fig_contact_maps(found, cutoff, fig_dir, formats):
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

    both = all(f in mats for f in FF_ORDER) and mats[FF_ORDER[0]][0].shape == mats[FF_ORDER[1]][0].shape
    ncols = 3 if both else len(mats)
    fig, axes = plt.subplots(1, ncols, figsize=(4.3 * ncols, 4.2))
    axes = np.atleast_1d(axes)
    for ax, (ff, (m, x, y)) in zip(axes, mats.items()):
        ext = [x[0] - 0.5, x[-1] + 0.5, y[0] - 0.5, y[-1] + 0.5]
        im = ax.imshow(m, origin="lower", extent=ext, cmap=SEQ_CMAP, vmin=0, vmax=cutoff,
                       interpolation="nearest", aspect="equal")
        ax.set_title(f"{FF_NAME[ff]}: mean min. distance", loc="left", color=INK)
        ax.set_xlabel("Residue")
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("Distance (nm)")
    axes[0].set_ylabel("Residue")
    if both:
        m1, x, y = mats[FF_ORDER[0]]
        diff = m1 - mats[FF_ORDER[1]][0]
        lim = max(np.nanmax(np.abs(diff)), 1e-3)
        ax = axes[2]
        ext = [x[0] - 0.5, x[-1] + 0.5, y[0] - 0.5, y[-1] + 0.5]
        im = ax.imshow(diff, origin="lower", extent=ext, cmap=DIV_CMAP,
                       norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim),
                       interpolation="nearest", aspect="equal")
        ax.set_title("CHARMM36m minus GROMOS 53A6", loc="left", color=INK)
        ax.set_xlabel("Residue")
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("$\\Delta$ distance (nm); <0: closer in CHARMM")
    fig.suptitle(f"Residue contact maps, heavy atoms (gmx mdmat, cutoff {cutoff:g} nm)",
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
                    ax.plot(d[:, 0] * 1000, d[:, 1], color=FF_COLOR[ff], lw=1.6)
                    ax.set_xlabel("EM step")
                continue
            for rep in reps:
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
                ax.plot(t, y, color=FF_COLOR[ff], lw=0.8, alpha=0.9)
                ax.set_xlabel("Time (ps)" if key.startswith(("nvt", "npt")) else "Time (ns)")
        ax.set_title(title, loc="left", color=INK)
        ax.set_ylabel(ylabel)
    ff_legend(axes[0, 1], list(found), loc="lower right")
    fig.suptitle("Equilibration and production quality checks",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    fig.tight_layout()
    save(fig, fig_dir, "fig6_quality_checks", formats)


def write_summary(rows, path):
    fields = ["observable", "force_field", "replica", "mean", "sd", "block_sem", "n_frames"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.4f}" if isinstance(r[k], float) else r[k]) for k in fields})
    print(f"  wrote {path}")

    # Console overview: per force field, pooled over replicas
    print("\n  observable      " + "".join(f"{FF_NAME[f]:>24s}" for f in FF_ORDER))
    for key in dict.fromkeys(r["observable"] for r in rows):
        line = f"  {key:<15s} "
        for ff in FF_ORDER:
            reps = [r for r in rows if r["observable"] == key and r["force_field"] == ff and r["replica"] != "all"]
            if not reps:
                line += f"{'-':>24s}"
                continue
            if len(reps) == 1:
                r = reps[0]
                line += f"{r['mean']:>12.3f} +/- {r['block_sem']:<7.3f}"
            else:
                agg = [r for r in rows if r["observable"] == key and r["force_field"] == ff and r["replica"] == "all"][0]
                line += f"{agg['mean']:>12.3f} +/- {agg['sd']:<7.3f}"
        print(line)
    print("  (one replica: +/- block SEM; several replicas: +/- SD of replica means)")


# -----------------------------------------------------------------------------
def main():
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description="CHARMM36m vs GROMOS 53A6 comparison plots")
    ap.add_argument("--analysis-dir", type=Path, default=root / "analysis")
    ap.add_argument("--fig-dir", type=Path, default=root / "figures")
    ap.add_argument("--skip-ns", type=float, default=None,
                    help="exclude this initial time from averages (default: value used by 03_analysis.sh)")
    ap.add_argument("--smooth-ns", type=float, default=1.0, help="running-mean window (ns)")
    ap.add_argument("--mdmat-cutoff", type=float, default=None,
                    help="color scale limit for contact maps (default: cutoff used by 03_analysis.sh)")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    found = discover(args.analysis_dir)
    if not found:
        sys.exit(f"no analysed replicas under {args.analysis_dir}; run scripts/03_analysis.sh first")
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    first = read_meta(next(iter(found.values()))[0])
    if args.skip_ns is None:
        args.skip_ns = float(first.get("skip_ps", 0)) / 1000.0
    if args.mdmat_cutoff is None:
        args.mdmat_cutoff = float(first.get("mdmat_cutoff_nm", 1.5))

    print("replicas found: " + ", ".join(f"{FF_NAME[f]} x{len(r)}" for f, r in found.items()))
    if len(found) < 2:
        print("  note: only one force field found; difference panels are skipped")

    stats_rows = []
    fig_timeseries(found, args.skip_ns, args.smooth_ns, args.fig_dir, args.formats)
    fig_distributions(found, args.skip_ns, args.fig_dir, args.formats, stats_rows)
    fig_per_residue(found, "rmsf_ca.xvg", 1, "RMSF (nm)", "C-alpha RMSF per residue",
                    "$\\Delta$RMSF (nm)\nCHARMM - GROMOS", "fig3_rmsf", args.fig_dir, args.formats)
    fig_contact_maps(found, args.mdmat_cutoff, args.fig_dir, args.formats)
    fig_per_residue(found, "sasa_heavy_res.xvg", 1, "SASA (nm$^2$)",
                    "Heavy-atom SASA per residue", "$\\Delta$SASA (nm$^2$)\nCHARMM - GROMOS",
                    "fig5_sasa_per_residue", args.fig_dir, args.formats)
    fig_quality(found, args.analysis_dir, args.fig_dir, args.formats)
    write_summary(stats_rows, args.fig_dir / "summary_statistics.csv")


if __name__ == "__main__":
    main()
