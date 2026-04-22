#!/usr/bin/env python3
"""
nmax_study.py
=============

Empirical study for choosing (dx, dz, nmax) for per-plane clustering of
photon showers.

For each clustering config (dx, dz):
  1. Randomly sample N events per ROOT file matching --input-glob.
  2. Cluster each event into per-plane (dx x dx) cells with plane_idx = floor(z/dz).
  3. For each candidate nmax, compute:
       - n_cells distribution (before truncation)
       - n_planes distribution (how many distinct z-planes are occupied)
       - fraction of events where n_cells > nmax (truncation triggered)
       - fractional energy loss from top-nmax truncation

Writes:
  <outdir>/nmax_study_results.csv      one row per (dx, dz, nmax, material)
  <outdir>/nmax_study_report.txt       human-readable summary
  <outdir>/plot_ncells_hist.png        histograms of n_cells per config
  <outdir>/plot_nplanes_hist.png       histograms of n_planes per config
  <outdir>/plot_energy_loss_vs_nmax.png curves of energy loss vs nmax per config

Usage:
  python /project/ctb-stelzer/hamza95/photons_gen/nmax_study.py \
      --input-glob "/project/ctb-stelzer/hamza95/photo_gen_files_*/*.root" \
      --events-per-file 20 \
      --outdir /scratch/hamza95/nmax_study_v2 \
      --seed 42
"""

import argparse
import glob
import os
import re
import sys
import time

import numpy as np
import uproot

import matplotlib
matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt


# =============================================================================
# Filename parsing — same regex as create_h5_file.py
# =============================================================================

FILENAME_RE = re.compile(
    r"photons_(\d+)x(\d+)x(\d+)cm_.*?_(PbWO4|PbF2)\.root$"
)


def parse_material(filepath):
    name = os.path.basename(filepath)
    match = FILENAME_RE.match(name)
    if not match:
        return None
    return match.group(4)


# =============================================================================
# Per-plane clustering (matches the proposed new cluster_event_arrays)
# =============================================================================

def cluster_event_per_plane(x, y, z, e, dx, dz):
    """
    Bin z into planes of thickness dz, cluster (x, y) with cell size (dx, dx)
    within each plane independently.

    Returns a (n_cells, 4) array of (x_center, y_center, plane_idx, energy_sum).
    Cells are NOT sorted or truncated — this is the full clustered output.
    """
    mask = e > 0
    x, y, z, e = x[mask], y[mask], z[mask], e[mask]
    if len(e) == 0:
        return np.zeros((0, 4), dtype=np.float32)

    plane_idx = np.floor(z / dz).astype(np.int32)

    out_chunks = []
    for p in np.unique(plane_idx):
        m = plane_idx == p
        xp, yp, ep = x[m], y[m], e[m]

        ix = np.floor(xp / dx).astype(np.int32)
        iy = np.floor(yp / dx).astype(np.int32)  # dy = dx

        x_min = int(ix.min())
        y_min = int(iy.min())
        stride = int(ix.max()) - x_min + 2

        keys = (ix.astype(np.int64) - x_min) * stride + (iy.astype(np.int64) - y_min)
        uniq, inv = np.unique(keys, return_inverse=True)

        e_sum = np.zeros(len(uniq), dtype=np.float32)
        np.add.at(e_sum, inv, ep)

        ux = (uniq // stride).astype(np.int32) + x_min
        uy = (uniq  % stride).astype(np.int32) + y_min

        xc = (ux + 0.5) * dx
        yc = (uy + 0.5) * dx  # dy = dx
        pc = np.full(len(uniq), p, dtype=np.float32)

        out_chunks.append(np.column_stack([xc, yc, pc, e_sum]).astype(np.float32))

    return np.concatenate(out_chunks, axis=0)


# =============================================================================
# Per-event evaluation
# =============================================================================

def evaluate_event(clustered, nmax_list):
    """
    Given a clustered (n_cells, 4) array, compute per-event metrics for each
    candidate nmax.

    Returns dict:
        n_cells       : int
        n_planes      : int
        total_energy  : float
        per_nmax[nm]  : {'truncated': bool, 'energy_loss_frac': float}
    """
    n_cells = len(clustered)
    if n_cells == 0:
        return None

    n_planes = len(np.unique(clustered[:, 2]))
    total_e = float(clustered[:, 3].sum())

    # Sort cells by energy descending, once
    order = np.argsort(clustered[:, 3])[::-1]
    e_sorted = clustered[order, 3]
    cum_e = np.cumsum(e_sorted)
    total_sorted = cum_e[-1]  # == total_e modulo fp error

    per_nmax = {}
    for nm in nmax_list:
        if n_cells <= nm:
            per_nmax[nm] = {"truncated": False, "energy_loss_frac": 0.0}
        else:
            kept_e = float(cum_e[nm - 1])
            loss = 1.0 - kept_e / total_sorted if total_sorted > 0 else 0.0
            per_nmax[nm] = {"truncated": True, "energy_loss_frac": float(loss)}

    return {
        "n_cells": n_cells,
        "n_planes": n_planes,
        "total_energy": total_e,
        "per_nmax": per_nmax,
    }


# =============================================================================
# ROOT file reader — random event sampling
# =============================================================================

def sample_events_from_root(path, n_events, rng):
    """
    Read up to n_events randomly-sampled events from one ROOT file.
    Yields dict with x, y, z, dE arrays per event.
    """
    try:
        tree = uproot.open(f"{path}:photon_sim")
    except Exception as e:
        print(f"  [WARN] cannot open {path}: {e}", file=sys.stderr)
        return

    try:
        all_ids = tree["EventID"].array(library="np")
    except Exception as e:
        print(f"  [WARN] cannot read EventID from {path}: {e}", file=sys.stderr)
        return

    n_hits = len(all_ids)
    if n_hits == 0:
        return

    # Find event boundaries (row indices where a new event starts)
    boundaries = np.flatnonzero(all_ids[1:] != all_ids[:-1]) + 1
    event_starts = np.concatenate(([0], boundaries, [n_hits])).astype(np.int64)
    n_events_total = len(event_starts) - 1

    if n_events_total == 0:
        return

    # Randomly choose which events to take
    take = min(n_events, n_events_total)
    chosen = rng.choice(n_events_total, size=take, replace=False)
    chosen.sort()  # sorted reads are faster

    # Read each chosen event's row range individually.
    # For 100 events this is fine — each event is a small contiguous slice.
    for ev in chosen:
        r0 = int(event_starts[ev])
        r1 = int(event_starts[ev + 1])
        if r1 <= r0:
            continue
        try:
            arr = tree.arrays(
                ["x", "y", "z", "dE"],
                entry_start=r0,
                entry_stop=r1,
                library="np",
            )
        except Exception as e:
            print(f"  [WARN] read error in {path} event {ev}: {e}", file=sys.stderr)
            continue
        yield {
            "x":  arr["x"],
            "y":  arr["y"],
            "z":  arr["z"],
            "dE": arr["dE"],
        }


# =============================================================================
# Statistics helpers
# =============================================================================

def percentiles(arr, ps=(50, 90, 95, 99, 99.9)):
    if len(arr) == 0:
        return {p: float("nan") for p in ps}
    return {p: float(np.percentile(arr, p)) for p in ps}


def summarize(values):
    """Return mean/median/p90/p95/p99/p99.9/max for a 1-D numeric iterable."""
    a = np.asarray(values, dtype=np.float64)
    if a.size == 0:
        return dict(mean=float("nan"), median=float("nan"),
                    p90=float("nan"), p95=float("nan"),
                    p99=float("nan"), p999=float("nan"),
                    max=float("nan"))
    return dict(
        mean   = float(a.mean()),
        median = float(np.median(a)),
        p90    = float(np.percentile(a, 90)),
        p95    = float(np.percentile(a, 95)),
        p99    = float(np.percentile(a, 99)),
        p999   = float(np.percentile(a, 99.9)),
        max    = float(a.max()),
    )


# =============================================================================
# Main study
# =============================================================================

def run_study(input_glob, events_per_file, configs, nmax_list, outdir, seed):
    """
    configs   : list of (dx, dz) tuples
    nmax_list : list of nmax values to evaluate for each config
    """
    os.makedirs(outdir, exist_ok=True)

    files = sorted(glob.glob(input_glob))
    if not files:
        raise RuntimeError(f"No ROOT files matched: {input_glob}")

    print(f"[INFO] Found {len(files)} ROOT files")
    print(f"[INFO] Sampling {events_per_file} events per file")
    print(f"[INFO] Clustering configs: {configs}")
    print(f"[INFO] nmax values: {nmax_list}")
    print(f"[INFO] Output dir: {outdir}")
    print()

    rng = np.random.default_rng(seed)

    # Storage: records[(dx, dz)] = list of per-event dicts
    records = {cfg: [] for cfg in configs}

    t0 = time.time()
    total_events = 0

    for fi, path in enumerate(files):
        material = parse_material(path)
        if material is None:
            print(f"  [SKIP] cannot parse material from: {os.path.basename(path)}")
            continue

        n_this_file = 0
        events = list(sample_events_from_root(path, events_per_file, rng))

        for ev in events:
            x, y, z, dE = ev["x"], ev["y"], ev["z"], ev["dE"]

            for (dx, dz) in configs:
                cl = cluster_event_per_plane(x, y, z, dE, dx, dz)
                res = evaluate_event(cl, nmax_list)
                if res is None:
                    continue
                res["material"] = material
                records[(dx, dz)].append(res)

            n_this_file += 1

        total_events += n_this_file
        dt = time.time() - t0
        print(f"  [{fi+1:3d}/{len(files)}] {os.path.basename(path):60s}  "
              f"{material}  {n_this_file:3d} events  (total {total_events}, {dt:.1f}s)")

    print()
    print(f"[INFO] Processed {total_events} events total in {time.time()-t0:.1f}s")
    print()

    # --- Aggregate into CSV rows ---
    csv_rows = []
    header = [
        "dx", "dz", "nmax", "material", "n_events",
        "ncells_mean", "ncells_median", "ncells_p90",
        "ncells_p95", "ncells_p99", "ncells_p999", "ncells_max",
        "nplanes_mean", "nplanes_median", "nplanes_p90",
        "nplanes_p95", "nplanes_p99", "nplanes_max",
        "frac_events_truncated",
        "eloss_mean", "eloss_median", "eloss_p90",
        "eloss_p95", "eloss_p99", "eloss_max",
    ]

    # Keep aggregates for plotting
    plot_data = {}  # plot_data[(dx, dz, material)] = {ncells: [...], nplanes: [...], eloss: {nm: [...]}}

    materials = sorted({r["material"] for recs in records.values() for r in recs})
    for (dx, dz) in configs:
        recs = records[(dx, dz)]
        for mat in materials + ["ALL"]:
            sub = [r for r in recs if mat == "ALL" or r["material"] == mat]
            if not sub:
                continue

            ncells_arr  = np.array([r["n_cells"]  for r in sub])
            nplanes_arr = np.array([r["n_planes"] for r in sub])

            plot_data.setdefault((dx, dz, mat), {
                "ncells": ncells_arr,
                "nplanes": nplanes_arr,
                "eloss": {},
            })

            for nm in nmax_list:
                trunc = np.array([r["per_nmax"][nm]["truncated"] for r in sub])
                eloss = np.array([r["per_nmax"][nm]["energy_loss_frac"] for r in sub])

                plot_data[(dx, dz, mat)]["eloss"][nm] = eloss

                nc = summarize(ncells_arr)
                npl = summarize(nplanes_arr)
                el = summarize(eloss)

                csv_rows.append([
                    dx, dz, nm, mat, len(sub),
                    nc["mean"], nc["median"], nc["p90"], nc["p95"], nc["p99"], nc["p999"], nc["max"],
                    npl["mean"], npl["median"], npl["p90"], npl["p95"], npl["p99"], npl["max"],
                    float(trunc.mean()),
                    el["mean"], el["median"], el["p90"], el["p95"], el["p99"], el["max"],
                ])

    # --- Write CSV ---
    csv_path = os.path.join(outdir, "nmax_study_results.csv")
    with open(csv_path, "w") as f:
        f.write(",".join(header) + "\n")
        for row in csv_rows:
            # Format floats to reasonable precision
            out = []
            for v in row:
                if isinstance(v, float):
                    out.append(f"{v:.6g}")
                else:
                    out.append(str(v))
            f.write(",".join(out) + "\n")
    print(f"[WRITE] {csv_path}")

    # --- Write human-readable report ---
    report_path = os.path.join(outdir, "nmax_study_report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("nmax / dx / dz study\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Input glob       : {input_glob}\n")
        f.write(f"Files found      : {len(files)}\n")
        f.write(f"Events / file    : {events_per_file} (random sample)\n")
        f.write(f"Total events     : {total_events}\n")
        f.write(f"Random seed      : {seed}\n")
        f.write(f"Configs          : {configs}\n")
        f.write(f"nmax values      : {nmax_list}\n\n")

        for (dx, dz) in configs:
            recs = records[(dx, dz)]
            if not recs:
                continue
            f.write("-" * 80 + "\n")
            f.write(f"Config: dx={dx} cm, dy={dx} cm, dz={dz} cm   "
                    f"(n_events={len(recs)})\n")
            f.write("-" * 80 + "\n")

            ncells  = np.array([r["n_cells"]  for r in recs])
            nplanes = np.array([r["n_planes"] for r in recs])

            f.write("\n  n_cells (before truncation):\n")
            nc = summarize(ncells)
            f.write(f"    mean={nc['mean']:.1f}  median={nc['median']:.0f}  "
                    f"p90={nc['p90']:.0f}  p95={nc['p95']:.0f}  "
                    f"p99={nc['p99']:.0f}  p99.9={nc['p999']:.0f}  "
                    f"max={nc['max']:.0f}\n")

            f.write("\n  n_planes (occupied z-planes per event):\n")
            npl = summarize(nplanes)
            f.write(f"    mean={npl['mean']:.1f}  median={npl['median']:.0f}  "
                    f"p90={npl['p90']:.0f}  p95={npl['p95']:.0f}  "
                    f"p99={npl['p99']:.0f}  max={npl['max']:.0f}\n")

            # Per-material breakdown
            for mat in materials:
                sub = [r for r in recs if r["material"] == mat]
                if not sub:
                    continue
                f.write(f"\n  [{mat}] n_events={len(sub)}\n")
                nc_m = summarize([r["n_cells"]  for r in sub])
                np_m = summarize([r["n_planes"] for r in sub])
                f.write(f"    n_cells  : mean={nc_m['mean']:.1f}  "
                        f"median={nc_m['median']:.0f}  "
                        f"p99={nc_m['p99']:.0f}  max={nc_m['max']:.0f}\n")
                f.write(f"    n_planes : mean={np_m['mean']:.1f}  "
                        f"median={np_m['median']:.0f}  max={np_m['max']:.0f}\n")

            f.write("\n  Energy loss from top-nmax truncation:\n")
            f.write(f"    {'nmax':>6s}  {'%trunc':>8s}  "
                    f"{'mean':>10s}  {'median':>10s}  {'p95':>10s}  "
                    f"{'p99':>10s}  {'max':>10s}\n")
            for nm in nmax_list:
                trunc = np.array([r["per_nmax"][nm]["truncated"] for r in recs])
                eloss = np.array([r["per_nmax"][nm]["energy_loss_frac"] for r in recs])
                s = summarize(eloss)
                f.write(f"    {nm:>6d}  {100*trunc.mean():>7.1f}%  "
                        f"{100*s['mean']:>9.3f}%  "
                        f"{100*s['median']:>9.3f}%  "
                        f"{100*s['p95']:>9.3f}%  "
                        f"{100*s['p99']:>9.3f}%  "
                        f"{100*s['max']:>9.3f}%\n")

            # Per-material energy loss
            for mat in materials:
                sub = [r for r in recs if r["material"] == mat]
                if not sub:
                    continue
                f.write(f"\n  [{mat}] Energy loss:\n")
                f.write(f"    {'nmax':>6s}  {'%trunc':>8s}  "
                        f"{'mean':>10s}  {'median':>10s}  {'p99':>10s}\n")
                for nm in nmax_list:
                    trunc = np.array([r["per_nmax"][nm]["truncated"] for r in sub])
                    eloss = np.array([r["per_nmax"][nm]["energy_loss_frac"] for r in sub])
                    s = summarize(eloss)
                    f.write(f"    {nm:>6d}  {100*trunc.mean():>7.1f}%  "
                            f"{100*s['mean']:>9.3f}%  "
                            f"{100*s['median']:>9.3f}%  "
                            f"{100*s['p99']:>9.3f}%\n")
            f.write("\n")

    print(f"[WRITE] {report_path}")

    # --- Plots ---
    make_plots(plot_data, configs, nmax_list, materials, outdir)


# =============================================================================
# Plotting
# =============================================================================

def make_plots(plot_data, configs, nmax_list, materials, outdir):

    colors = {"PbF2": "tab:blue", "PbWO4": "tab:orange", "ALL": "black"}

    # ---- Plot 1: n_cells histogram per config ----
    fig, axes = plt.subplots(1, len(configs), figsize=(6 * len(configs), 4.5),
                             squeeze=False)
    for i, (dx, dz) in enumerate(configs):
        ax = axes[0, i]
        for mat in materials:
            key = (dx, dz, mat)
            if key not in plot_data:
                continue
            arr = plot_data[key]["ncells"]
            if len(arr) == 0:
                continue
            ax.hist(arr, bins=50, histtype="step", linewidth=1.5,
                    label=f"{mat} (n={len(arr)})", color=colors.get(mat))

        # Reference lines for nmax values
        for nm in nmax_list:
            ax.axvline(nm, linestyle="--", linewidth=0.8, alpha=0.6,
                       color="gray")
            ax.text(nm, ax.get_ylim()[1] * 0.95, f"nmax={nm}",
                    rotation=90, va="top", ha="right",
                    fontsize=8, color="gray")

        ax.set_xlabel("n_cells (before truncation)")
        ax.set_ylabel("events")
        ax.set_title(f"dx={dx} cm,  dz={dz} cm")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle("Distribution of n_cells per event (before truncation)",
                 fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "plot_ncells_hist.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT]  {p}")

    # ---- Plot 2: n_planes histogram per config ----
    fig, axes = plt.subplots(1, len(configs), figsize=(6 * len(configs), 4.5),
                             squeeze=False)
    for i, (dx, dz) in enumerate(configs):
        ax = axes[0, i]
        for mat in materials:
            key = (dx, dz, mat)
            if key not in plot_data:
                continue
            arr = plot_data[key]["nplanes"]
            if len(arr) == 0:
                continue
            # Integer-spaced bins
            bins = np.arange(arr.min(), arr.max() + 2) - 0.5
            ax.hist(arr, bins=bins, histtype="step", linewidth=1.5,
                    label=f"{mat} (n={len(arr)})", color=colors.get(mat))
        ax.set_xlabel("n_planes (occupied z-slabs per event)")
        ax.set_ylabel("events")
        ax.set_title(f"dx={dx} cm,  dz={dz} cm")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle("Distribution of occupied planes per event", fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "plot_nplanes_hist.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT]  {p}")

    # ---- Plot 3: energy loss curves vs nmax ----
    fig, axes = plt.subplots(1, len(configs), figsize=(6 * len(configs), 4.5),
                             squeeze=False, sharey=True)
    for i, (dx, dz) in enumerate(configs):
        ax = axes[0, i]
        for mat in materials + ["ALL"]:
            key = (dx, dz, mat)
            if key not in plot_data:
                continue
            losses = plot_data[key]["eloss"]
            mean_curve = [100 * losses[nm].mean()  for nm in nmax_list]
            p99_curve  = [100 * np.percentile(losses[nm], 99) for nm in nmax_list]

            ls_mean = "-"   if mat != "ALL" else "--"
            ls_p99  = ":"
            lbl = f"{mat} mean" if mat != "ALL" else "ALL mean"
            ax.plot(nmax_list, mean_curve, linestyle=ls_mean, marker="o",
                    color=colors.get(mat), label=lbl)
            if mat != "ALL":
                ax.plot(nmax_list, p99_curve, linestyle=ls_p99, marker="x",
                        color=colors.get(mat), alpha=0.7,
                        label=f"{mat} p99")
        ax.set_xlabel("nmax")
        ax.set_ylabel("Energy loss from truncation (%)")
        ax.set_title(f"dx={dx} cm,  dz={dz} cm")
        ax.set_xticks(nmax_list)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("Fractional energy loss vs nmax (mean + 99th percentile)",
                 fontsize=13)
    fig.tight_layout()
    p = os.path.join(outdir, "plot_energy_loss_vs_nmax.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT]  {p}")


# =============================================================================
# CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Study how (dx, dz, nmax) affect per-plane shower clustering."
    )
    ap.add_argument("--input-glob", required=True,
                    help="Glob pattern for ROOT files "
                         "(e.g. '/project/ctb-stelzer/hamza95/photo_gen_files_*/*.root')")
    ap.add_argument("--events-per-file", type=int, default=100,
                    help="Random events sampled per ROOT file (default: 100)")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for CSV, TXT, and PNG files")
    ap.add_argument("--seed", type=int, default=42,
                    help="Random seed (default: 42)")
    args = ap.parse_args()

    # Fixed study grid per user request
    configs = [
        (5., 10.0),   # fine
        (10.0, 20.0),   # coarse
    ]
    nmax_list = [2048, 4096, 6016]

    run_study(
        input_glob      = args.input_glob,
        events_per_file = args.events_per_file,
        configs         = configs,
        nmax_list       = nmax_list,
        outdir          = args.outdir,
        seed            = args.seed,
    )

    print()
    print("[DONE] See outputs in:", args.outdir)


if __name__ == "__main__":
    main()