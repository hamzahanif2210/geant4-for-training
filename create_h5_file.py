#!/usr/bin/env python3
"""
create_h5_file.py
=================

Convert photon-shower ROOT files into per-event HDF5 using per-plane
clustering with variable-length (vlen) storage per event.

For each event:
  1. Bin z into planes of thickness `dz`         → plane_idx = floor(z/dz)
  2. Within each plane, bin (x, y) into (dx, dx) cells; sum energy per cell
  3. If more than `nmax` cells, keep the top-`nmax` by energy; drop the rest
  4. Sort kept cells by plane_idx  (shallow → deep)
  5. Store flat in a vlen float32 dataset; no padding

Per-cell features (4 columns, flattened per event):
    x_center    (cm)
    y_center    (cm)
    plane_idx   (integer, stored as float32 for uniform dtype)
    energy_sum

Physical z is reconstructable as (plane_idx + 0.5) * dz.

Usage
-----
Single-file:
    python create_h5_file.py --input-root file.root --output file.h5 \
        --dx 10 --dz 20 --nmax 4096

Combine per-file H5s:
    python create_h5_file.py --combine-only \
        --input-glob "/path/per_file/*.h5" --output combined.h5

Slurm submission (array of mappers + dependent combine):
    python /project/ctb-stelzer/hamza95/photons_gen/create_h5_file.py \
        --submit-slurm \
        --input-glob "/project/ctb-stelzer/hamza95/photo_gen_files_*/*.root" \
        --output /scratch/hamza95/all_clustered_showers.h5 \
        --work-dir /project/ctb-stelzer/hamza95/photo_gen_slurm_h5_work \
        --per-file-dir /scratch/hamza95/per_root_h5 \
        --dx 10 --dz 20 --nmax 4096 \
        --account def-mdanning
"""

import argparse
import fcntl
import glob
import os
import re
import subprocess

import h5py
import numpy as np
import uproot


# =============================================================================
# Defaults
# =============================================================================

DEFAULT_DX   = 10.0    # cm  (also used for dy)
DEFAULT_DZ   = 20.0    # cm
DEFAULT_NMAX = 4096

N_FEAT = 4             # (x, y, plane_idx, energy)


# =============================================================================
# Filename parsing — only material is used; cell size comes from CLI
# =============================================================================

FILENAME_RE = re.compile(
    r"photons_(\d+)x(\d+)x(\d+)cm_.*?_(PbWO4|PbF2)\.root$"
)


def parse_filename_material(filepath):
    name = os.path.basename(filepath)
    match = FILENAME_RE.match(name)
    if not match:
        raise ValueError(f"Bad filename (material not parseable): {name}")
    return match.group(4)


# =============================================================================
# Per-plane clustering
# =============================================================================

def cluster_event_per_plane(x, y, z, e, dx, dz):
    """
    Bin z into planes of thickness dz, cluster (x, y) with cell size (dx, dx)
    within each plane independently.

    Returns (n_cells, 4) float32:
        col 0: x_center (cm)
        col 1: y_center (cm)
        col 2: plane_idx (integer, stored as float32)
        col 3: energy_sum

    Cells are NOT sorted or truncated here — caller handles that.
    """
    mask = e > 0
    x, y, z, e = x[mask], y[mask], z[mask], e[mask]
    if len(e) == 0:
        return np.zeros((0, N_FEAT), dtype=np.float32)

    plane_idx = np.floor(z / dz).astype(np.int32)

    out_chunks = []
    # Iterate over each occupied plane separately — hits in different planes
    # never merge, which is the whole point of "per-plane" clustering.
    for p in np.unique(plane_idx):
        m = plane_idx == p
        xp, yp, ep = x[m], y[m], e[m]

        ix = np.floor(xp / dx).astype(np.int32)
        iy = np.floor(yp / dx).astype(np.int32)  # dy = dx

        # Fast 2D unique via a 1-D key (TAMBO trick): encode (ix, iy) into a
        # single int64 so np.unique runs on a 1-D array.
        x_min = int(ix.min())
        y_min = int(iy.min())
        stride = int(ix.max()) - x_min + 2  # +1 for inclusive range, +1 margin
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


def truncate_and_sort(clustered, nmax):
    """
    Keep top-`nmax` cells by energy (col 3), then sort by plane_idx (col 2)
    so the output is in shallow-to-deep z order.

    Returns array of shape (n_kept, N_FEAT) with n_kept == min(len(clustered), nmax).
    No padding — caller uses n_kept to know the real size.
    """
    n = len(clustered)
    if n == 0:
        return np.zeros((0, N_FEAT), dtype=np.float32)

    if n > nmax:
        # Partial partition is O(n) and faster than a full argsort.
        top = np.argpartition(clustered[:, 3], -nmax)[-nmax:]
        kept = clustered[top]
    else:
        kept = clustered

    # Sort by plane_idx (stable, so cells within a plane keep deterministic order)
    order = np.argsort(kept[:, 2], kind="stable")
    return kept[order].astype(np.float32)


# =============================================================================
# Write helpers — buffered per-event append
# =============================================================================

def _flush_events(datasets, buffers, total):
    """Append buffered events to all datasets with one resize each."""
    showers_ds, n_cells_ds, event_id_ds, primaryE_ds, material_ds = datasets
    b_show, b_nc, b_eid, b_E0, b_mat = buffers

    k = len(b_nc)
    if k == 0:
        return total

    new = total + k

    showers_ds.resize((new,))
    n_cells_ds.resize((new,))
    event_id_ds.resize((new,))
    primaryE_ds.resize((new,))
    material_ds.resize((new,))

    # vlen dataset requires per-element assignment.
    for i, flat in enumerate(b_show):
        showers_ds[total + i] = flat

    n_cells_ds[total:new]  = np.asarray(b_nc,  dtype=np.int32)
    event_id_ds[total:new] = np.asarray(b_eid, dtype=np.int64)
    primaryE_ds[total:new] = np.asarray(b_E0,  dtype=np.float32)
    material_ds[total:new] = np.asarray(b_mat, dtype=object)

    return new


# =============================================================================
# Per-ROOT processing
# =============================================================================

def _write_attrs(f, dx, dz, nmax, max_planes):
    """Write self-describing attributes on an h5 file."""
    f.attrs["dx"]              = float(dx)
    f.attrs["dy"]              = float(dx)
    f.attrs["dz"]              = float(dz)
    f.attrs["nmax"]            = int(nmax)
    f.attrs["n_planes_max"]    = int(max_planes)
    f.attrs["feature_columns"] = "x,y,plane_idx,energy"


# =============================================================================
# Per-file summary line (lock-protected, safe for parallel Slurm tasks)
# =============================================================================

SUMMARY_HEADER = (
    "# Per-file clustering summary. Columns:\n"
    "#   root_file              : basename of the ROOT input\n"
    "#   material               : PbF2 or PbWO4\n"
    "#   n_events               : events processed\n"
    "#   n_clusters_raw         : total cells across all events BEFORE truncation\n"
    "#   n_clusters_kept        : total cells across all events AFTER truncation (≤ nmax*n_events)\n"
    "#   n_events_truncated     : events where n_raw > nmax (truncation triggered)\n"
    "#   max_n_cells_raw        : max cells-per-event before truncation\n"
    "#   max_n_cells_kept       : max cells-per-event after truncation\n"
    "#   energy_raw             : sum of all cell energies before truncation\n"
    "#   energy_kept            : sum of all cell energies after truncation\n"
    "#   energy_retention_pct   : 100 * energy_kept / energy_raw\n"
    "#\n"
    f"{'root_file':<80s} {'material':>8s} {'n_events':>10s} "
    f"{'n_clust_raw':>14s} {'n_clust_kept':>14s} {'n_ev_trunc':>11s} "
    f"{'max_nc_raw':>11s} {'max_nc_kept':>12s} "
    f"{'energy_raw':>14s} {'energy_kept':>14s} {'retention%':>11s}\n"
)


def _append_summary_line(summary_path, line):
    """
    Append one line to a shared summary text file, protected by fcntl.flock
    so parallel Slurm array tasks don't interleave their writes.

    If the file doesn't exist yet (we acquire the lock and find it empty),
    write the header first.
    """
    if summary_path is None:
        return

    os.makedirs(os.path.dirname(os.path.abspath(summary_path)), exist_ok=True)

    # Open in append mode; create if absent. fcntl.flock on the same fd
    # guarantees mutual exclusion across processes on the same host. On
    # shared filesystems (Lustre, NFSv4) this is honored by Slurm workers
    # writing to the same file.
    with open(summary_path, "a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            # If file is empty (we were the first to grab the lock), prepend header.
            f.seek(0, os.SEEK_END)
            if f.tell() == 0:
                f.write(SUMMARY_HEADER)
            f.write(line)
            f.flush()
            os.fsync(f.fileno())
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def process_root(input_root, output_h5, dx, dz, nmax,
                 events_per_chunk=100, summary_file=None):
    """
    Stream a ROOT file in chunks of `events_per_chunk` events, cluster each
    event, and write a vlen-based HDF5 where each event has exactly n_kept cells.

    If `summary_file` is given, append one summary line for this ROOT file,
    protected by fcntl.flock so parallel Slurm tasks can share the same file.
    """

    material = parse_filename_material(input_root)

    out_dir = os.path.dirname(os.path.abspath(output_h5))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    vlen   = h5py.vlen_dtype(np.float32)
    str_dt = h5py.string_dtype()

    tree = uproot.open(f"{input_root}:photon_sim")

    # 1) Read only EventID to find event boundaries.
    all_ids = tree["EventID"].array(library="np")
    n_hits = len(all_ids)

    if n_hits == 0:
        # Empty file → write schema-compatible empty file so combine step doesn't choke.
        with h5py.File(output_h5, "w") as f:
            f.create_dataset("showers", (0,), maxshape=(None,), dtype=vlen, chunks=(1024,))
            f.create_dataset("n_cells", (0,), maxshape=(None,), dtype=np.int32, chunks=(4096,))
            f.create_dataset("EventID", (0,), maxshape=(None,), dtype=np.int64, chunks=(4096,))
            f.create_dataset("primaryE", (0,), maxshape=(None,), dtype=np.float32, chunks=(4096,))
            f.create_dataset("material", (0,), maxshape=(None,), dtype=str_dt, chunks=(4096,))
            f.create_dataset("shape", data=np.array([0, 0, N_FEAT], dtype=np.int64))
            _write_attrs(f, dx, dz, nmax, max_planes=0)
        print(f"[DONE] 0 events → {output_h5}")

        _append_summary_line(summary_file,
            f"{os.path.basename(input_root):<80s} {material:>8s} "
            f"{0:>10d} {0:>14d} {0:>14d} {0:>11d} "
            f"{0:>11d} {0:>12d} "
            f"{0.0:>14.4f} {0.0:>14.4f} {0.0:>10.2f}%\n"
        )
        return

    # Event-start row indices, with a sentinel at n_hits.
    boundaries = np.flatnonzero(all_ids[1:] != all_ids[:-1]) + 1
    event_starts = np.concatenate(([0], boundaries, [n_hits])).astype(np.int64)
    n_events = len(event_starts) - 1

    del all_ids

    with h5py.File(output_h5, "w") as f:

        showers = f.create_dataset(
            "showers", (0,), maxshape=(None,),
            dtype=vlen, chunks=(1024,),
        )
        n_cells = f.create_dataset(
            "n_cells", (0,), maxshape=(None,),
            dtype=np.int32, chunks=(4096,),
        )
        event_id = f.create_dataset(
            "EventID", (0,), maxshape=(None,),
            dtype=np.int64, chunks=(4096,),
        )
        primaryE = f.create_dataset(
            "primaryE", (0,), maxshape=(None,),
            dtype=np.float32, chunks=(4096,),
        )
        material_ds = f.create_dataset(
            "material", (0,), maxshape=(None,),
            dtype=str_dt, chunks=(4096,),
        )

        datasets = (showers, n_cells, event_id, primaryE, material_ds)

        total = 0
        max_cells     = 0   # biggest n_kept across all events in this file
        max_planes    = 0   # biggest (plane_idx + 1) across all events
        max_raw_cells = 0   # biggest n_raw (before truncation) — for summary
        n_clusters_raw_total  = 0
        n_clusters_kept_total = 0
        n_events_truncated    = 0
        energy_raw_total      = 0.0
        energy_kept_total     = 0.0

        # 2) Walk events in groups.
        for ev0 in range(0, n_events, events_per_chunk):
            ev1 = min(ev0 + events_per_chunk, n_events)

            row_start = int(event_starts[ev0])
            row_stop  = int(event_starts[ev1])
            if row_stop <= row_start:
                continue

            arr = tree.arrays(
                ["EventID", "primaryE", "x", "y", "z", "dE"],
                entry_start=row_start,
                entry_stop=row_stop,
                library="np",
            )
            ids = arr["EventID"]
            E0  = arr["primaryE"]
            x, y, z, dE = arr["x"], arr["y"], arr["z"], arr["dE"]

            # Local event boundaries within this chunk.
            local_bd     = np.flatnonzero(ids[1:] != ids[:-1]) + 1
            local_splits = np.concatenate(([0], local_bd, [len(ids)]))

            b_show, b_nc, b_eid, b_E0, b_mat = [], [], [], [], []

            for s, t in zip(local_splits[:-1], local_splits[1:]):
                cl = cluster_event_per_plane(x[s:t], y[s:t], z[s:t], dE[s:t], dx, dz)
                n_raw  = cl.shape[0]
                e_raw  = float(cl[:, 3].sum()) if n_raw > 0 else 0.0

                kept = truncate_and_sort(cl, nmax)     # (n_kept, N_FEAT), n_kept ≤ nmax
                n_kept = kept.shape[0]
                e_kept = float(kept[:, 3].sum()) if n_kept > 0 else 0.0

                # Flatten to 1-D for vlen storage; reshape on read as (n_kept, N_FEAT).
                flat = kept.reshape(-1).astype(np.float32)

                b_show.append(flat)
                b_nc.append(n_kept)
                b_eid.append(ids[s])
                b_E0.append(E0[s])
                b_mat.append(material)

                # Running totals for summary
                n_clusters_raw_total  += n_raw
                n_clusters_kept_total += n_kept
                energy_raw_total      += e_raw
                energy_kept_total     += e_kept
                if n_raw > nmax:
                    n_events_truncated += 1
                if n_raw > max_raw_cells:
                    max_raw_cells = n_raw

                if n_kept > max_cells:
                    max_cells = n_kept
                if n_kept > 0:
                    mp = int(kept[:, 2].max()) + 1  # max plane_idx + 1 = #planes covered
                    if mp > max_planes:
                        max_planes = mp

            total = _flush_events(
                datasets,
                (b_show, b_nc, b_eid, b_E0, b_mat),
                total,
            )

        f.create_dataset("shape", data=np.array([total, max_cells, N_FEAT], dtype=np.int64))
        _write_attrs(f, dx, dz, nmax, max_planes)

    retention_pct = (100.0 * energy_kept_total / energy_raw_total
                     if energy_raw_total > 0 else 0.0)

    print(f"[DONE] {total} events → {output_h5}  "
          f"(dx={dx}, dz={dz}, nmax={nmax}, "
          f"max n_cells={max_cells}, max planes={max_planes}, "
          f"energy retention={retention_pct:.2f}%)")

    _append_summary_line(summary_file,
        f"{os.path.basename(input_root):<80s} {material:>8s} "
        f"{total:>10d} "
        f"{n_clusters_raw_total:>14d} {n_clusters_kept_total:>14d} "
        f"{n_events_truncated:>11d} "
        f"{max_raw_cells:>11d} {max_cells:>12d} "
        f"{energy_raw_total:>14.4f} {energy_kept_total:>14.4f} "
        f"{retention_pct:>10.2f}%\n"
    )


# =============================================================================
# Combine many per-file H5s into one
# =============================================================================

def combine_h5(files, output, chunk_size=5000):
    """
    Combine per-file H5s produced by process_root(). All files must share the
    same (dx, dz, nmax); combine enforces this by reading each file's attrs.
    """

    if not files:
        raise RuntimeError("No per-file h5 files found to combine.")

    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    vlen   = h5py.vlen_dtype(np.float32)
    str_dt = h5py.string_dtype()
    files  = sorted(files)

    # First pass: validate consistency + count events + track max n_cells/max planes.
    N          = 0
    max_cells  = 0
    max_planes = 0
    nmax = dx = dz = None

    for fpath in files:
        with h5py.File(fpath, "r") as fin:
            file_nmax = int(fin.attrs.get("nmax", -1))
            file_dx   = float(fin.attrs.get("dx", -1))
            file_dz   = float(fin.attrs.get("dz", -1))

            if nmax is None:
                nmax, dx, dz = file_nmax, file_dx, file_dz
            elif (file_nmax, file_dx, file_dz) != (nmax, dx, dz):
                raise RuntimeError(
                    f"Inconsistent clustering params across files. "
                    f"First file had (dx={dx}, dz={dz}, nmax={nmax}); "
                    f"{fpath} has (dx={file_dx}, dz={file_dz}, nmax={file_nmax})."
                )

            nc_ds = fin["n_cells"]
            n = nc_ds.shape[0]
            N += n
            for s in range(0, n, chunk_size):
                t = min(s + chunk_size, n)
                blk = nc_ds[s:t]
                if blk.size:
                    m = int(blk.max())
                    if m > max_cells:
                        max_cells = m

            fmp = int(fin.attrs.get("n_planes_max", 0))
            if fmp > max_planes:
                max_planes = fmp

    print(f"[COMBINE] {len(files)} files → {N} total events  "
          f"(dx={dx}, dz={dz}, nmax={nmax}, max n_cells seen={max_cells})")

    with h5py.File(output, "w") as f:

        showers = f.create_dataset(
            "showers", (0,), maxshape=(None,),
            dtype=vlen, chunks=(1024,),
        )
        n_cells = f.create_dataset(
            "n_cells", (0,), maxshape=(None,),
            dtype=np.int32, chunks=(4096,),
        )
        event_id = f.create_dataset(
            "EventID", (0,), maxshape=(None,),
            dtype=np.int64, chunks=(4096,),
        )
        primaryE = f.create_dataset(
            "primaryE", (0,), maxshape=(None,),
            dtype=np.float32, chunks=(4096,),
        )
        material_ds = f.create_dataset(
            "material", (0,), maxshape=(None,),
            dtype=str_dt, chunks=(4096,),
        )

        # Pre-size once.
        showers.resize((N,))
        n_cells.resize((N,))
        event_id.resize((N,))
        primaryE.resize((N,))
        material_ds.resize((N,))

        offset = 0
        for fpath in files:
            with h5py.File(fpath, "r") as fin:
                n = fin["showers"].shape[0]
                if n == 0:
                    continue

                for s in range(0, n, chunk_size):
                    t = min(s + chunk_size, n)
                    o0 = offset + s
                    o1 = offset + t

                    # vlen showers must be copied element-by-element.
                    src_showers = fin["showers"][s:t]
                    for i in range(t - s):
                        showers[o0 + i] = src_showers[i]

                    n_cells[o0:o1]     = fin["n_cells"][s:t]
                    event_id[o0:o1]    = fin["EventID"][s:t]
                    primaryE[o0:o1]    = fin["primaryE"][s:t]
                    material_ds[o0:o1] = fin["material"][s:t]

                offset += n

        f.create_dataset("shape", data=np.array([offset, max_cells, N_FEAT], dtype=np.int64))
        _write_attrs(f, dx, dz, nmax, max_planes)

    print(f"[COMBINED] {offset} events → {output}")


# =============================================================================
# Slurm submission
# =============================================================================

def submit_slurm(args):

    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise RuntimeError(f"No ROOT files matched: {args.input_glob}")

    work   = os.path.abspath(args.work_dir)
    logs   = os.path.join(work, "logs")
    sbatch = os.path.join(work, "sbatch")

    os.makedirs(work, exist_ok=True)
    os.makedirs(logs, exist_ok=True)
    os.makedirs(sbatch, exist_ok=True)
    os.makedirs(args.per_file_dir, exist_ok=True)

    # Default summary file lives next to the per-file H5s unless overridden.
    summary_file = args.summary_file or os.path.join(
        os.path.abspath(args.per_file_dir), "clustering_summary.txt"
    )
    # Start fresh so old runs' lines don't carry over.
    if os.path.exists(summary_file):
        os.remove(summary_file)

    root_list = os.path.join(work, "roots.txt")
    with open(root_list, "w") as f:
        f.write("\n".join(files))

    proc = os.path.join(sbatch, "proc.sbatch")
    comb = os.path.join(sbatch, "comb.sbatch")

    with open(proc, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --job-name=root2h5
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output={logs}/%A_%a.out
#SBATCH --error={logs}/%A_%a.err
#SBATCH --array=0-{len(files)-1}
#SBATCH --account={args.account}

set -euo pipefail

module load python/3.12
# Activate your venv here if needed, e.g.:
# source /project/ctb-stelzer/hamza95/venv/bin/activate

ROOT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" {root_list})
OUT={args.per_file_dir}/"$(basename "$ROOT" .root)".h5

python {os.path.abspath(__file__)} \\
    --input-root "$ROOT" --output "$OUT" \\
    --dx {args.dx} --dz {args.dz} --nmax {args.nmax} \\
    --events-per-chunk {args.events_per_chunk} \\
    --summary-file "{summary_file}"
""")

    with open(comb, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output={logs}/combine_%j.out
#SBATCH --error={logs}/combine_%j.err
#SBATCH --account={args.account}

set -euo pipefail

module load python/3.12
# source /project/ctb-stelzer/hamza95/venv/bin/activate

python {os.path.abspath(__file__)} --combine-only \\
    --input-glob "{args.per_file_dir}/*.h5" \\
    --output "{args.output}" \\
    --combine-chunk-size {args.combine_chunk_size}
""")

    jid = subprocess.check_output(["sbatch", proc]).decode().split()[-1]
    subprocess.run(["sbatch", f"--dependency=afterok:{jid}", comb], check=True)

    print(f"Submitted array job: {jid}  ({len(files)} tasks)")
    print(f"Logs:     {logs}")
    print(f"Sbatch:   {sbatch}")
    print(f"Per-file: {args.per_file_dir}")
    print(f"Summary:  {summary_file}")
    print(f"Params:   dx={args.dx} cm, dz={args.dz} cm, nmax={args.nmax}")


# =============================================================================
# CLI
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="ROOT → per-plane clustered HDF5 with vlen (variable-length) storage."
    )
    p.add_argument("--input-root")
    p.add_argument("--input-glob")
    p.add_argument("--output", required=True)
    p.add_argument("--combine-only", action="store_true")
    p.add_argument("--submit-slurm", action="store_true")
    p.add_argument("--work-dir",     default="./slurm_work")
    p.add_argument("--per-file-dir", default="./per_file_h5")
    p.add_argument("--account",      default="def-mdanning")

    # Clustering parameters.
    p.add_argument("--dx",   type=float, default=DEFAULT_DX,
                   help=f"Transverse cell size in cm, also used for dy (default: {DEFAULT_DX})")
    p.add_argument("--dz",   type=float, default=DEFAULT_DZ,
                   help=f"Plane (z) slab thickness in cm (default: {DEFAULT_DZ})")
    p.add_argument("--nmax", type=int, default=DEFAULT_NMAX,
                   help=f"Max cells per event (top-N by energy, default: {DEFAULT_NMAX})")

    p.add_argument("--events-per-chunk",   type=int, default=100,
                   help="Events per streaming chunk in process_root (default: 100)")
    p.add_argument("--combine-chunk-size", type=int, default=5000,
                   help="Events per streaming window in combine_h5 (default: 5000)")
    p.add_argument("--summary-file",       default=None,
                   help="Optional path to a shared text summary file. Each processed "
                        "ROOT file appends one line (lock-protected with fcntl.flock). "
                        "Defaults to <per-file-dir>/clustering_summary.txt when using "
                        "--submit-slurm.")

    args = p.parse_args()

    if args.submit_slurm:
        if not args.input_glob:
            p.error("--submit-slurm requires --input-glob")
        submit_slurm(args)
    elif args.combine_only:
        if not args.input_glob:
            p.error("--combine-only requires --input-glob")
        combine_h5(glob.glob(args.input_glob), args.output,
                   chunk_size=args.combine_chunk_size)
    else:
        if not args.input_root:
            p.error("processing mode requires --input-root")
        process_root(args.input_root, args.output,
                     dx=args.dx, dz=args.dz, nmax=args.nmax,
                     events_per_chunk=args.events_per_chunk,
                     summary_file=args.summary_file)


if __name__ == "__main__":
    main()