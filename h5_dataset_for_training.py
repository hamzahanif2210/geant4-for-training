#!/usr/bin/env python3
"""
Build a training subset (and optionally a disjoint test subset) from a combined h5.

Selects a random N events per material from a combined h5 produced by
create_h5_file.py and writes a new h5 with the same schema. If test-set
arguments are given, also writes a second h5 drawn from the events NOT
selected for training — so the two sets are guaranteed disjoint.

Schema (matches create_h5_file.py):
  Datasets:
    showers     vlen float32, each entry is (n_cells * 4,) — reshape to (n_cells, 4)
                columns: x, y, plane_idx, energy
    n_cells     int32   (N,)     actual cells per event
    EventID     int64   (N,)
    primaryE    float32 (N,)
    material    str     (N,)     "PbF2" or "PbWO4"
    shape       int64   (3,)     [N, max_cells, 4]
  Attributes on the root group:
    dx, dy, dz           cm — transverse/longitudinal cell sizes
    nmax                 top-N-by-energy cap used during clustering
    n_planes_max         highest plane_idx + 1 seen in this file
    feature_columns      "x,y,plane_idx,energy"

Usage (train only):
    python h5_dataset_for_training.py \
        --input  /scratch/hamza95/all_clustered_showers.h5 \
        --output /scratch/hamza95/train_130k.h5 \
        --per-material PbF2=65000 PbWO4=65000 \
        --seed 42

Usage (train + test, equal counts):
    python /project/ctb-stelzer/hamza95/photons_gen/h5_dataset_for_training.py \
        --input  /scratch/hamza95/all_clustered_showers.h5 \
        --output /scratch/hamza95/train_130k.h5 \
        --n-per-material 65000 \
        --test-output /scratch/hamza95/test_50k.h5 \
        --n-per-material-test 25000 \
        --seed 42

Usage (train + test, explicit per-material counts for each):
    python h5_dataset_for_training.py \
        --input  all_clustered_showers.h5 \
        --output train.h5      --per-material      PbF2=65000 PbWO4=65000 \
        --test-output test.h5  --per-material-test PbF2=15000 PbWO4=15000 \
        --seed 42
"""

import os
import argparse

import h5py
import numpy as np


# Attributes we expect on input files and propagate verbatim to output files.
# Keeping this as a named list so adding new per-file metadata later is a one-line change.
_PROPAGATED_ATTRS = (
    "dx", "dy", "dz", "nmax", "n_planes_max", "feature_columns",
)


def parse_per_material(items):
    """Parse ['PbF2=65000', 'PbWO4=65000'] into {'PbF2': 65000, 'PbWO4': 65000}."""
    out = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Bad entry: {item!r} (expected MATERIAL=COUNT)")
        mat, n = item.split("=", 1)
        out[mat] = int(n)
    return out


def load_materials(input_h5, chunk_size=100_000):
    """Read the material dataset in chunks, return as a numpy array of str."""
    with h5py.File(input_h5, "r") as f:
        ds = f["material"]
        n = ds.shape[0]

        # Read in chunks so we don't balloon memory on huge files. The material
        # column itself is tiny (one short string per event) so this is cheap.
        out = np.empty(n, dtype=object)
        for s in range(0, n, chunk_size):
            t = min(s + chunk_size, n)
            blk = ds[s:t]
            # h5py returns bytes for string_dtype datasets — decode.
            if n > 0 and isinstance(blk[0], (bytes, bytearray)):
                blk = np.array([b.decode("utf-8") for b in blk], dtype=object)
            out[s:t] = blk
        return out


def select_indices(materials, per_material, rng, allow_short=False,
                   exclude=None, label="SELECT"):
    """
    Return a dict {material: sorted np.int64 array of selected indices}.

    If `exclude` is provided, it must be a dict {material: 1-D int array} of
    indices to remove from the candidate pool before sampling — used to make
    the test set disjoint from the train set.
    """
    exclude = exclude or {}
    selected = {}

    for mat, n_want in per_material.items():
        all_idx = np.flatnonzero(materials == mat)

        if mat in exclude and len(exclude[mat]) > 0:
            # np.setdiff1d is fine here — both inputs are int arrays and
            # len(all_idx) is bounded by total events per material.
            pool = np.setdiff1d(all_idx, exclude[mat], assume_unique=True)
        else:
            pool = all_idx

        n_avail = len(pool)

        if n_avail == 0:
            msg = f"Material {mat!r} has 0 events available ({label})."
            if allow_short:
                print(f"[WARN] {msg} Skipping.")
                continue
            raise RuntimeError(msg)

        if n_avail < n_want:
            msg = (f"Material {mat!r} ({label}): requested {n_want} events "
                   f"but only {n_avail} available.")
            if allow_short:
                print(f"[WARN] {msg} Taking all {n_avail}.")
                take = pool
            else:
                raise RuntimeError(msg + " Use --allow-short to take what's available.")
        else:
            take = rng.choice(pool, size=n_want, replace=False)

        # HDF5 fancy indexing requires sorted, unique indices for efficient reads.
        take = np.sort(take.astype(np.int64))
        selected[mat] = take
        print(f"[{label}] {mat}: {len(take)} / {n_avail} available")

    return selected


def write_subset(input_h5, output_h5, selected, rng, read_chunk=5000, tag="OUT"):
    """
    Stream selected events from input_h5 into output_h5, shuffling the final
    output order so materials are interleaved. Propagates clustering attributes
    (dx, dy, dz, nmax, ...) from input to output.
    """

    # Build the final (shuffled) output order and the corresponding input
    # indices to read from. We copy material-by-material (contiguous within
    # the input file) for cache-friendly reads, then write into shuffled
    # output positions.
    #
    # Plan: assemble one flat array `input_idx` of shape (N_total,) containing
    # the input row indices we want, in input-sorted order per material, then
    # generate `output_pos` = a shuffled permutation telling us where each
    # input row lands in the output. Reads stay sorted for speed; the shuffle
    # happens via the output_pos indirection.

    mats_in_order = list(selected.keys())
    per_mat_idx = [selected[m] for m in mats_in_order]
    per_mat_sizes = [len(a) for a in per_mat_idx]
    n_total = sum(per_mat_sizes)

    if n_total == 0:
        raise RuntimeError(f"Nothing to write for {tag} — aborting.")

    # output_pos[i] = destination row in the output file for input row
    # input_idx[i]. A random permutation of 0..N-1 gives us material shuffling.
    output_pos = rng.permutation(n_total)

    # Walk through per-material blocks with a running input-offset cursor.
    cursors = np.cumsum([0] + per_mat_sizes)

    vlen = h5py.vlen_dtype(np.float32)
    str_dt = h5py.string_dtype()

    out_dir = os.path.dirname(os.path.abspath(output_h5))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with h5py.File(input_h5, "r") as fin, h5py.File(output_h5, "w") as fout:

        # Pre-size output datasets to N_total. Schema matches create_h5_file.py:
        # showers (vlen), n_cells, EventID, primaryE, material  — no cellsize anymore.
        showers_o = fout.create_dataset("showers", (n_total,), maxshape=(None,),
                                        dtype=vlen, chunks=(1024,))
        n_cells_o = fout.create_dataset("n_cells", (n_total,), maxshape=(None,),
                                        dtype=np.int32, chunks=(4096,))
        event_id_o = fout.create_dataset("EventID", (n_total,), maxshape=(None,),
                                         dtype=np.int64, chunks=(4096,))
        primaryE_o = fout.create_dataset("primaryE", (n_total,), maxshape=(None,),
                                         dtype=np.float32, chunks=(4096,))
        material_o = fout.create_dataset("material", (n_total,), maxshape=(None,),
                                         dtype=str_dt, chunks=(4096,))

        showers_i  = fin["showers"]
        n_cells_i  = fin["n_cells"]
        event_id_i = fin["EventID"]
        primaryE_i = fin["primaryE"]
        material_i = fin["material"]

        max_cells = 0
        max_planes_seen = 0  # track plane_idx max for the output attr

        for mi, mat in enumerate(mats_in_order):
            idx = per_mat_idx[mi]               # sorted input indices for this material
            base = cursors[mi]                  # position within the flat input_idx array

            if len(idx) == 0:
                continue

            # Stream this material's selected rows in read_chunk-sized windows.
            for s in range(0, len(idx), read_chunk):
                t = min(s + read_chunk, len(idx))
                sub = idx[s:t]                  # still sorted, still unique

                # Corresponding destination rows in the output file.
                dst = output_pos[base + s : base + t]

                # Fixed-size branches: bulk read with fancy indexing.
                nc_block  = n_cells_i[sub]
                eid_block = event_id_i[sub]
                E0_block  = primaryE_i[sub]
                mat_block = material_i[sub]

                # vlen showers: must read element-by-element (h5py limitation).
                show_block = [showers_i[int(k)] for k in sub]

                # Writes: dst is NOT sorted (shuffled). HDF5 supports writing
                # to an arbitrary index list, so this is fine — it's just
                # slower than sorted writes, which is unavoidable if we want
                # a shuffled output.
                for j, d in enumerate(dst):
                    showers_o[int(d)]  = show_block[j]
                    n_cells_o[int(d)]  = nc_block[j]
                    event_id_o[int(d)] = eid_block[j]
                    primaryE_o[int(d)] = E0_block[j]
                    material_o[int(d)] = mat_block[j]

                if nc_block.size:
                    m = int(nc_block.max())
                    if m > max_cells:
                        max_cells = m

                    # Update plane_idx tracking from the shower payloads.
                    # Each vlen entry is a flat (n_cells * 4,) array; col 2 is plane_idx.
                    for flat, nc in zip(show_block, nc_block):
                        if nc > 0:
                            arr = np.asarray(flat, dtype=np.float32).reshape(int(nc), 4)
                            mp = int(arr[:, 2].max()) + 1
                            if mp > max_planes_seen:
                                max_planes_seen = mp

                print(f"[{tag}] {mat}: {t}/{len(idx)}")

        fout.create_dataset("shape", data=np.array([n_total, max_cells, 4], dtype=np.int64))

        # Propagate clustering attributes from input → output.
        for k in _PROPAGATED_ATTRS:
            if k in fin.attrs:
                fout.attrs[k] = fin.attrs[k]
        # Override n_planes_max with what we actually saw in the subset
        # (the full-file value from the input may be larger).
        if max_planes_seen > 0:
            fout.attrs["n_planes_max"] = int(max_planes_seen)

    print(f"[DONE {tag}] {n_total} events → {output_h5}")


def resolve_per_material(per_material_arg, n_per_material_arg, materials, kind):
    """
    Return {material: count} from either --(...)per-material or
    --n-(...)per-material. `kind` is a human label for error messages,
    e.g. 'train' or 'test'.
    """
    if per_material_arg is None and n_per_material_arg is None:
        return None
    if per_material_arg is not None and n_per_material_arg is not None:
        raise ValueError(
            f"Specify only one of --per-material{'-test' if kind == 'test' else ''}"
            f" or --n-per-material{'-test' if kind == 'test' else ''} (got both for {kind})."
        )

    if n_per_material_arg is not None:
        unique = sorted(set(materials.tolist()))
        out = {m: n_per_material_arg for m in unique}
        print(f"[INFO] {kind}: --n-per-material{'-test' if kind == 'test' else ''} "
              f"{n_per_material_arg} applied to: {unique}")
        return out

    return parse_per_material(per_material_arg)


def main():
    p = argparse.ArgumentParser(
        description="Build a random per-material training subset (and optional "
                    "disjoint test subset) from a combined h5."
    )
    p.add_argument("--input", required=True,
                   help="Combined h5 file produced by create_h5_file.py")

    # Train set
    p.add_argument("--output", required=True,
                   help="Output h5 training subset")
    p.add_argument("--per-material", nargs="+", default=None,
                   help="Per-material train counts, e.g. PbF2=65000 PbWO4=65000")
    p.add_argument("--n-per-material", type=int, default=None,
                   help="Shortcut: same train count for every material present")

    # Test set (optional, disjoint from train)
    p.add_argument("--test-output", default=None,
                   help="Optional output h5 test subset (drawn from events NOT "
                        "selected for training).")
    p.add_argument("--per-material-test", nargs="+", default=None,
                   help="Per-material test counts, e.g. PbF2=15000 PbWO4=15000")
    p.add_argument("--n-per-material-test", type=int, default=None,
                   help="Shortcut: same test count for every material present")

    # Common
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility (default: 42)")
    p.add_argument("--allow-short", action="store_true",
                   help="If a material has fewer events than requested, take all "
                        "available and continue (instead of erroring).")
    p.add_argument("--read-chunk", type=int, default=5000,
                   help="Events per read/write window when streaming (default: 5000)")

    args = p.parse_args()

    # Train spec must be present.
    if (args.per_material is None) == (args.n_per_material is None):
        p.error("Specify exactly one of --per-material or --n-per-material.")

    # Test: spec and output path must be either both set or both absent.
    test_spec_given = (args.per_material_test is not None) or (args.n_per_material_test is not None)
    if test_spec_given and args.test_output is None:
        p.error("--per-material-test / --n-per-material-test requires --test-output.")
    if args.test_output is not None and not test_spec_given:
        p.error("--test-output requires --per-material-test or --n-per-material-test.")

    rng = np.random.default_rng(args.seed)

    print(f"[READ] materials from {args.input}")
    materials = load_materials(args.input)
    print(f"[READ] {len(materials)} total events")

    # --- Train selection ---
    train_counts = resolve_per_material(
        args.per_material, args.n_per_material, materials, kind="train"
    )
    train_selected = select_indices(
        materials, train_counts, rng,
        allow_short=args.allow_short, label="TRAIN"
    )
    if not train_selected:
        raise RuntimeError("Nothing selected for training — aborting.")

    write_subset(args.input, args.output, train_selected, rng,
                 read_chunk=args.read_chunk, tag="TRAIN")

    # --- Test selection (disjoint from train) ---
    if test_spec_given:
        test_counts = resolve_per_material(
            args.per_material_test, args.n_per_material_test, materials, kind="test"
        )
        test_selected = select_indices(
            materials, test_counts, rng,
            allow_short=args.allow_short,
            exclude=train_selected,          # <-- guarantees disjointness
            label="TEST",
        )
        if not test_selected:
            raise RuntimeError("Nothing selected for test — aborting.")

        write_subset(args.input, args.test_output, test_selected, rng,
                     read_chunk=args.read_chunk, tag="TEST")


if __name__ == "__main__":
    main()