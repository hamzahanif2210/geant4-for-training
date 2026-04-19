#!/usr/bin/env python3

'''
python /project/ctb-stelzer/hamza95/photons_gen/generate_jobs.py \
  --exe /project/ctb-stelzer/hamza95/photons_gen/build/exampleB4a \
  --output-dir /project/ctb-stelzer/hamza95/photo_gen_files \
  --job-dir /project/ctb-stelzer/hamza95/photo_jobs_sbatch_files \
  --total-events 100000 \
  --threads 1 \
  --time 01:00:00 \
  --mem 4 \
  --account def-mdanning \
  --material PbF2 \
  --emin 1 \
  --emax 5
'''
import argparse
import os
import random

CELL_CONFIGS = [
    (1, 1, 5),
    (2, 2, 4),
    (3, 3, 8),
    (4, 4, 10),
    (5, 5, 15),
]

DEFAULT_EMIN = 1.0   # GeV
DEFAULT_EMAX = 5.0   # GeV
TOTAL_EVENTS = 100_000


def fmt_value(v):
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def make_macro(n_events, n_threads, macro_path):
    lines = [
        f"/run/numberOfThreads {n_threads}",
        "/run/initialize",
        "/gun/particle gamma",
        f"/run/beamOn {n_events}",
    ]
    with open(macro_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_slurm_script(
    cx, cy, cz, material, emin, emax, n_events, macro_path, output_file,
    exe_path, sbatch_dir, log_dir, n_threads, time_limit, mem_gb, account, partition
):
    tag = f"cx{fmt_value(cx)}_cy{fmt_value(cy)}_cz{fmt_value(cz)}_{material}_E{fmt_value(emin)}to{fmt_value(emax)}GeV"
    script_path = os.path.join(sbatch_dir, f"job_{tag}.sh")
    log_path = os.path.join(log_dir, f"log_{tag}.txt")

    script = f"""#!/bin/bash
#SBATCH --job-name=g4_{tag}
#SBATCH --output={log_path}
#SBATCH --time={time_limit}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={n_threads}
#SBATCH --mem={mem_gb}G
"""
    if account:
        script += f"#SBATCH --account={account}\n"
    if partition:
        script += f"#SBATCH --partition={partition}\n"

    script += f"""
set -e
echo "Starting job: {tag}"
echo "Cell size: {cx} x {cy} x {cz} cm^3, Material: {material}, Energy: {emin}-{emax} GeV, Events: {n_events}"
date

{exe_path} \\
    -m {macro_path} \\
    -cx {cx} -cy {cy} -cz {cz} \\
    -mat {material} \\
    -emin {emin} -emax {emax} \\
    -o {output_file} \\
    -t {n_threads}

echo "Done."
date
"""
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate SLURM job scripts for Geant4 calorimeter simulation."
    )
    parser.add_argument("--exe", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--job-dir", default="jobs")
    parser.add_argument("--total-events", type=int, default=TOTAL_EVENTS)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--time", default="02:00:00")
    parser.add_argument("--mem", type=int, default=8)
    parser.add_argument("--account", default="")
    parser.add_argument("--partition", default="")
    parser.add_argument("--material", default="PbF2",
                        choices=["PbF2", "PbWO4"],
                        help="Detector material (default: PbF2)")
    parser.add_argument("--emin", type=float, default=DEFAULT_EMIN,
                        metavar="E_GEV",
                        help="Minimum incident energy in GeV (default: 1)")
    parser.add_argument("--emax", type=float, default=DEFAULT_EMAX,
                        metavar="E_GEV",
                        help="Maximum incident energy in GeV (default: 5)")
    parser.add_argument("--cell-configs", nargs="+", metavar="CX,CY,CZ")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.emin >= args.emax:
        parser.error(f"--emin ({args.emin}) must be less than --emax ({args.emax})")

    if args.cell_configs:
        cell_configs = []
        for spec in args.cell_configs:
            parts = spec.split(",")
            if len(parts) != 3:
                parser.error(f"Cell config must be CX,CY,CZ — got: {spec!r}")
            try:
                cell_configs.append(tuple(float(p) for p in parts))
            except ValueError:
                parser.error(f"Invalid numeric cell config: {spec!r}")
    else:
        cell_configs = CELL_CONFIGS

    os.makedirs(args.job_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    sbatch_dir = os.path.join(args.job_dir, "sbatch_files")
    log_dir = os.path.join(args.job_dir, "logs")
    macro_dir = os.path.join(args.job_dir, "macros")

    os.makedirs(sbatch_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(macro_dir, exist_ok=True)

    jobs = list(cell_configs)

    rng = random.Random(args.seed)
    rng.shuffle(jobs)

    total_jobs = len(jobs)
    base_events = args.total_events // total_jobs
    leftover = args.total_events % total_jobs

    job_event_counts = {
        job: base_events + (1 if i < leftover else 0)
        for i, job in enumerate(jobs)
    }

    all_scripts = []

    for cx, cy, cz in cell_configs:
        n_events = job_event_counts[(cx, cy, cz)]

        tag = f"cx{fmt_value(cx)}_cy{fmt_value(cy)}_cz{fmt_value(cz)}_{args.material}_E{fmt_value(args.emin)}to{fmt_value(args.emax)}GeV"
        macro_path = os.path.join(macro_dir, f"run_{tag}.mac")
        output_file = os.path.join(
            args.output_dir,
            f"photons_{fmt_value(cx)}x{fmt_value(cy)}x{fmt_value(cz)}cm_"
            f"{fmt_value(args.emin)}to{fmt_value(args.emax)}GeV_{args.material}.root"
        )

        make_macro(
            n_events=n_events,
            n_threads=args.threads,
            macro_path=macro_path,
        )

        script_path = make_slurm_script(
            cx=cx, cy=cy, cz=cz,
            material=args.material,
            emin=args.emin,
            emax=args.emax,
            n_events=n_events,
            macro_path=macro_path,
            output_file=output_file,
            exe_path=args.exe,
            sbatch_dir=sbatch_dir,
            log_dir=log_dir,
            n_threads=args.threads,
            time_limit=args.time,
            mem_gb=args.mem,
            account=args.account,
            partition=args.partition,
        )
        all_scripts.append(script_path)

    submit_script = os.path.join(args.job_dir, "submit_all.sh")
    with open(submit_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Submit all generated Geant4 simulation jobs\n")
        f.write(f"# Total jobs: {len(all_scripts)}\n\n")
        for s in all_scripts:
            f.write(f"sbatch {s}\n")
    os.chmod(submit_script, 0o755)

    print(f"Generated {len(all_scripts)} job scripts")
    print(f"  Material:     {args.material}")
    print(f"  Energy range: {args.emin} - {args.emax} GeV (uniform random per event)")
    print(f"  Sbatch files: {sbatch_dir}/")
    print(f"  Logs:         {log_dir}/")
    print(f"  Macros:       {macro_dir}/")
    print(f"  Outputs:      {args.output_dir}/")
    print(f"  Submit all:   bash {submit_script}")


if __name__ == "__main__":
    main()
