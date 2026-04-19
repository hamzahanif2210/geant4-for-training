#!/usr/bin/env python3
"""
Generate SLURM bash scripts for Geant4 calorimeter simulation jobs.

Cell-size families and energies follow the paper dataset:
  - 5 cell configurations: 1x1x5, 2x2x4, 3x3x8, 4x4x10, 5x5x15 cm^3
  - 11 photon energies: 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 GeV
  - ~9091 events per job -> ~100k events per cell configuration
"""

import argparse
import os
import math

# Default dataset parameters from the paper
CELL_CONFIGS = [
    (1, 1, 5),
    (2, 2, 4),
    (3, 3, 8),
    (4, 4, 10),
    (5, 5, 15),
]

ENERGIES_GEV = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

EVENTS_PER_CONFIG = 100_000


def make_macro(energy_gev, n_events, n_threads, macro_path):
    lines = [
        f"/run/numberOfThreads {n_threads}",
        "/run/initialize",
        f"/gun/energy {energy_gev} GeV",
        f"/run/beamOn {n_events}",
    ]
    with open(macro_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_slurm_script(
    cx, cy, cz, energy_gev, n_events, macro_path, output_file,
    exe_path, job_dir, n_threads, time_limit, mem_gb, account, partition
):
    tag = f"cx{cx}_cy{cy}_cz{cz}_E{energy_gev}GeV"
    script_path = os.path.join(job_dir, f"job_{tag}.sh")
    log_path = os.path.join(job_dir, f"log_{tag}.txt")

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
echo "Cell size: {cx} x {cy} x {cz} cm^3,  Energy: {energy_gev} GeV,  Events: {n_events}"
date

{exe_path} \\
    -m {macro_path} \\
    -cx {cx} -cy {cy} -cz {cz} \\
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
    parser.add_argument(
        "--exe", required=True,
        help="Path to the compiled exampleB4a executable"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory where output .root files will be written"
    )
    parser.add_argument(
        "--job-dir", default="jobs",
        help="Directory to write generated job scripts and macros (default: jobs/)"
    )
    parser.add_argument(
        "--events-per-config", type=int, default=EVENTS_PER_CONFIG,
        help=f"Total events per cell configuration (default: {EVENTS_PER_CONFIG})"
    )
    parser.add_argument(
        "--threads", type=int, default=8,
        help="Number of CPU threads per job (default: 8)"
    )
    parser.add_argument(
        "--time", default="02:00:00",
        help="SLURM time limit per job, e.g. 02:00:00 (default: 02:00:00)"
    )
    parser.add_argument(
        "--mem", type=int, default=8,
        help="Memory in GB per job (default: 8)"
    )
    parser.add_argument(
        "--account", default="",
        help="SLURM account (e.g. ctb-stelzer)"
    )
    parser.add_argument(
        "--partition", default="",
        help="SLURM partition (e.g. cpu)"
    )
    parser.add_argument(
        "--cell-configs", nargs="+", metavar="CX,CY,CZ",
        help="Override cell configurations as CX,CY,CZ in cm (e.g. 1,1,5 2,2,4). "
             "Default: all five paper configurations."
    )
    parser.add_argument(
        "--energies", nargs="+", type=float, metavar="E_GEV",
        default=ENERGIES_GEV,
        help="Photon energies in GeV (default: 1 10 20 30 40 50 60 70 80 90 100)"
    )
    args = parser.parse_args()

    # Parse cell configs
    if args.cell_configs:
        cell_configs = []
        for spec in args.cell_configs:
            parts = spec.split(",")
            if len(parts) != 3:
                parser.error(f"Cell config must be CX,CY,CZ — got: {spec!r}")
            cell_configs.append(tuple(float(p) for p in parts))
    else:
        cell_configs = CELL_CONFIGS

    energies = args.energies
    events_per_energy = math.ceil(args.events_per_config / len(energies))

    os.makedirs(args.job_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    macro_dir = os.path.join(args.job_dir, "macros")
    os.makedirs(macro_dir, exist_ok=True)

    all_scripts = []

    for cx, cy, cz in cell_configs:
        for energy in energies:
            tag = f"cx{cx}_cy{cy}_cz{cz}_E{energy}GeV"

            # Format sizes: strip trailing .0 for clean filenames
            def fmt(v):
                return str(int(v)) if v == int(v) else str(v)

            macro_path = os.path.join(macro_dir, f"run_{tag}.mac")
            output_file = os.path.join(
                args.output_dir,
                f"photons_{fmt(cx)}x{fmt(cy)}x{fmt(cz)}cm_{fmt(energy)}GeV.root"
            )

            make_macro(
                energy_gev=energy,
                n_events=events_per_energy,
                n_threads=args.threads,
                macro_path=macro_path,
            )

            script_path = make_slurm_script(
                cx=cx, cy=cy, cz=cz,
                energy_gev=energy,
                n_events=events_per_energy,
                macro_path=macro_path,
                output_file=output_file,
                exe_path=args.exe,
                job_dir=args.job_dir,
                n_threads=args.threads,
                time_limit=args.time,
                mem_gb=args.mem,
                account=args.account,
                partition=args.partition,
            )
            all_scripts.append(script_path)

    # Master submission script
    submit_script = os.path.join(args.job_dir, "submit_all.sh")
    with open(submit_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Submit all generated Geant4 simulation jobs\n")
        f.write(f"# Total jobs: {len(all_scripts)}\n\n")
        for s in all_scripts:
            f.write(f"sbatch {s}\n")
    os.chmod(submit_script, 0o755)

    print(f"Generated {len(all_scripts)} job scripts in: {args.job_dir}/")
    print(f"  Macros:   {macro_dir}/")
    print(f"  Outputs:  {args.output_dir}/")
    print(f"  Submit all with:  bash {submit_script}")
    print(f"\nJob breakdown:")
    print(f"  Cell configs : {len(cell_configs)}")
    print(f"  Energies     : {len(energies)}")
    print(f"  Events/job   : {events_per_energy}")
    print(f"  Events/config: ~{events_per_energy * len(energies)}")


if __name__ == "__main__":
    main()
