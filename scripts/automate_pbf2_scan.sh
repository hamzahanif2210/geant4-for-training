#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/automate_pbf2_scan.sh [options]

Automate PbF2 scans for 5 cell-size families and 11 energies.
Total events are split as evenly as possible across all 55 scan points.

Options:
  --exec PATH          Path to Geant4 executable (default: ./exampleB4a)
  --mode MODE          local or slurm (default: local)
  --output-dir DIR     Output directory (default: scan_outputs)
  --total-events N     Total events across all runs (default: 100000)
  --seed-base N        Base random seed offset (default: 1000)
  --job-time TIME      Slurm time limit (default: 00:20:00)
  --partition NAME     Slurm partition/queue (optional)
  --threads N          Number of threads passed to -t (optional)
  -h, --help           Show this help message

Examples:
  scripts/automate_pbf2_scan.sh --exec ./build/exampleB4a
  scripts/automate_pbf2_scan.sh --mode slurm --partition short --exec ./build/exampleB4a
USAGE
}

EXECUTABLE="./exampleB4a"
MODE="local"
OUTPUT_DIR="scan_outputs"
TOTAL_EVENTS=100000
SEED_BASE=1000
JOB_TIME="00:20:00"
PARTITION=""
THREADS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exec)
      EXECUTABLE="$2"
      shift 2
      ;;
    --mode)
      MODE="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --total-events)
      TOTAL_EVENTS="$2"
      shift 2
      ;;
    --seed-base)
      SEED_BASE="$2"
      shift 2
      ;;
    --job-time)
      JOB_TIME="$2"
      shift 2
      ;;
    --partition)
      PARTITION="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ "$MODE" != "local" && "$MODE" != "slurm" ]]; then
  echo "--mode must be 'local' or 'slurm'" >&2
  exit 1
fi

if [[ ! -x "$EXECUTABLE" ]]; then
  echo "Executable not found or not executable: $EXECUTABLE" >&2
  exit 1
fi

if ! [[ "$TOTAL_EVENTS" =~ ^[0-9]+$ ]] || (( TOTAL_EVENTS <= 0 )); then
  echo "--total-events must be a positive integer" >&2
  exit 1
fi

CELL_FAMILIES=(
  "1 1 5"
  "2 2 4"
  "3 3 8"
  "4 4 10"
  "5 5 15"
)
ENERGIES_GEV=(1 10 20 30 40 50 60 70 80 90 100)

N_FAMILIES=${#CELL_FAMILIES[@]}
N_ENERGIES=${#ENERGIES_GEV[@]}
TOTAL_POINTS=$((N_FAMILIES * N_ENERGIES))
BASE_EVENTS=$((TOTAL_EVENTS / TOTAL_POINTS))
REMAINDER=$((TOTAL_EVENTS % TOTAL_POINTS))

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/macros" "$OUTPUT_DIR/logs" "$OUTPUT_DIR/results"

echo "=== PbF2 scan configuration ==="
echo "Mode:             $MODE"
echo "Executable:       $EXECUTABLE"
echo "Output directory: $OUTPUT_DIR"
echo "Families x energies: ${N_FAMILIES} x ${N_ENERGIES} = ${TOTAL_POINTS} runs"
echo "Total events:     $TOTAL_EVENTS"
echo "Events/run:       $BASE_EVENTS (+1 for first $REMAINDER runs)"
echo

run_index=0
for family in "${CELL_FAMILIES[@]}"; do
  read -r cell_x_cm cell_y_cm cell_z_cm <<< "$family"

  for energy in "${ENERGIES_GEV[@]}"; do
    run_index=$((run_index + 1))

    events_this_run=$BASE_EVENTS
    if (( run_index <= REMAINDER )); then
      events_this_run=$((events_this_run + 1))
    fi

    run_name="PbF2_${cell_x_cm}x${cell_y_cm}x${cell_z_cm}cm_E${energy}GeV"
    macro_file="$OUTPUT_DIR/macros/${run_name}.mac"
    log_file="$OUTPUT_DIR/logs/${run_name}.log"
    out_root="$OUTPUT_DIR/results/${run_name}.root"
    seed=$((SEED_BASE + run_index))

    cat > "$macro_file" <<MACRO
/run/initialize
/gun/particle gamma
/gun/energy ${energy} GeV
/analysis/setFileName ${out_root}
/run/printProgress 1000
/run/beamOn ${events_this_run}
MACRO

    if [[ "$MODE" == "local" ]]; then
      echo "[${run_index}/${TOTAL_POINTS}] Running $run_name (${events_this_run} events)"
      export B4_CELL_SIZE_X_CM="$cell_x_cm"
      export B4_CELL_SIZE_Y_CM="$cell_y_cm"
      export B4_CELL_SIZE_Z_CM="$cell_z_cm"
      export B4_NUM_CELLS_X=1
      export B4_NUM_CELLS_Y=1
      export B4_NUM_CELLS_Z=1

      if [[ -n "$THREADS" ]]; then
        "$EXECUTABLE" -m "$macro_file" -s "$seed" -t "$THREADS" > "$log_file" 2>&1
      else
        "$EXECUTABLE" -m "$macro_file" -s "$seed" > "$log_file" 2>&1
      fi
    else
      sbatch_file="$OUTPUT_DIR/${run_name}.sbatch"
      {
        echo "#!/usr/bin/env bash"
        echo "#SBATCH --job-name=$run_name"
        echo "#SBATCH --time=$JOB_TIME"
        echo "#SBATCH --output=$log_file"
        if [[ -n "$PARTITION" ]]; then
          echo "#SBATCH --partition=$PARTITION"
        fi
        echo ""
        echo "set -euo pipefail"
        echo "export B4_CELL_SIZE_X_CM=$cell_x_cm"
        echo "export B4_CELL_SIZE_Y_CM=$cell_y_cm"
        echo "export B4_CELL_SIZE_Z_CM=$cell_z_cm"
        echo "export B4_NUM_CELLS_X=1"
        echo "export B4_NUM_CELLS_Y=1"
        echo "export B4_NUM_CELLS_Z=1"
        if [[ -n "$THREADS" ]]; then
          echo "\"$EXECUTABLE\" -m \"$macro_file\" -s \"$seed\" -t \"$THREADS\""
        else
          echo "\"$EXECUTABLE\" -m \"$macro_file\" -s \"$seed\""
        fi
      } > "$sbatch_file"

      sbatch "$sbatch_file"
      echo "[${run_index}/${TOTAL_POINTS}] Submitted $run_name (${events_this_run} events)"
    fi
  done
done

echo "Completed setup for ${TOTAL_POINTS} runs."
