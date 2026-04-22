#!/usr/bin/env bash
#SBATCH --job-name=combine
#SBATCH --time=10:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=/project/ctb-stelzer/hamza95/photons_gen/logs/log_%x_%A_%a.out
#SBATCH --error=/project/ctb-stelzer/hamza95/photons_gen/logs/log_%x_%A_%a.err
#SBATCH --account=def-mdanning
#SBATCH --cpus-per-task=8

module load python/3.12

python /project/ctb-stelzer/hamza95/photons_gen/create_h5_file.py \
  --input-glob "/project/ctb-stelzer/hamza95/photo_gen_files_*/*.root" \
  --tree-name photon_sim \
  --output /scratch/hamza95/all_clustered_showers.h5