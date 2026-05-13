#!/bin/bash
#SBATCH --job-name=fsqc_cohort
#SBATCH --partition=cpu20
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/fsqc/slurm_logs/%x_%j.out
#SBATCH --error=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/fsqc/slurm_logs/%x_%j.err

source ~/.bashrc
micromamba activate fsqc

COHORT="$1"
SUFFIX="$2"
echo "[fsqc] Running cohort '$COHORT' with suffix '$SUFFIX'"

/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/scripts/run_fsqc_on_cohort.sh "$COHORT" "$SUFFIX"
