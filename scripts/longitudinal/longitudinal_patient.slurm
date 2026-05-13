#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --output=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/slurm_logs/%x_%j.out
#SBATCH --error=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/slurm_logs/%x_%j.err

# FastSurfer longitudinal per-patient launcher pinned to 8 GB GPUs (view aggregation on CPU).
set -e
uenv verbose singularity-4.2.2

SITE=$SITE
PATIENT=$PATIENT

echo "Running longitudinal FastSurfer"
echo "Site: $SITE"
echo "Patient: $PATIENT"

PW=/nfs/br1_prosjekt/ParkWest
DATA=$PW/user/2026vae/AmundEspen/pd_thesis/data/$SITE/$PATIENT
OUT=$PW/user/2026vae/AmundEspen/pd_thesis/outputs/fastsurfer_longitudinal/$SITE/$PATIENT
CTR=$PW/user/2026vae/AmundEspen/pd_thesis/containers/fastsurfer_sm120_sandbox
LICENSE=$PW/user/2026vae/AmundEspen/pd_thesis/fastsurfer/license/license.txt

mkdir -p "$OUT"

BL=$(find $DATA/BL -name "*.nii*" | head -n 1)
Y3=$(find $DATA/3Y -name "*.nii*" | head -n 1)
Y5=$(find $DATA/5Y -name "*.nii*" | head -n 1)

echo "BL: $BL"
echo "3Y: $Y3"
echo "5Y: $Y5"

singularity exec --nv --env FS_LICENSE=$LICENSE \
  --bind $PW:$PW \
  $CTR \
  /fastsurfer/long_fastsurfer.sh \
  --tid ${PATIENT}_template \
  --sd $OUT \
  --tpids BL 3Y 5Y \
  --t1s $BL $Y3 $Y5 \
  --parallel 4 \
  --parallel_surf 4 \
  --parallel_seg 1 \
  --viewagg_device cpu
