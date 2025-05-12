#!/bin/bash
#SBATCH -p mit_normal
#SBATCH -a 1-5                     # 5 array tasks = 5 chunks
#SBATCH -c 32                      # reserve 32 cores on ONE node
#SBATCH --mem=128G
#SBATCH -J hmmnod_2025
#SBATCH -e %x_%A_%a.err
#SBATCH -t 04:00:00

module purge
module load miniforge/23.11.0-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prodigal_env        # contains hmmer 3.x

# -------------------------------------------------------------------------
# Pick the N-th chunk file for this array task.
# Use 'printf %q' if filenames may contain spaces.
# -------------------------------------------------------------------------

CHUNKDIR=/orcd/pool/004/imgm/chunks

chunk=$(ls "${CHUNKDIR}"/*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z $chunk ]]; then
    echo "No chunk found for task $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

out=${chunk%.fasta}3nod.tbl

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID scanning $chunk -> $out"

hmmscan --cpu 32 -E 1e-20 \
        --domtblout "$out" \
        /orcd/pool/004/imgm/try3nod.hmm  "$chunk"

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID finished"
