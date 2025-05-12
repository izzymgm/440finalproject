#!/bin/bash
#SBATCH -p mit_normal             # partition/queue
#SBATCH -N 1                      # 1 node
#SBATCH -c 32                     # 32 logical cores
#SBATCH --mem=128G                # (4 GB × 32); adjust as you like
#SBATCH -J hmmscan_caricao
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -t 08:00:00               # wall time

module purge
module load miniforge/23.11.0-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prodigal_env       # env that has hmmer

hmmscan --cpu 32 -E 1e-20 --domtblout caricaohits.domtbl \
        try3nod.hmm  combined_caricao_orfs.faa
