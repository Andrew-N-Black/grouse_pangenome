#!/bin/bash
# Genome-wide chromosome QC, run BEFORE graph building.
# Usage:
#   ls /path/to/haps/*.fa > haps.txt
#   sbatch -A dewoody -p cpu -q standby -t 04:00:00 -n 16 --mem=32G \
#          --array=1-$(wc -l < haps.txt) run_chrom_qc.sh map anchor.fa haps.txt
#   bash run_chrom_qc.sh summarize anchor.fa haps.txt
#SBATCH -A dewoody
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 04:00:00
#SBATCH --mem=32G
#SBATCH -J chromqc
#SBATCH -o chromqc_%A_%a.log

set -euo pipefail
MODE=$1; ANCHOR=$2; LIST=$3
OUT=chromqc_paf; mkdir -p "$OUT"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if [[ $MODE == map ]]; then
    module load biocontainers minimap2 samtools
    HAP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
    base=$(basename "$HAP" | sed 's/\.f[a-z]*\(\.gz\)\?$//')
    [[ -s ${HAP}.fai ]] || samtools faidx "$HAP"
    minimap2 -x asm10 --secondary=no -t "${SLURM_NTASKS:-16}" "$ANCHOR" "$HAP" > "$OUT/${base}.paf"
elif [[ $MODE == summarize ]]; then
    python3 "$SCRIPT_DIR/chrom_qc.py" --paf "$OUT"/*.paf --fai $(sed 's/$/.fai/' "$LIST") --out chrom_qc.tsv
fi
