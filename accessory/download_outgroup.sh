#!/bin/bash
#SBATCH --job-name=sra_outgroup
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=12:00:00
#SBATCH --output=sra_outgroup_%j.out
#SBATCH -A dewoody
#SBATCH -p cpu

# Download Gunnison sage-grouse (Centrocercus minimus) outgroup reads from
# PRJNA531321 with sra-tools and place them alongside the grouse shotgun
# fastqs, named to match the *_R1_001.fastq.gz / *_R2_001.fastq.gz convention.
# If compute nodes have no internet access, run this on a login node with bash.

set -euo pipefail

module load biocontainers sra-tools/2.11.0-pl5262
module load pigz 2>/dev/null || true
command -v pigz >/dev/null || pigz() { gzip "${@/#-p*/}"; }   # fall back to gzip

OUTDIR=/scratch/gautschi/blackan/GROUSE/grouse_asm/raw/Shotgun
WORK=/scratch/gautschi/blackan/GROUSE/grouse_asm/raw/sra_tmp
THREADS=${SLURM_NTASKS:-8}

mkdir -p "$WORK"
cd "$WORK"

# Run accession -> sample name (from SraRunTable.csv)
declare -A RUNS=(
    [SRR8863303]=GB8
    [SRR8863311]=CR2
    [SRR8863315]=GB56
)

for srr in "${!RUNS[@]}"; do
    s=${RUNS[$srr]}
    r1="$OUTDIR/${s}_${srr}_R1_001.fastq.gz"
    r2="$OUTDIR/${s}_${srr}_R2_001.fastq.gz"

    if [[ -s "$r1" && -s "$r2" ]]; then
        echo "[$s] already done, skipping"
        continue
    fi

    echo "[$s] prefetch $srr"
    prefetch --max-size 50g -O "$WORK" "$srr"
    # 2.11 may save to $WORK/SRR.sra, $WORK/SRR/SRR.sra, or *.sralite
    sra=$(find "$WORK" /scratch/gautschi/blackan/sra_cache -maxdepth 3 -type f \
              -name "${srr}.sra*" 2>/dev/null | head -n 1)
    if [[ -z "$sra" ]]; then
        echo "ERROR: [$s] prefetch output for $srr not found under $WORK" >&2
        ls -lR "$WORK" >&2
        exit 1
    fi
    echo "[$s] using $sra"

    echo "[$s] fasterq-dump"
    fasterq-dump --split-files --threads "$THREADS" \
        --temp "$WORK" -O "$WORK" "$sra"

    # Paired reads should have identical line counts
    n1=$(wc -l < "${srr}_1.fastq")
    n2=$(wc -l < "${srr}_2.fastq")
    if [[ "$n1" -ne "$n2" ]]; then
        echo "ERROR: [$s] R1 ($n1 lines) and R2 ($n2 lines) differ" >&2
        exit 1
    fi
    echo "[$s] $((n1 / 4)) read pairs"

    pigz -p "$THREADS" -c "${srr}_1.fastq" > "${r1}.part" && mv "${r1}.part" "$r1"
    pigz -p "$THREADS" -c "${srr}_2.fastq" > "${r2}.part" && mv "${r2}.part" "$r2"

    rm -rf "${srr}_1.fastq" "${srr}_2.fastq" "$WORK/$srr" "$sra"
    echo "[$s] wrote $r1 and $r2"
done

ls -lh "$OUTDIR"/{GB8,CR2,GB56}_SRR*_R?_001.fastq.gz
