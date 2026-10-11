#!/bin/bash
# =============================================================================
# PANGENOME INPUT PREP 01 — break Hi-C-unsupported RagTag joins
#
# Reads the ptarmigan-scaffolded assemblies (grouse_asm/final_ptarmigan) and
# writes one pangenome-input FASTA per haplotype in which every RagTag join
# that the assembly repo's step 13 (13_join_support_ptarmigan.sh) called 'unsupported' is cut at its 100 bp gap. Everything else
# is copied through unchanged. Details and safety checks are in
# 01_break_unsupported_joins.py; in short, it proves the AGP coordinates still
# match each FASTA (object length, and the cut gap must be all N) before it
# cuts anything, and stops on any mismatch.
#
# Run AFTER the assembly repo's relabel_male_W.sh (male chr_W -> scaffold_N), which this reads
# from W_relabel_log.tsv so male W pieces are recognised under their new name.
#
# USAGE (from the grouse_pangenome repo root)
#   DRY_RUN=true bash 01_break_unsupported_joins.sh     # login node: plan only, ~1 min
#   sbatch 01_break_unsupported_joins.sh                 # write the FASTAs
#
#   BREAK=unsupported,ambiguous sbatch --export=ALL 01_break_unsupported_joins.sh
#   BREAK=all                   sbatch --export=ALL 01_break_unsupported_joins.sh
#       'all' returns every chromosome to its yahs scaffolds (no joins at all).
#
# OUTPUT  ${PROJECT_DIR}/pangenome_inputs/   (PROJECT_DIR = $CLUSTER_SCRATCH/GROUSE/grouse_asm)
#   <SP>_<S>_<H>.pg_input.fasta (+ .fai)   one per haplotype
#   breaks.tsv  pieces.tsv  summary.tsv   what was cut, and the bp accounting
# =============================================================================
#SBATCH --job-name=pg01_break_joins
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# 8G: the script holds one chromosome at a time (largest ~195 Mb) in memory.
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

show_tsv() { if command -v column >/dev/null 2>&1; then column -t -s$'\t' "$1"; else cat "$1"; fi; }

# ---- paths ------------------------------------------------------------------
ASM_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"           # assembly/QC repo outputs
PG_WORK="${PG_WORK:-${ASM_DIR}}"                         # shared PROJECT_DIR, as steps 03-04

FINAL_DIR="${ASM_DIR}/final_ptarmigan"
RAGTAG_DIR="${ASM_DIR}/ragtag_ptarmigan"
JOINS_DIR="${ASM_DIR}/qc_ptarmigan/join_support"
W_LOG="${FINAL_DIR}/W_relabel_log.tsv"
OUT_DIR="${PG_WORK}/pangenome_inputs"

BREAK="${BREAK:-unsupported}"
# Joins RagTag made on the reference's UNPLACED scaffolds (objects named
# scaffold_N) were never scored by step 13. true = cut them all (default).
CUT_UNPLACED="${CUT_UNPLACED:-true}"
DRY_RUN="${DRY_RUN:-false}"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY_SCRIPT="${SCRIPT_DIR}/01_break_unsupported_joins.py"

# ---- environment --------------------------------------------------------------
set +u
module unload anaconda >/dev/null 2>&1 || true
module load biocontainers >/dev/null 2>&1 || true
module load samtools/1.22.1 >/dev/null 2>&1 || true
set -u
unset LD_PRELOAD    # RCAC XALT tracker; can kill subshells on some nodes

mkdir -p logs
for X in "$PY_SCRIPT"; do [[ -f "$X" ]] || { echo "ERROR: missing $X"; exit 1; }; done
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH"; exit 1; }
python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)' \
    || { echo "ERROR: python >= 3.6 required"; exit 1; }
for D in "$FINAL_DIR" "$RAGTAG_DIR"; do [[ -d "$D" ]] || { echo "ERROR: missing $D"; exit 1; }; done
if [[ "$BREAK" != "all" && ! -d "$JOINS_DIR" ]]; then
    echo "ERROR: missing $JOINS_DIR (step 13 output). Use BREAK=all to cut without it."
    exit 1
fi
if [[ ! -s "$W_LOG" ]]; then
    echo "WARNING: ${W_LOG} not found — run relabel_male_W.sh first if males still carry chr_W."
fi

ARGS=(--final-dir "$FINAL_DIR" --ragtag-dir "$RAGTAG_DIR" --out-dir "$OUT_DIR"
      --break "$BREAK" --w-log "$W_LOG")
[[ "$BREAK" != "all" ]] && ARGS+=(--joins-dir "$JOINS_DIR")
[[ -n "${ONLY:-}" ]]     && ARGS+=(--only "$ONLY")
[[ "$CUT_UNPLACED" == true ]] || ARGS+=(--keep-unplaced-joins)

echo ">>> break mode: ${BREAK}"
echo ">>> cut joins on unplaced placements: ${CUT_UNPLACED}"
echo ">>> input     : ${FINAL_DIR}"
echo ">>> output    : ${OUT_DIR}"

if [[ "$DRY_RUN" == true ]]; then
    python3 "$PY_SCRIPT" "${ARGS[@]}" --dry-run
    exit 0
fi

python3 "$PY_SCRIPT" "${ARGS[@]}"

# ---- index and reconcile ------------------------------------------------------
echo ""
echo ">>> indexing outputs"
if command -v samtools >/dev/null 2>&1; then
    for F in "${OUT_DIR}"/*.pg_input.fasta; do samtools faidx "$F"; done
else
    echo "  samtools not available — index later with: samtools faidx <file>"
fi

# bp accounting: output must equal input minus the N runs of the cut gaps
awk -F'\t' 'NR>1 { if ($7 != $3 - $8) { bad++; print "  MISMATCH: " $1 } }
            END   { if (bad) exit 1 }' "${OUT_DIR}/summary.tsv" \
    || { echo "ERROR: bp accounting failed — see above"; exit 1; }

echo ""
echo ">>> summary (${OUT_DIR}/summary.tsv)"
show_tsv "${OUT_DIR}/summary.tsv" | sed 's/^/  /'
echo ""
echo ">>> cuts by chromosome across all haplotypes"
awk -F'\t' 'NR>1{c[$2]++} END{for (k in c) printf "  %-12s %d\n", k, c[k]}' \
    "${OUT_DIR}/breaks.tsv" | sort -k2,2nr
echo ""
echo ">>> next: sbatch --array=0-22%6 02_assembly_qc_array.sh   (one task per manifest row)"
