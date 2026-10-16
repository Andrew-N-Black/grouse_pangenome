#!/bin/bash
# =============================================================================
# STEP 05d: RESCUE THEx WINDOWS LOST TO ONE OR TWO MISSING HAPLOTYPES
#
# Run after 05_thex_phylogeny.sh has finished. THExBuilder drops a whole window
# when any single sequence in it is missing: trimal deletes all-gap sequences
# and THEx then discards windows whose sequence count differs from the rest of
# the chromosome. So one haplotype assembly that lacks a stretch (e.g. STGR
# F5503 hap2 on chr_Z 11.4-19.2 Mb) removes that stretch for all 47 taxa.
#
# This script retries every window that has no tree, without the haplotypes that
# cover < WIN_MIN_COV of it:
#   prep       prune those haplotypes (never the outgroup/reference; at most
#              MAX_DROP per window; every species keeps >= MIN_PER_SPECIES),
#              then thexb --trimal (sequence-count check off, since counts now
#              differ by window) and thexb --pw_filter with step 05's settings,
#              and merge the windows that pass into step 05's
#              pairwise_filtered_windows/
#   iqtree     step 05's own iqtree stage, on the rescued windows only
#   summarize  step 05's own summarize: rewrites work/thex_inputs/ with the
#              rescued windows included (their Taxa column shows the lower count)
#
# Windows that already have trees are not touched, so their trees stay valid.
# Log of what was dropped where: work/rescue/prune_log.tsv
#
# USAGE (Negishi, repo root):  bash 05d_thex_rescue.sh submit
#   then re-run 05b_thex_stats.py on the new thex_inputs/ files.
# =============================================================================
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail
unset LD_PRELOAD

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SCRIPT="${REPO_DIR}/05d_thex_rescue.sh"
STEP05="${REPO_DIR}/05_thex_phylogeny.sh"
PYHELP="${REPO_DIR}/05d_thex_rescue.py"

PROJECT_DIR="${PROJECT_DIR:-${CLUSTER_SCRATCH}/GROUSE/grouse_thex}"
WORK="${PROJECT_DIR}/work"
THEXB_DIR="${WORK}/thexb"
RESCUE="${WORK}/rescue"
IQ_DIR="${WORK}/iqtree"
ENV_PREFIX="${ENV_PREFIX:-${PROJECT_DIR}/envs/thex_env}"
OUTGROUP_NAME="${OUTGROUP_NAME:-Lmuta}"
SPECIES="${SPECIES:-STGR,LEPC,GRPC}"

ACCOUNT="${ACCOUNT:-fnrdewoody}"
PARTITION="${PARTITION:-cpu}"
RESCUE_CPUS="${RESCUE_CPUS:-64}"
IQ_ARRAY="${IQ_ARRAY:-4}"
IQ_CPUS="${IQ_CPUS:-64}"

# ---- pruning rules ----
WIN_MIN_COV="${WIN_MIN_COV:-0.5}"        # drop a haplotype from a window if < 50% of it is aligned
MAX_DROP="${MAX_DROP:-4}"                # retry only if at most this many haplotypes go
MIN_PER_SPECIES="${MIN_PER_SPECIES:-4}"  # and every species keeps at least this many

# ---- same filter settings as step 05 ----
TRIMAL_GT="${TRIMAL_GT:-0.9}"
TRIMAL_MINLEN="${TRIMAL_MINLEN:-1kb}"
PW_SUBWIN="${PW_SUBWIN:-100bp}"
PW_STEP="${PW_STEP:-10bp}"
PW_MAX_PDIST="${PW_MAX_PDIST:-0.15}"
PW_ZSCORE="${PW_ZSCORE:-2}"
PW_COV="${PW_COV:-0.9}"

MODE="${1:-}"
die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

stage_prep() {
    [[ -x "${ENV_PREFIX}/bin/thexb" ]] || die "THEx env not found at ${ENV_PREFIX}"
    export PATH="${ENV_PREFIX}/bin:${PATH}"
    [[ -s "${WORK}/ref_id.txt" ]] || die "step 05 has not run (no ${WORK}/ref_id.txt)"
    local REF_ID; REF_ID="$(cat "${WORK}/ref_id.txt")"
    local CPUS="${SLURM_CPUS_PER_TASK:-$RESCUE_CPUS}"

    rm -rf "$RESCUE"; mkdir -p "$RESCUE"
    python "$PYHELP" prune --windowed "${THEXB_DIR}/windowed_fastas" \
        --pw "${THEXB_DIR}/pairwise_filtered_windows" --iqtree "$IQ_DIR" \
        --out "${RESCUE}/windowed_fastas" --log "${RESCUE}/prune_log.tsv" \
        --keep "${REF_ID},${OUTGROUP_NAME}" --species "$SPECIES" \
        --min-cov "$WIN_MIN_COV" --max-drop "$MAX_DROP" --min-per-species "$MIN_PER_SPECIES"
    awk -F'\t' 'NR>1 {n[$3]++} END {for (k in n) printf "  %-20s %d windows\n", k, n[k]}' "${RESCUE}/prune_log.tsv"
    [[ -d "${RESCUE}/windowed_fastas" ]] || { log "nothing to retry"; : > "${WORK}/iqtree_windows.txt"; return 0; }

    cd "$RESCUE"
    log "thexb --trimal on the pruned windows (sequence-count check off)"
    # --trimal_drop_windows is a store_false flag: passing it turns the count check OFF
    thexb --trimal -i "${RESCUE}/windowed_fastas" -o "$RESCUE" --trimal_gap_threshold "$TRIMAL_GT" \
        --trimal_min_seq_len "$TRIMAL_MINLEN" --trimal_drop_windows --cpu "$CPUS"

    local maxp="$PW_MAX_PDIST" subwin_bp
    subwin_bp=$(python -c "
s='${PW_SUBWIN}'.lower(); m={'bp':1,'kb':1000,'mb':1000000}
u=[k for k in ('kb','mb','bp') if s.endswith(k)][0]; print(int(float(s[:-2])*m[u]))")
    if python -c "import inspect, sys, thexb.STAGE_pairwise_filter as m; sys.exit(0 if 'countmismatch + countmissing / denominator' in inspect.getsource(m) else 1)"; then
        maxp=$(python -c "print(${PW_MAX_PDIST} * ${subwin_bp})")
    fi
    log "thexb --pw_filter (as step 05: max p-dist $maxp, z $PW_ZSCORE, coverage $PW_COV)"
    thexb --pw_filter -i "$RESCUE" -o "$RESCUE" -r "$REF_ID" --pw_exclude "${REF_ID},${OUTGROUP_NAME}" \
        --pw_subwindow_size "$PW_SUBWIN" --pw_step "$PW_STEP" --pw_max_pdist "$maxp" \
        --pw_zscore_cutoff "$PW_ZSCORE" --pw_seq_coverage "$PW_COV" --pw_missing_char N --cpu "$CPUS"

    python "$PYHELP" merge --rescued-pw "${RESCUE}/pairwise_filtered_windows" \
        --pw "${THEXB_DIR}/pairwise_filtered_windows"
    # step 05's iqtree stage reads this list; point it at the rescued windows only
    [[ -s "${WORK}/iqtree_windows.step05.txt" ]] || cp "${WORK}/iqtree_windows.txt" "${WORK}/iqtree_windows.step05.txt"
    ( cd "$THEXB_DIR" && find pairwise_filtered_windows -name '*.fasta' ! -name '*-DROPPED*' | sort ) \
        | while read -r f; do
            c="$(basename "$(dirname "$f")")"; s="$(basename "$f" .fasta)"
            [[ -s "${IQ_DIR}/${c}/${s}.contree" ]] || echo "$f"
          done > "${WORK}/iqtree_windows.txt"
    log "$(wc -l < "${WORK}/iqtree_windows.txt") rescued windows go to IQ-TREE"
}

stage_submit() {
    mkdir -p "${REPO_DIR}/logs"
    [[ -s "$STEP05" && -s "$PYHELP" ]] || die "run from the repo root (needs 05_thex_phylogeny.sh and 05d_thex_rescue.py)"
    if squeue -u "${USER:-$(id -un)}" -h -o '%j' 2>/dev/null | grep -q '^thex_'; then
        die "thex_* jobs are still queued; wait for them or scancel them"
    fi
    local J1 J2 J3
    J1=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q standby -t 04:00:00 --export=ALL \
        -N 1 -c "$RESCUE_CPUS" -J thex_rescue -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" prep) \
        || die "rescue submission refused"
    log "rescue prep  $J1"
    J2=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q standby -t 04:00:00 --export=ALL \
        --dependency=afterok:"$J1" --array=0-$((IQ_ARRAY-1)) -N 1 -c "$IQ_CPUS" \
        -J thex_iqtree -o "${REPO_DIR}/logs/%x_%A_%a.out" "$STEP05" iqtree) || die "iqtree submission refused"
    log "iqtree       $J2 ($IQ_ARRAY tasks)"
    J3=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q standby -t 04:00:00 --export=ALL \
        --dependency=afterany:"$J2" -c 8 --mem=16G \
        -J thex_summary -o "${REPO_DIR}/logs/%x_%j.out" "$STEP05" summarize) || die "summary submission refused"
    log "summarize    $J3  -> rewrites ${WORK}/thex_inputs/"
}

case "$MODE" in
    submit) stage_submit ;;
    prep)   stage_prep ;;
    *) sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
