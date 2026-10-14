#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: 3-SPECIES PANGENOME + COMPARATIVE ANALYSIS
# Step 03 — requires steps 01-02 of this repo, which in turn read the ASSEMBLY
# pipeline's output (separate repo: prairie_grouse). The GRAPH is built from
# step 01's pangenome_inputs/ FASTAs (Hi-C-unsupported RagTag joins cut);
# the TE and gene comparisons read the assembly pipeline's final_ptarmigan/
# directly (hifiasm + yahs + RagTag pseudo-chromosomes, with Liftoff gene
# annotations transferred from the rock ptarmigan reference), since cutting
# a join removes only a 100 bp gap and the Liftoff GFF3 seqids name the
# uncut chromosomes.
#
# Builds ONE combined pangenome graph across all haplotype assemblies from
# all three species (PGGB), then uses it — plus the existing Liftoff/TE
# annotations — to compare Sharp-tailed Grouse, Lesser Prairie-Chicken, and
# Greater Prairie-Chicken on:
#   - assembly length (per-haplotype path length in the graph)
#   - transposable element content (RepeatModeler + RepeatMasker)
#   - gene content (Liftoff-mapped ptarmigan gene counts)
#
# This is NOT an array job: a pangenome is inherently a joint analysis across
# all samples, so everything runs in one job. The one genuinely
# embarrassingly-parallel sub-step (RepeatMasker per genome) is parallelized
# internally with `xargs -P`.
#
# HONEST SCOPE NOTE: whole-genome PGGB across ~46 haplotypes of ~1.1 Gb bird
# genomes is a heavy, potentially multi-day-to-multi-week alignment job. The
# resource/time requests below are rough starting points, not a guarantee —
# check them against what your Purdue RCAC allocation actually offers, and
# consider per-chromosome partitioning (pggb's `partition-before-pggb`
# helper) if the whole-genome run doesn't finish in your walltime. Because
# your assemblies are already RagTag-oriented against the same rock ptarmigan
# reference, chromosome-level partitioning should work cleanly if needed.
#
# RUN MODES (all from the repo root)
#   Whole genome, one job (original behaviour):
#       sbatch 03_pangenome_build.sh
#   Per chromosome (recommended; graphs are built independently per ptarmigan
#   chromosome, which is valid because every haplotype is RagTag-oriented to
#   the same reference, and pieces cut in step 01 keep their chromosome name):
#       bash 03_pangenome_build.sh list            # writes pangenome/chromosomes.txt, prints the array command
#       CHROM=chr_20 sbatch --export=ALL -t 2-00:00:00 --mem=128G 03_pangenome_build.sh   # one-chromosome test
#       sbatch --export=ALL,CHROM=array --array=0-<N-1>%8 -t 4-00:00:00 --mem=250G 03_pangenome_build.sh
#       STAGE=summary sbatch --export=ALL -t 4-00:00:00 03_pangenome_build.sh
#   The summary stage merges the per-chromosome graphs into one
#   (pangenome/merged/pangenome.final.{og,gfa}, used by step 04) and runs the
#   genome-wide steps (graph stats, panacus, TE annotation, species tables).
#   Per-chromosome jobs run steps 1-5 only, into pangenome/<CHROM>/ and
#   comparison/<CHROM>/. CHROM=unplaced collects every sequence not named
#   chr_* (unplaced scaffolds, relabelled male W pieces).
#
# SOFTWARE: Gautschi has no pggb, odgi, vg or RepeatMasker module. Each tool is
#   used from a module when one is on PATH, otherwise from an Apptainer image
#   in ${PROJECT_DIR}/containers (fetch once with accessory/pull_containers.sh).
#
# A NOTE ON CHROMOSOME NAMES
#   The assemblies carry the ptarmigan reference's own chromosome numbering,
#   which is a DESCENDING-LENGTH RANK, not homology to the chicken karyotype:
#   chr_6 here is not chicken chromosome 6. Nothing in this script depends on
#   that, because every assembly uses the same convention — but any result
#   reported per chromosome has to be translated through the correspondence
#   table the assembly repo's step 04 produces before it can be compared
#   against the galliform literature.
#
# USAGE:
#   sbatch 03_pangenome_build.sh
# =============================================================================
#SBATCH --job-name=grouse_pangenome
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
# Literal address: SBATCH directives are comments to the shell, so ${USER}
# would NOT be expanded here, and the mail would go to a user named "${USER}".
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

# unset LD_PRELOAD: RCAC's XALT library is injected this way and fails on some
# nodes (GLIBC mismatch), which can kill subshells under `set -e`.
unset LD_PRELOAD

ml biocontainers
ml samtools/1.22.1
ml repeatmodeler/2.0.3 2>/dev/null || true   # bioconda build also ships RepeatMasker/BuildDatabase
# (no python module: the system python3 is used)
# No panacus module and no cargo/Rust toolchain on this cluster, so it's
# fetched below as a precompiled static (musl) binary. See PANACUS_BIN.

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR:-$PWD}/assembly_manifest.tsv"

# Same PROJECT_DIR as the assembly pipeline — this reads its output in place.
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"

# Which reference the assemblies were ordered against. The assembly pipeline
# suffixes its output trees with this, so a run against a different reference
# produces a parallel tree and both can be kept side by side. Must match
# REF_TAG in that pipeline's 05_ragtag_liftoff_array.sh.
REF_TAG="ptarmigan"

FINAL_DIR="${PROJECT_DIR}/final_${REF_TAG}"      # input: assemblies + Liftoff GFF3 (TE/gene tables)
INPUT_DIR="${PROJECT_DIR}/pangenome_inputs"     # input: graph FASTAs from step 01
LIFTOFF_DIR="${PROJECT_DIR}/liftoff_${REF_TAG}"  # input: unmapped-feature lists
PANGENOME_DIR="${PROJECT_DIR}/pangenome"
TE_DIR="${PROJECT_DIR}/te_annotation"
COMPARE_DIR="${PROJECT_DIR}/comparison"

# ---- run mode -----------------------------------------------------------------
#   CHROM=""        whole genome in one graph (default)
#   CHROM=chr_N     one chromosome; CHROM=unplaced: everything not named chr_*
#   CHROM=array     take the chromosome from pangenome/chromosomes.txt by
#                   SLURM_ARRAY_TASK_ID (write it first with: bash $0 list)
#   STAGE=summary   merge per-chromosome graphs, then genome-wide steps 4-7
#   STAGE=te        TE annotation only (step 6); can run while the graphs build
CHROM="${CHROM:-}"
STAGE="${STAGE:-}"
MODE="${1:-}"
CHROM_LIST="${PANGENOME_DIR}/chromosomes.txt"

# ---- tools: module if available, else Apptainer image ---------------------------
CONTAINER_DIR="${PROJECT_DIR}/containers"
PGGB_SIF="${PGGB_SIF:-${CONTAINER_DIR}/pggb.sif}"        # pggb, odgi, bgzip, samtools
TETOOLS_SIF="${TETOOLS_SIF:-${CONTAINER_DIR}/tetools.sif}"   # RepeatModeler, RepeatMasker
APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
export APPTAINERENV_LD_PRELOAD="" SINGULARITYENV_LD_PRELOAD=""
# --cleanenv: the host shell exports wrapper functions (pggb, odgi, bgzip, cx); without it they
# leak into the container, so pggb's own internal odgi/bgzip calls would try to re-launch apptainer.
cx() { local sif="$1"; shift; "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$sif" "$@"; }
use_tool() {   # use_tool <command> <image>: define a wrapper if <command> is not on PATH
    local t="$1" sif="$2"
    command -v "$t" >/dev/null 2>&1 && return 0
    if [[ -z "$APPTAINER_BIN" || ! -s "$sif" ]]; then
        echo "ERROR: '$t' is not on PATH and image ${sif} is missing."
        echo "       Run accessory/pull_containers.sh on a login node first."; exit 1
    fi
    eval "${t}() { cx \"${sif}\" ${t} \"\$@\"; }"; export -f "$t"
}
export -f cx; export APPTAINER_BIN CLUSTER_SCRATCH

if [[ "$CHROM" == "array" ]]; then
    [[ -s "$CHROM_LIST" && -n "${SLURM_ARRAY_TASK_ID:-}" ]] \
        || { echo "ERROR: CHROM=array needs ${CHROM_LIST} (bash $0 list) and an --array job"; exit 1; }
    CHROM=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CHROM_LIST")
    [[ -n "$CHROM" ]] || { echo "ERROR: no chromosome at line $((SLURM_ARRAY_TASK_ID + 1)) of ${CHROM_LIST}"; exit 1; }
fi
if [[ -n "$CHROM" ]]; then
    PANGENOME_DIR="${PANGENOME_DIR}/${CHROM}"
    COMPARE_DIR="${COMPARE_DIR}/${CHROM}"
fi

COMBINED_FASTA="${PANGENOME_DIR}/combined_haplotypes.fa"

# panacus — precompiled static (musl) binary, cached under PROJECT_DIR so
# this only downloads once across reruns. No cargo/Rust toolchain required.
PANACUS_VERSION="0.5.1"
PANACUS_BIN="${PROJECT_DIR}/panacus-${PANACUS_VERSION}/bin/panacus"

# PGGB alignment parameters — moderate cross-species divergence within the
# genus Tympanuchus. See https://pggb.readthedocs.io/en/latest/rst/essential_parameters.html
PGGB_MIN_IDENTITY=90
PGGB_SEGMENT_LEN="5k"
# Extra pggb options for one job, e.g. for a very large chromosome:
#   PGGB_EXTRA="-B 2000000"   smaller seqwish transitive-closure batches (less
#                             memory, slower; default 10,000,000 bp)
#   PGGB_EXTRA="--resume"     reuse finished stage outputs. pggb writes each
#                             stage in place, so first delete the file of the
#                             stage that was running when the job died.
PGGB_EXTRA="${PGGB_EXTRA:-}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
# Concurrency for the per-genome RepeatMasker loop (each call still gets
# multiple threads via -pa); tune alongside THREADS.
TE_PARALLEL_JOBS=4
TE_THREADS_PER_JOB=$(( THREADS / TE_PARALLEL_JOBS > 0 ? THREADS / TE_PARALLEL_JOBS : 1 ))

mkdir -p logs "$PANGENOME_DIR" "$TE_DIR" "$COMPARE_DIR"

# ---- `bash 03_pangenome_build.sh list`: chromosome list for the array ----------
if [[ "$MODE" == "list" ]]; then
    ls "${INPUT_DIR}"/*.pg_input.fasta.fai >/dev/null 2>&1 \
        || { echo "ERROR: no ${INPUT_DIR}/*.pg_input.fasta.fai (run step 01)"; exit 1; }
    cat "${INPUT_DIR}"/*.pg_input.fasta.fai | cut -f1 | sed -E 's/_p[0-9]+$//' \
        | grep '^chr_' | sort -u -V > "$CHROM_LIST"
    echo "unplaced" >> "$CHROM_LIST"
    N=$(wc -l < "$CHROM_LIST")
    echo "Wrote ${CHROM_LIST} (${N} partitions):"; tr '\n' ' ' < "$CHROM_LIST"; echo
    echo "Submit: sbatch --export=ALL,CHROM=array --array=0-$((N-1))%8 -t 4-00:00:00 --mem=250G 03_pangenome_build.sh"
    exit 0
fi

for t in pggb odgi bgzip; do use_tool "$t" "$PGGB_SIF"; done
if [[ -z "$CHROM" ]]; then
    for t in BuildDatabase RepeatModeler RepeatMasker; do use_tool "$t" "$TETOOLS_SIF"; done
fi
echo ">>> Mode: ${CHROM:+chromosome }${CHROM:-${STAGE:-whole genome}}"

# =============================================================================
# INSTALL panacus (once — skipped on reruns if already present)
# No module and no cargo/Rust toolchain on this cluster, so we pull the
# precompiled musl-linked binary from GitHub releases instead of building
# from source. musl = statically linked, so it runs on any x86_64 Linux
# node regardless of glibc version — no runtime dependencies to manage.
# =============================================================================
if [[ ! -x "$PANACUS_BIN" ]]; then
    echo ">>> Installing panacus v${PANACUS_VERSION} (precompiled binary)"
    PANACUS_TARBALL="panacus-${PANACUS_VERSION}_x86_64-unknown-linux-musl.tar.gz"
    PANACUS_INSTALL_DIR="${PROJECT_DIR}/panacus-${PANACUS_VERSION}"
    PANACUS_TMP="${PROJECT_DIR}/.panacus_extract_tmp"

    rm -rf "$PANACUS_TMP"
    mkdir -p "$PANACUS_TMP" "$PANACUS_INSTALL_DIR"
    wget -q -O "${PANACUS_TMP}/${PANACUS_TARBALL}" \
        "https://github.com/marschall-lab/panacus/releases/download/v${PANACUS_VERSION}/${PANACUS_TARBALL}"
    tar -xzf "${PANACUS_TMP}/${PANACUS_TARBALL}" -C "$PANACUS_TMP"

    # Move the archive's inner contents up into PANACUS_INSTALL_DIR rather
    # than relying on `tar --strip-components` — tested interactively and it
    # did not flatten the release's top-level directory as expected, so this
    # is deliberately independent of that archive's exact internal layout.
    PANACUS_INNER_DIR=$(find "$PANACUS_TMP" -mindepth 1 -maxdepth 1 -type d | head -n1)
    if [[ -z "$PANACUS_INNER_DIR" ]]; then
        echo "ERROR: panacus archive did not extract as expected — check ${PANACUS_TMP}"
        exit 1
    fi
    mv "$PANACUS_INNER_DIR"/* "$PANACUS_INSTALL_DIR"/
    rm -rf "$PANACUS_TMP"
fi

if [[ ! -x "$PANACUS_BIN" ]]; then
    echo "ERROR: panacus install failed — expected binary not found at ${PANACUS_BIN}"
    exit 1
fi

echo ">>> 03_pangenome_build.sh"
echo ">>> Node      : $(hostname)"
echo ">>> CPUs      : ${THREADS}"
echo ">>> Start time: $(date)"

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

# =============================================================================
# MANIFEST SANITY CHECK
# Every assembly path below is built as <SPECIES>_<SAMPLE>_<HAP>.pg_input.fasta
# (graph) or .pseudo_chr.fasta (TE/gene tables)
# from the manifest's own fields. A stray leading or trailing space in either
# field therefore changes the path it looks for, and the only symptom is a
# "WARNING: missing assembly ... skipping" line buried in the log — the
# pangenome then quietly gets built from a subset of the panel. Catch it here
# instead. This does NOT trim the fields: if the assembly pipeline ran with the
# same whitespace, the files on disk carry it too, and silently trimming would
# break paths that currently work. Fix the manifest and re-run both pipelines,
# or leave it alone deliberately.
# =============================================================================
WS_ROWS=$(awk -F'\t' '
    !/^#/ && FNR > 1 && NF > 1 {
        if ($1 ~ /^[ \t]|[ \t]$/ || $2 ~ /^[ \t]|[ \t]$/)
            printf "    line %d: sample=[%s] species=[%s]\n", FNR, $1, $2
    }' "$MANIFEST")
if [[ -n "$WS_ROWS" ]]; then
    echo ""
    echo "  !! WARNING: leading/trailing whitespace in sample or species fields:"
    echo "$WS_ROWS"
    echo "     Assembly filenames are built from these fields verbatim, so these"
    echo "     rows will only be found if the assembly pipeline ran with the same"
    echo "     whitespace. Check the 'Collected N haplotype assemblies' count"
    echo "     below against the number of rows in the manifest before trusting"
    echo "     the graph."
    echo ""
fi

# =============================================================================
# STEP 1: Collect phased assemblies, rename headers to PanSN-spec
# (sample#haplotype#contig), and build one combined multi-FASTA.
# Species and sample are read from the manifest (not parsed back out of
# filenames) since either field may itself contain underscores.
# =============================================================================
echo ">>> Step 1: Collecting phased assemblies -> ${COMBINED_FASTA}${CHROM:+ (chromosome ${CHROM} only)}"

sanitize() {
    # printf, not echo: echo's trailing newline would become a trailing '_'
    printf '%s' "$1" | tr -c 'A-Za-z0-9_-' '_'
}

N_HAPS=0
mapfile -t MANIFEST_ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2)

# STAGE=te skips the graph (steps 1-5) and runs only the TE annotation (step 6)
if [[ "$STAGE" != "te" ]]; then
if [[ "$STAGE" != "summary" ]]; then
> "$COMBINED_FASTA"

for LINE in "${MANIFEST_ROWS[@]}"; do
    IFS=$'\t' read -r SAMPLE SPECIES _ _ _ _ <<< "$LINE"
    SPECIES="${SPECIES//[[:space:]]/}"   # manifest species may carry stray spaces ("LEPC ")
    [[ -z "$SAMPLE" ]] && continue

    SAMPLE_TAG=$(sanitize "$SAMPLE")
    SPECIES_TAG=$(sanitize "$SPECIES")
    PANSN_SAMPLE="${SPECIES_TAG}.${SAMPLE_TAG}"

    for HAP in hap1 hap2; do
        ASSEMBLY="${INPUT_DIR}/${SPECIES}_${SAMPLE}_${HAP}.pg_input.fasta"

        if [[ ! -f "$ASSEMBLY" ]]; then
            echo "  WARNING: missing assembly for ${SAMPLE} (${SPECIES}) ${HAP} — skipping: ${ASSEMBLY}" >&2
            continue
        fi

        # CHROM set: keep only that chromosome (and its step-01 pieces, <chrom>_pN);
        # CHROM=unplaced: keep every sequence not named chr_*.
        # PanSN haplotype field must be an integer (hap1 -> 1) or pggb warns
        awk -v pre="${PANSN_SAMPLE}#${HAP#hap}" -v chrom="$CHROM" '
            /^>/ { split($0, a, " "); name = substr(a[1], 2); base = name; sub(/_p[0-9]+$/, "", base)
                   keep = (chrom == "") || (chrom == "unplaced" ? base !~ /^chr_/ : base == chrom)
                   if (keep) { print ">" pre "#" name; n++ }
                   next }
            keep { print }
            END { exit (n > 0 ? 0 : 3) }
        ' "$ASSEMBLY" >> "$COMBINED_FASTA" && N_HAPS=$((N_HAPS + 1)) || true
    done
done

echo "  Collected ${N_HAPS} haplotype assemblies"

if [[ -n "$CHROM" ]]; then
    echo "  (${N_HAPS} haplotypes carry ${CHROM}; e.g. chr_W is present only in females)"
fi

if [[ "$N_HAPS" -lt 2 ]]; then
    echo "ERROR: fewer than 2 haplotype assemblies found — nothing to build a pangenome from."
    echo "       Looked in ${INPUT_DIR} for <SPECIES>_<SAMPLE>_hap{1,2}.pg_input.fasta"
    echo "       Check that step 01 (01_break_unsupported_joins.sh) has"
    echo "       completed for at least one sample, and that REF_TAG matches."
    exit 1
fi

# =============================================================================
# STEP 2: bgzip + index the combined FASTA (pggb's expected input format)
# =============================================================================
echo ">>> Step 2: Compressing and indexing combined FASTA"

if [[ ! -f "${COMBINED_FASTA}.gz" ]]; then
    bgzip -@ "$THREADS" -f "$COMBINED_FASTA"
fi
samtools faidx "${COMBINED_FASTA}.gz"

# =============================================================================
# STEP 3: Build the pangenome graph with pggb
# Skipped if a prior run already produced a final GFA (heavy step — safe to
# resubmit without redoing it).
# =============================================================================
echo ">>> Step 3: pggb pangenome construction (n=${N_HAPS} haplotypes)"

EXISTING_GFA=$(find "$PANGENOME_DIR" -maxdepth 1 -name "*.final.gfa" 2>/dev/null | head -n1)

if [[ -n "$EXISTING_GFA" ]]; then
    echo "  Existing pangenome graph found — skipping pggb run: ${EXISTING_GFA}"
else
    pggb \
        -i "${COMBINED_FASTA}.gz" \
        -o "$PANGENOME_DIR" \
        -n "$N_HAPS" \
        -p "$PGGB_MIN_IDENTITY" \
        -s "$PGGB_SEGMENT_LEN" \
        -t "$THREADS" \
        ${PGGB_EXTRA:-}
fi

GFA=$(find "$PANGENOME_DIR" -maxdepth 1 -name "*.final.gfa" 2>/dev/null | head -n1)
OG=$(find "$PANGENOME_DIR" -maxdepth 1 -name "*.final.og" 2>/dev/null | head -n1)

if [[ -z "$GFA" || -z "$OG" ]]; then
    echo "ERROR: pggb did not produce the expected *.final.gfa / *.final.og in ${PANGENOME_DIR}"
    exit 1
fi

echo "  Pangenome graph: ${GFA}"

else
# =============================================================================
# STAGE=summary: merge the per-chromosome graphs into one graph
# =============================================================================
    echo ">>> Summary stage: merging per-chromosome graphs"
    [[ -s "$CHROM_LIST" ]] || { echo "ERROR: ${CHROM_LIST} not found"; exit 1; }
    MERGED_DIR="${PANGENOME_DIR}/merged"; mkdir -p "$MERGED_DIR"
    OG_LIST="${MERGED_DIR}/graphs.txt"; > "$OG_LIST"; MISSING=0
    while read -r C; do
        G=$(find "${PANGENOME_DIR}/${C}" -maxdepth 1 -name "*.final.og" 2>/dev/null | head -n1)
        if [[ -n "$G" ]]; then echo "$G" >> "$OG_LIST"; else echo "  MISSING graph for ${C}"; MISSING=$((MISSING+1)); fi
    done < "$CHROM_LIST"
    (( MISSING == 0 )) || { echo "ERROR: ${MISSING} chromosome graph(s) missing — rerun those array tasks first"; exit 1; }
    OG="${MERGED_DIR}/pangenome.final.og"; GFA="${MERGED_DIR}/pangenome.final.gfa"
    [[ -s "$OG" ]]  || odgi squeeze -f "$OG_LIST" -o "$OG" -t "$THREADS"
    [[ -s "$GFA" ]] || odgi view -i "$OG" -g > "$GFA"
    echo "  Merged graph: ${GFA} ($(wc -l < "$OG_LIST") chromosome graphs)"
fi

# =============================================================================
# STEP 4: Graph-level summary + per-haplotype path (assembly) lengths
# =============================================================================
echo ">>> Step 4: Pangenome graph summary and per-haplotype lengths"

odgi stats -i "$OG" -S > "${COMPARE_DIR}/pangenome_graph_summary.yaml"

# `odgi paths -H -D '#'` prints a header line, then columns:
#   group.name  path.name  path.length  path.step.count  [node.1 ... node.N]
# group.name is everything before the first '#' in path.name — for our
# PanSN names (species.sample#1|2#contig) that's "species.sample".
odgi paths -i "$OG" -H -D '#' -t "$THREADS" > "${COMPARE_DIR}/path_lengths_raw.tsv"
echo -e "species\tsample\thaplotype\tcontig\tpath_length_bp" > "${COMPARE_DIR}/genome_length_by_haplotype.tsv"
awk -F'\t' 'NR>1 {
    split($1, sp, ".")          # $1=group.name=species.sample -> sp[1]=species sp[2]=sample
    split($2, p, "#")           # $2=path.name=species.sample#hap#contig -> p[2]=hap p[3]=contig
    print sp[1]"\t"sp[2]"\thap"p[2]"\t"p[3]"\t"$3   # PanSN hap "1" -> "hap1" to match manifest-built keys
}' "${COMPARE_DIR}/path_lengths_raw.tsv" >> "${COMPARE_DIR}/genome_length_by_haplotype.tsv"

echo "  Wrote ${COMPARE_DIR}/genome_length_by_haplotype.tsv"

# =============================================================================
# STEP 5: Core/accessory sequence content by species (panacus)
# Wrapped so an unexpected panacus grouping-format issue doesn't take down
# the rest of the pipeline — TE and gene comparisons below don't depend on it.
# =============================================================================
echo ">>> Step 5: Core/accessory pangenome growth by species (panacus)"

odgi paths -i "$OG" -L -t "$THREADS" > "${COMPARE_DIR}/path_names.txt"
# Growth curves per haplotype (-H merges a haplotype's pieces, e.g. chr_2_p1
# and chr_2_p2), counted in bp, for quorum 0 (pangenome: in any haplotype) and
# quorum 1 (core: in every haplotype). One table for all 46 haplotypes and one
# per species. Per-species curves use a copy of the GFA that keeps only that
# species' P-lines: segments no kept path visits have coverage 0 and drop out.
# (panacus -s/--subset was tried and does not restrict the haplotype groups;
# the earlier two-column -s file made panacus panic, exit 101.)
cut -d'.' -f1 "${COMPARE_DIR}/path_names.txt" | sort -u > "${COMPARE_DIR}/panacus_species.txt"
rm -f "${COMPARE_DIR}/panacus_groups.tsv" "${COMPARE_DIR}/pangenome_histgrowth_by_species.tsv"

set +e
PANACUS_STATUS=0
"$PANACUS_BIN" histgrowth -t "$THREADS" -c bp -l 1,1 -q 0,1 -H "$GFA" \
    > "${COMPARE_DIR}/pangenome_histgrowth_all.tsv" || PANACUS_STATUS=$?
while read -r SP; do
    SP_GFA="${COMPARE_DIR}/.${SP}.paths_only.gfa"
    awk -F'\t' -v sp="${SP}." '$1 != "P" || index($2, sp) == 1' "$GFA" > "$SP_GFA"
    "$PANACUS_BIN" histgrowth -t "$THREADS" -c bp -l 1,1 -q 0,1 -H "$SP_GFA" \
        > "${COMPARE_DIR}/pangenome_histgrowth_${SP}.tsv" || PANACUS_STATUS=$?
    rm -f "$SP_GFA"
done < "${COMPARE_DIR}/panacus_species.txt"
set -e

if [[ "$PANACUS_STATUS" -ne 0 ]]; then
    echo "  WARNING: panacus histgrowth failed (exit ${PANACUS_STATUS})."
    echo "  See '${PANACUS_BIN} histgrowth --help'"
    echo "  and rerun manually — this does not block the TE/gene comparisons below."
fi

fi   # end of: STAGE != te

if [[ -n "$CHROM" ]]; then
    echo ""
    echo ">>> Chromosome ${CHROM} done: graph ${GFA}"
    echo "    TE/gene tables and the merged graph come from: STAGE=summary sbatch 03_pangenome_build.sh"
    echo ">>> End time: $(date)"
    exit 0
fi

# =============================================================================
# STEP 6: Transposable element annotation
# One de novo RepeatModeler2 library per species (built from one
# representative haplotype), merged, then RepeatMasker run against every
# haplotype assembly so both species-level and within-species TE content
# can be reported.
# =============================================================================
echo ">>> Step 6: Transposable element annotation"

MERGED_TE_LIB="${TE_DIR}/merged_TE_library.fa"
> "$MERGED_TE_LIB"

declare -A SEEN_SPECIES

for LINE in "${MANIFEST_ROWS[@]}"; do
    IFS=$'\t' read -r SAMPLE SPECIES _ _ _ _ <<< "$LINE"
    SPECIES="${SPECIES//[[:space:]]/}"   # manifest species may carry stray spaces ("LEPC ")
    [[ -z "$SAMPLE" ]] && continue
    [[ -n "${SEEN_SPECIES[$SPECIES]:-}" ]] && continue

    REP_ASSEMBLY="${FINAL_DIR}/${SPECIES}_${SAMPLE}_hap1.pseudo_chr.fasta"
    [[ ! -f "$REP_ASSEMBLY" ]] && continue

    SEEN_SPECIES[$SPECIES]=1
    SPECIES_TAG=$(sanitize "$SPECIES")
    DB_NAME="${TE_DIR}/${SPECIES_TAG}_repeatdb"
    FAMILIES_FA="${DB_NAME}-families.fa"

    echo "  Building RepeatModeler library for ${SPECIES} (representative: ${SAMPLE})"

    if [[ ! -f "$FAMILIES_FA" ]]; then
        # run inside TE_DIR: RepeatModeler writes its RM_<pid>.<date> work folder to the cwd
        ( cd "$TE_DIR" && BuildDatabase -name "$DB_NAME" "$REP_ASSEMBLY" \
          && RM_HELP=$(RepeatModeler -help 2>&1 || true) \
          && if [[ "$RM_HELP" == *"-threads"* ]]; then
                 RepeatModeler -database "$DB_NAME" -threads "$THREADS" -LTRStruct
             else
                 # RepeatModeler < 2.0.4 (e.g. the repeatmodeler/2.0.3 module) has no -threads;
                 # -pa N runs N parallel searches of ~4 threads each
                 RepeatModeler -database "$DB_NAME" -pa $(( THREADS / 4 > 0 ? THREADS / 4 : 1 )) -LTRStruct
             fi )
    else
        echo "    Library already built — skipping"
    fi

    cat "$FAMILIES_FA" >> "$MERGED_TE_LIB"
done

echo "  Merged TE library: ${MERGED_TE_LIB} ($(grep -c '^>' "$MERGED_TE_LIB") consensus sequences)"

echo ">>> Step 6b: RepeatMasker on every haplotype assembly (parallel x${TE_PARALLEL_JOBS})"

RM_JOBLIST="${TE_DIR}/repeatmasker_jobs.txt"
> "$RM_JOBLIST"

for LINE in "${MANIFEST_ROWS[@]}"; do
    IFS=$'\t' read -r SAMPLE SPECIES _ _ _ _ <<< "$LINE"
    SPECIES="${SPECIES//[[:space:]]/}"   # manifest species may carry stray spaces ("LEPC ")
    [[ -z "$SAMPLE" ]] && continue
    for HAP in hap1 hap2; do
        ASSEMBLY="${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.pseudo_chr.fasta"
        [[ -f "$ASSEMBLY" ]] && echo "$ASSEMBLY" >> "$RM_JOBLIST"
    done
done

run_repeatmasker() {
    local assembly="$1"
    local outdir="${TE_DIR}/$(basename "${assembly%.fasta}")"
    mkdir -p "$outdir"
    if [[ -f "${outdir}/$(basename "$assembly").tbl" ]]; then
        echo "    Already masked — skipping: ${assembly}"
        return 0
    fi
    RepeatMasker -pa "$TE_THREADS_PER_JOB" -lib "$MERGED_TE_LIB" -gff -dir "$outdir" "$assembly"
}
export -f run_repeatmasker sanitize
declare -F RepeatMasker >/dev/null && export -f RepeatMasker
export MERGED_TE_LIB TE_DIR TE_THREADS_PER_JOB

xargs -a "$RM_JOBLIST" -I{} -P "$TE_PARALLEL_JOBS" bash -c 'run_repeatmasker "$@"' _ {}

if [[ "$STAGE" == "te" ]]; then
    echo ">>> TE stage done: libraries and RepeatMasker output in ${TE_DIR}"
    echo "    STAGE=summary reuses them (finished species/haplotypes are skipped) and builds the tables."
    echo ">>> End time: $(date)"
    exit 0
fi

# =============================================================================
# STEP 7: Summarize TE content and gene content (Liftoff) per haplotype,
# then aggregate everything by species.
# =============================================================================
echo ">>> Step 7: Aggregating TE + gene content + length by species"

python3 - "$COMPARE_DIR" "$FINAL_DIR" "$LIFTOFF_DIR" "$TE_DIR" "$MANIFEST" << 'PYEOF'
import csv, os, re, statistics, sys

compare_dir, final_dir, liftoff_dir, te_dir, manifest_path = sys.argv[1:6]

# --- read manifest: ordered list of (species, sample) ---
with open(manifest_path) as fh:
    data_rows = [l for l in fh if not l.startswith('#')]
    reader = csv.DictReader(data_rows, delimiter='\t')
    samples = [(row['species'], row['sample_id']) for row in reader if row.get('sample_id')]

records = {}  # (species, sample, hap) -> dict of metrics

def rec(species, sample, hap):
    return records.setdefault((species, sample, hap), {
        'species': species, 'sample': sample, 'haplotype': hap,
        'length_bp': None, 'genes_mapped': None, 'genes_unmapped': None,
        'te_pct_masked': None,
    })

# --- genome length, from the pangenome path table written in Step 4 ---
length_tsv = os.path.join(compare_dir, 'genome_length_by_haplotype.tsv')
per_hap_len = {}
if os.path.exists(length_tsv):
    with open(length_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = (row['species'], row['sample'], row['haplotype'])
            per_hap_len[key] = per_hap_len.get(key, 0) + int(row['path_length_bp'])
for key, total in per_hap_len.items():
    rec(*key)['length_bp'] = total

# All paths below are constructed forward from (species, sample, hap) taken
# straight from the manifest — never reverse-parsed out of filenames, since
# species/sample values may themselves contain underscores.
te_re = re.compile(r'Total interspersed repeats.*?([\d.]+)\s*%')

for species, sample in samples:
    for hap in ('hap1', 'hap2'):
        # --- gene content, from Liftoff GFF3 + unmapped-feature lists ---
        gff = os.path.join(final_dir, f"{species}_{sample}_{hap}.liftoff.gff3")
        if os.path.exists(gff):
            n_genes = 0
            with open(gff) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    if len(fields) > 2 and fields[2] == 'gene':
                        n_genes += 1
            rec(species, sample, hap)['genes_mapped'] = n_genes

        unmapped = os.path.join(liftoff_dir, f"{sample}.{hap}.unmapped_features.txt")
        if os.path.exists(unmapped):
            with open(unmapped) as fh:
                rec(species, sample, hap)['genes_unmapped'] = sum(1 for _ in fh)

        # --- TE content, from the RepeatMasker .tbl summary ---
        # Must match run_repeatmasker()'s naming exactly:
        #   outdir = TE_DIR/<assembly basename minus .fasta>
        #   tbl    = outdir/<assembly basename>.tbl
        assembly_stem = f"{species}_{sample}_{hap}.pseudo_chr"
        assembly_name = f"{assembly_stem}.fasta"
        tbl = os.path.join(te_dir, assembly_stem, f"{assembly_name}.tbl")
        if os.path.exists(tbl):
            with open(tbl) as fh:
                for line in fh:
                    m = te_re.search(line)
                    if m:
                        rec(species, sample, hap)['te_pct_masked'] = float(m.group(1))
                        break

# --- write per-haplotype table ---
per_hap_out = os.path.join(compare_dir, 'per_haplotype_summary.tsv')
fields = ['species', 'sample', 'haplotype', 'length_bp', 'genes_mapped', 'genes_unmapped', 'te_pct_masked']
with open(per_hap_out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    for r in sorted(records.values(), key=lambda r: (r['species'], r['sample'], r['haplotype'])):
        w.writerow(r)

# --- aggregate by species (mean +/- sd, skipping missing values) ---
by_species = {}
for r in records.values():
    by_species.setdefault(r['species'], []).append(r)

species_out = os.path.join(compare_dir, 'species_comparison_summary.tsv')
with open(species_out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['species', 'n_haplotypes',
                'mean_length_bp', 'sd_length_bp',
                'mean_genes_mapped', 'sd_genes_mapped',
                'mean_genes_unmapped', 'sd_genes_unmapped',
                'mean_te_pct_masked', 'sd_te_pct_masked'])
    for species, recs in sorted(by_species.items()):
        def col(name):
            vals = [r[name] for r in recs if r[name] is not None]
            if not vals:
                return ('NA', 'NA')
            mean = statistics.mean(vals)
            sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
            return (f'{mean:.2f}', f'{sd:.2f}')
        length_m, length_sd = col('length_bp')
        gm_m, gm_sd = col('genes_mapped')
        gu_m, gu_sd = col('genes_unmapped')
        te_m, te_sd = col('te_pct_masked')
        w.writerow([species, len(recs), length_m, length_sd, gm_m, gm_sd, gu_m, gu_sd, te_m, te_sd])

print(f"  Wrote {per_hap_out}")
print(f"  Wrote {species_out}")
PYEOF

echo ""
echo ">>> Pangenome + comparative analysis complete."
echo "    Pangenome graph        : ${GFA}"
echo "    Per-haplotype summary  : ${COMPARE_DIR}/per_haplotype_summary.tsv"
echo "    Species comparison     : ${COMPARE_DIR}/species_comparison_summary.tsv"
echo "    Core/accessory growth  : ${COMPARE_DIR}/pangenome_histgrowth_{all,<SPECIES>}.tsv"
echo ">>> End time: $(date)"
