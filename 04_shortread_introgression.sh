#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: SHORT-READ MAPPING TO PANGENOME + INTROGRESSION
# Step 04 — requires 03_pangenome_build.sh to have completed (consumes its
# *.final.gfa pangenome graph).
#
# Pipeline:
#   1. Convert the PGGB graph to a giraffe-ready GBZ + indexes.
#   2. Per short-read sample: vg giraffe (map) -> vg pack (coverage) ->
#      vg call (genotype), looped over each contig of a chosen reference
#      haplotype, then concatenated into one whole-genome VCF per sample.
#   3. bcftools merge across all samples into one joint, multi-sample VCF.
#   4. Dsuite Dtrios computes genome-wide D-statistics (ABBA-BABA) for the
#      3 focal species + outgroup; Dsuite Dinvestigate then computes
#      windowed D/f_d statistics for whichever pair Dtrios' data-driven
#      BBAA pattern identifies as showing excess allele sharing.
#
# WHY AN OUTGROUP IS REQUIRED: D-statistics need a 4th, more distantly
# related taxon to polarize ancestral vs. derived alleles. Without one,
# "quantify introgression" isn't a well-defined computation — it's a guess.
# See shortread_manifest.tsv; at least one sample must be labeled "Outgroup".
#
# WHY NO HARDCODED SPECIES TREE: with only 3 ingroup species there is
# exactly one possible trio, and we don't have independent knowledge of
# these species' true phylogeny or which pairs are geographically sympatric.
# Rather than fabricate either, Dsuite's default (untreed) mode is used: it
# infers the best-supported P1/P2/P3 arrangement directly from allele-sharing
# patterns in the data (the "BBAA" pattern), and that inferred arrangement's
# P2/P3 pair is the one showing the introgression/hybridization signal.
# Cross-reference the result against your own knowledge of species ranges to
# judge whether the statistically supported pair is also geographically
# sympatric (biologically plausible) or not.
#
# VERIFICATION NOTE: `vg gbwt -G`/`vg autoindex`/`vg call`/`vg giraffe`/
# `Dsuite Dtrios`/`Dsuite Dinvestigate` flags below were checked against
# current --help/README text. `vg pack` and `vg paths` flags are from
# stable, long-established usage rather than a freshly-quoted --help — if
# either errors, check `vg pack --help` / `vg paths --help` first.
#
# INPUT MANIFEST (tab-separated, header required, see shortread_manifest.tsv):
#   sample_id  species  fastq_r1  fastq_r2
#     - species: one of the 3 focal species labels from assembly_manifest.tsv,
#       or "Outgroup"
#
# USAGE:
#   sbatch 04_shortread_introgression.sh
# =============================================================================
#SBATCH --job-name=grouse_introgression
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
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
ml bcftools
ml dsuite/0.5.r44

# vg has no module on Gautschi: use the Apptainer image from
# accessory/pull_containers.sh unless a vg is already on PATH.
VG_SIF="${VG_SIF:-${CLUSTER_SCRATCH}/GROUSE/grouse_asm/containers/vg.sif}"
APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
export APPTAINERENV_LD_PRELOAD="" SINGULARITYENV_LD_PRELOAD=""
if ! command -v vg >/dev/null 2>&1; then
    [[ -n "$APPTAINER_BIN" && -s "$VG_SIF" ]] || { echo "ERROR: no vg on PATH and no ${VG_SIF} (run accessory/pull_containers.sh)"; exit 1; }
    vg() { "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$VG_SIF" vg "$@"; }
    export -f vg; export APPTAINER_BIN VG_SIF CLUSTER_SCRATCH
fi

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
ASSEMBLY_MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"
SHORTREAD_MANIFEST="${SLURM_SUBMIT_DIR}/shortread_manifest.tsv"

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
PANGENOME_DIR="${PROJECT_DIR}/pangenome"   # input: from 03_pangenome_build.sh
GRAPH_DIR="${PROJECT_DIR}/graph_index"
MAP_DIR="${PROJECT_DIR}/mapped"
CALL_DIR="${PROJECT_DIR}/called"
JOINT_DIR="${PROJECT_DIR}/joint_vcf"
INTROGRESSION_DIR="${PROJECT_DIR}/introgression"

THREADS=$SLURM_CPUS_PER_TASK
MAP_PARALLEL_JOBS=4
MAP_THREADS_PER_JOB=$(( THREADS / MAP_PARALLEL_JOBS > 0 ? THREADS / MAP_PARALLEL_JOBS : 1 ))

# Minimum mapping quality for pack/genotyping support (standard vg default)
PACK_MIN_MAPQ=5

# Dinvestigate windowed-statistics parameters (SNP count, step) — Dsuite's
# own documented default.
DINVESTIGATE_WINDOW="50,25"

mkdir -p logs "$GRAPH_DIR" "$MAP_DIR" "$CALL_DIR" "$JOINT_DIR" "$INTROGRESSION_DIR"

echo ">>> 04_shortread_introgression.sh"
echo ">>> Node      : $(hostname)"
echo ">>> CPUs      : ${THREADS}"
echo ">>> Start time: $(date)"

if [[ ! -f "$ASSEMBLY_MANIFEST" ]]; then
    echo "ERROR: Assembly manifest not found: ${ASSEMBLY_MANIFEST}"
    exit 1
fi
if [[ ! -f "$SHORTREAD_MANIFEST" ]]; then
    echo "ERROR: Short-read manifest not found: ${SHORTREAD_MANIFEST}"
    exit 1
fi

# Per-chromosome builds (03 with CHROM=...) are merged by 03's STAGE=summary
# into pangenome/merged/; a whole-genome build sits directly in pangenome/.
if [[ -s "${PANGENOME_DIR}/merged/pangenome.final.gfa" ]]; then
    GFA="${PANGENOME_DIR}/merged/pangenome.final.gfa"
else
    GFA=$(find "$PANGENOME_DIR" -maxdepth 1 -name "*.final.gfa" 2>/dev/null | head -n1)
fi
if [[ -z "$GFA" ]]; then
    echo "ERROR: No pangenome graph found in ${PANGENOME_DIR}"
    echo "Run 03_pangenome_build.sh first."
    exit 1
fi

sanitize() {
    # printf, not echo: echo's trailing newline would become a trailing '_'
    printf '%s' "$1" | tr -c 'A-Za-z0-9_-' '_'
}

mapfile -t ASM_ROWS < <(grep -v '^#' "$ASSEMBLY_MANIFEST" | tail -n +2)
mapfile -t SR_ROWS < <(grep -v '^#' "$SHORTREAD_MANIFEST" | tail -n +2)

if [[ "${#SR_ROWS[@]}" -eq 0 ]]; then
    echo "ERROR: No samples found in ${SHORTREAD_MANIFEST}"
    exit 1
fi

if ! printf '%s\n' "${SR_ROWS[@]}" | cut -f2 | grep -qx "Outgroup"; then
    echo "ERROR: No sample labeled 'Outgroup' in ${SHORTREAD_MANIFEST}."
    echo "D-statistics require an outgroup — see the header comment in this script."
    exit 1
fi

# Reference haplotype for VCF coordinates: defaults to the first sample in
# assembly_manifest.tsv (hap1), so this always points at a real path in the
# graph without guessing a sample ID. Override REF_HAP_TAG below if you'd
# rather anchor calling on a specific higher-quality individual.
IFS=$'\t' read -r REF_SAMPLE REF_SPECIES _ _ _ _ <<< "${ASM_ROWS[0]}"
REF_SPECIES="${REF_SPECIES//[[:space:]]/}"   # manifest species may carry stray spaces ("LEPC ")
# PanSN splits a path name into sample#haplotype#contig, and `vg paths -S`
# selects on the SAMPLE field alone — passing "sample#hap1" to it matches
# nothing. So keep the two apart: select the sample, then filter to the
# haplotype by name. REF_HAP_TAG is only used for display and error messages.
REF_SAMPLE_TAG="$(sanitize "$REF_SPECIES").$(sanitize "$REF_SAMPLE")"
REF_HAP="1"          # PanSN haplotype field is an integer (step 03 writes hap1 as #1#)
REF_HAP_TAG="${REF_SAMPLE_TAG}#${REF_HAP}"
echo ">>> Reference haplotype for VCF coordinates: ${REF_HAP_TAG}"

# =============================================================================
# STEP 1: Convert the PGGB GFA to a giraffe-ready GBZ, then build the
# distance/minimizer/zipcode indexes vg giraffe needs.
#
# `vg autoindex --workflow giraffe` cannot reliably read a raw PGGB GFA
# directly (PGGB writes all paths as P-lines; vg's GFA loader for autoindex
# expects W-lines for per-haplotype paths — see
# https://github.com/vgteam/vg/issues/4302). `vg gbwt -G` handles this path
# structure correctly and produces a GBZ, which autoindex then reads cleanly.
# =============================================================================
GBZ="${GRAPH_DIR}/pangenome.gbz"

if [[ ! -f "$GBZ" ]]; then
    echo ">>> Step 1a: Building GBZ from PGGB GFA (vg gbwt)"
    vg gbwt -G "$GFA" -g "$GBZ" -t "$THREADS"
else
    echo ">>> Step 1a: GBZ already exists — skipping: ${GBZ}"
fi

GIRAFFE_INDEX_PREFIX="${GRAPH_DIR}/pangenome"
if [[ ! -f "${GIRAFFE_INDEX_PREFIX}.dist" ]]; then
    echo ">>> Step 1b: Building giraffe indexes (vg autoindex)"
    vg autoindex --workflow giraffe -p "$GIRAFFE_INDEX_PREFIX" -G "$GBZ" -t "$THREADS"
else
    echo ">>> Step 1b: Giraffe indexes already exist — skipping"
fi

GIRAFFE_GBZ="${GIRAFFE_INDEX_PREFIX}.giraffe.gbz"
DIST_IDX="${GIRAFFE_INDEX_PREFIX}.dist"
MIN_IDX=$(find "$GRAPH_DIR" -name "pangenome*.min" 2>/dev/null | head -n1)
ZIP_IDX=$(find "$GRAPH_DIR" -name "pangenome*.zipcodes" 2>/dev/null | head -n1)

for f in "$GIRAFFE_GBZ" "$DIST_IDX" "$MIN_IDX" "$ZIP_IDX"; do
    if [[ -z "$f" || ! -f "$f" ]]; then
        echo "ERROR: Expected giraffe index file missing under ${GRAPH_DIR}"
        echo "Check 'vg autoindex --workflow giraffe' output above for errors."
        exit 1
    fi
done

# Contigs of the reference haplotype, used to loop `vg call` per chromosome.
REF_CONTIGS_FILE="${GRAPH_DIR}/ref_contigs.txt"
vg paths -x "$GIRAFFE_GBZ" -S "$REF_SAMPLE_TAG" -L \
    | grep -F "#${REF_HAP}#" > "${REF_CONTIGS_FILE}.part" || true
mv "${REF_CONTIGS_FILE}.part" "$REF_CONTIGS_FILE"
N_REF_CONTIGS=$(wc -l < "$REF_CONTIGS_FILE")
echo ">>> Reference haplotype has ${N_REF_CONTIGS} contigs/paths for calling"

if [[ "$N_REF_CONTIGS" -eq 0 ]]; then
    echo "ERROR: No paths found for REF_HAP_TAG=${REF_HAP_TAG} in the graph."
    echo "Check that this sample/haplotype was actually included when"
    echo "03_pangenome_build.sh built the combined FASTA."
    exit 1
fi

# =============================================================================
# STEP 2: Per-sample mapping (giraffe) -> coverage (pack) -> genotyping
# (call, looped per reference contig, then concatenated). Parallelized
# across samples with xargs -P (same pattern as the RepeatMasker loop in
# 03_pangenome_build.sh).
# =============================================================================
process_sample() {
    local line="$1"
    local sample species r1 r2
    IFS=$'\t' read -r sample species r1 r2 <<< "$line"

    echo ">>> [${sample}] Step 2a: vg giraffe mapping"
    local gam="${MAP_DIR}/${sample}.gam"
    vg giraffe -t "$MAP_THREADS_PER_JOB" \
        -Z "$GIRAFFE_GBZ" -d "$DIST_IDX" -z "$ZIP_IDX" -m "$MIN_IDX" \
        -f "$r1" -f "$r2" -b default > "$gam"

    echo ">>> [${sample}] Step 2b: vg pack (coverage, MAPQ >= ${PACK_MIN_MAPQ})"
    local pack="${MAP_DIR}/${sample}.pack"
    vg pack -x "$GIRAFFE_GBZ" -g "$gam" -o "$pack" -Q "$PACK_MIN_MAPQ" -t "$MAP_THREADS_PER_JOB"

    echo ">>> [${sample}] Step 2c: vg call (per reference contig)"
    local sample_vcf_parts=()
    while IFS= read -r contig; do
        local part="${CALL_DIR}/${sample}.${contig}.vcf.gz"
        vg call "$GIRAFFE_GBZ" -k "$pack" -s "$sample" -p "$contig" -t "$MAP_THREADS_PER_JOB" \
            | bgzip -c > "$part"
        tabix -p vcf "$part"
        sample_vcf_parts+=("$part")
    done < "$REF_CONTIGS_FILE"

    echo ">>> [${sample}] Step 2d: Concatenating per-contig VCFs"
    local sample_vcf="${CALL_DIR}/${sample}.vcf.gz"
    bcftools concat -Oz -o "$sample_vcf" "${sample_vcf_parts[@]}"
    tabix -p vcf "$sample_vcf"

    echo ">>> [${sample}] complete: ${sample_vcf}"
}
export -f process_sample
export MAP_DIR CALL_DIR GIRAFFE_GBZ DIST_IDX ZIP_IDX MIN_IDX REF_CONTIGS_FILE
export MAP_THREADS_PER_JOB PACK_MIN_MAPQ

printf '%s\n' "${SR_ROWS[@]}" | xargs -d '\n' -I{} -P "$MAP_PARALLEL_JOBS" bash -c 'process_sample "$@"' _ {}

# =============================================================================
# STEP 3: Joint multi-sample VCF (bcftools merge)
# =============================================================================
echo ">>> Step 3: bcftools merge across all samples"

JOINT_VCF="${JOINT_DIR}/joint.vcf.gz"
mapfile -t SAMPLE_VCFS < <(printf '%s\n' "${SR_ROWS[@]}" | cut -f1 | sed "s|.*|${CALL_DIR}/&.vcf.gz|")

bcftools merge -Oz -o "$JOINT_VCF" "${SAMPLE_VCFS[@]}"
tabix -p vcf "$JOINT_VCF"

echo "  Joint VCF: ${JOINT_VCF}"

# =============================================================================
# STEP 4: SETS.tsv for Dsuite — the manifest's sample/species columns already
# match Dsuite's required 2-column format (including the "Outgroup" label).
# =============================================================================
SETS_TSV="${INTROGRESSION_DIR}/SETS.tsv"
printf '%s\n' "${SR_ROWS[@]}" | cut -f1,2 > "$SETS_TSV"
echo ">>> Step 4: Wrote ${SETS_TSV}"

# =============================================================================
# STEP 5: Dsuite Dtrios — genome-wide D-statistics (ABBA-BABA), untreed
# (data-driven BBAA pattern; see header comment for why no tree is assumed).
# =============================================================================
echo ">>> Step 5: Dsuite Dtrios"

DTRIOS_PREFIX="${INTROGRESSION_DIR}/dtrios"
Dsuite Dtrios -o "$DTRIOS_PREFIX" "$JOINT_VCF" "$SETS_TSV"

BBAA_FILE="${DTRIOS_PREFIX}_BBAA.txt"
DMIN_FILE="${DTRIOS_PREFIX}_Dmin.txt"

if [[ ! -f "$BBAA_FILE" ]]; then
    echo "ERROR: Dsuite Dtrios did not produce expected output: ${BBAA_FILE}"
    exit 1
fi

echo "  BBAA (data-driven best arrangement): ${BBAA_FILE}"
echo "  Dmin (conservative bound, any arrangement): ${DMIN_FILE}"
echo "  --- BBAA result ---"
column -t "$BBAA_FILE"

# =============================================================================
# STEP 6: Dsuite Dinvestigate — windowed D/f_d statistics for the specific
# P1/P2/P3 arrangement Dtrios identified, localizing introgression signal
# across the genome. The P2/P3 pair in this output is the one showing
# excess allele sharing (the introgression/hybridization signal); cross-
# check against known species ranges to assess whether it's also sympatric.
# =============================================================================
echo ">>> Step 6: Dsuite Dinvestigate (windowed local statistics)"

TEST_TRIOS="${INTROGRESSION_DIR}/test_trios.txt"
awk 'NR==2 { print $1"\t"$2"\t"$3 }' "$BBAA_FILE" > "$TEST_TRIOS"

if [[ ! -s "$TEST_TRIOS" ]]; then
    echo "ERROR: Could not parse a P1/P2/P3 trio from ${BBAA_FILE}"
    exit 1
fi

echo "  Trio for windowed analysis (P1 P2 P3): $(cat "$TEST_TRIOS")"

Dsuite Dinvestigate -w "$DINVESTIGATE_WINDOW" "$JOINT_VCF" "$SETS_TSV" "$TEST_TRIOS"
mv ./*_localFstats_*.txt "$INTROGRESSION_DIR"/ 2>/dev/null || true

echo ""
echo ">>> Introgression analysis complete."
echo "    Joint VCF          : ${JOINT_VCF}"
echo "    Genome-wide D-stats : ${BBAA_FILE}"
echo "    Windowed local stats: ${INTROGRESSION_DIR}/*_localFstats_*.txt"
echo ">>> End time: $(date)"
