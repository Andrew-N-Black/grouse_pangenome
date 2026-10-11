#!/bin/bash
# =============================================================================
# PANGENOME INPUT PREP 02 — automated assembly QC, one array task per sample
#
# Runs on exactly the FASTAs that will go into the graph (step 01 output), so
# the QC describes the pangenome inputs, not an earlier version of them.
#
# Per sample (both haplotypes together):
#   [1] HiFi reads -> FASTA once (or reuse adapter-filtered reads if present)
#   [2] meryl k-mer database of the reads (k=21; ~1.1 Gb genome)
#   [3] Merqury on hap1 + hap2:
#         per-haplotype QV (consensus accuracy), k-mer completeness, and the
#         copy-number spectra plots, which show false duplication (read-1x
#         k-mers present twice in the assembly) and collapses
#   [4] HiFi depth: reads mapped to hap1+hap2 TOGETHER (so each read can go to
#       the haplotype it came from), depth in 100 kb windows (mosdepth), and
#       windows flagged at < 0.5x or > 1.5x the haplotype's median depth.
#         low  -> possible false duplication (reads split between two copies)
#                 or poorly supported sequence
#         high -> possible collapsed repeat or duplication
#       Flagged windows go to a BED for review in a browser / PretextView.
#   [5] optional: NCBI FCS-adaptor screen (set FCS_ADAPTOR_SCRIPT + FCS_ADAPTOR_SIF)
#
# NOT here: NCBI FCS-GX (contamination). Its database is ~470 GB and wants a
# ~512 GB-memory node, so it belongs in its own job; it is also required by NCBI
# at submission, which is the natural time to run it.
#
# USAGE (from the grouse_pangenome repo root, so assembly_manifest.tsv is found)
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%6 02_assembly_qc_array.sh
#   sbatch --array=7 02_assembly_qc_array.sh             # one sample
#   then: bash 02b_collate_qc.sh
#
# Resume-safe: each stage is skipped when its output exists.
#
# MOSDEPTH / MERYL: Gautschi has no mosdepth module, and its merqury module
# does not put meryl on PATH. For either tool missing from the modules the
# script uses (and, once, downloads) the official release binary into
# ${ASM_DIR}/tools/<tool>-<version>/. Do the download on a LOGIN node first,
# since compute nodes may have no outbound network:
#   INSTALL_ONLY=true bash 02_assembly_qc_array.sh
# =============================================================================
#SBATCH --job-name=pg02_asm_qc
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
# 128G: meryl counting is the peak (bounded by MERYL_MEM below); minimap2
# map-hifi against a 2.2 Gb dual-haplotype target needs ~20-30 GB.
#SBATCH --mem=128G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

show_tsv() { if command -v column >/dev/null 2>&1; then column -t -s$'\t' "$1"; else cat "$1"; fi; }

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
ASM_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
PG_WORK="${PG_WORK:-${ASM_DIR}}"     # shared PROJECT_DIR, as steps 01, 03, 04
# This repo's copy of the manifest (sample, species, hifi_bams, ...), as in step 03
MANIFEST="${MANIFEST:-${SLURM_SUBMIT_DIR:-$PWD}/assembly_manifest.tsv}"

IN_DIR="${PG_WORK}/pangenome_inputs"   # step 01 output
QC_DIR="${PG_WORK}/pangenome_qc"

# Adapter-filtered HiFi reads, if they were kept. A glob with SAMPLE in it;
# leave empty to convert the manifest's HiFi BAMs instead.
#   e.g. HIFI_READS_GLOB='${ASM_DIR}/hifi_filt/${SAMPLE}*.filt.fastq.gz'
HIFI_READS_GLOB="${HIFI_READS_GLOB:-}"

K=21                    # Merqury best_k.sh for ~1.1 Gb gives ~20.9
MERYL_MEM=96            # GB
WINDOW=100000           # depth window (bp)
LOW_FRAC=0.5            # flag windows below this x median depth
HIGH_FRAC=1.5           # flag windows above this x median depth
MIN_WIN_BP=50000        # ignore end-of-sequence windows shorter than this
KEEP_READS="${KEEP_READS:-false}"   # keep the converted reads FASTA afterwards?
KEEP_BAM="${KEEP_BAM:-false}"       # keep the read-to-assembly BAM afterwards?

# Optional NCBI FCS-adaptor (https://github.com/ncbi/fcs). Both must be set.
FCS_ADAPTOR_SCRIPT="${FCS_ADAPTOR_SCRIPT:-}"   # path to run_fcsadaptor.sh
FCS_ADAPTOR_SIF="${FCS_ADAPTOR_SIF:-}"         # path to fcs-adaptor.sif

# Static mosdepth release binary, used when no module provides mosdepth.
MOSDEPTH_VERSION="0.3.14"
MOSDEPTH_DIR="${ASM_DIR}/tools/mosdepth-${MOSDEPTH_VERSION}"
MOSDEPTH_URL="https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth"

# meryl release binary, used when no module provides meryl. Pinned to 1.3 to
# match merqury/1.3, so the k-mer database is written by the same meryl
# version merqury expects to read.
MERYL_VERSION="1.3"
MERYL_DIR="${ASM_DIR}/tools/meryl-${MERYL_VERSION}"
MERYL_URL="https://github.com/marbl/meryl/releases/download/v${MERYL_VERSION}/meryl-${MERYL_VERSION}.Linux-amd64.tar.xz"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# =============================================================================
# ENVIRONMENT
# =============================================================================
set +u
module unload anaconda >/dev/null 2>&1 || true
module load biocontainers  >/dev/null 2>&1 || true
module load samtools/1.22.1 >/dev/null 2>&1 || true
module load minimap2       >/dev/null 2>&1 || true
module load merqury/1.3    >/dev/null 2>&1 || true   # pinned: the version on Gautschi
module load meryl          >/dev/null 2>&1 || true   # may not exist; see below
module load mosdepth       >/dev/null 2>&1 || true   # absent on Gautschi; see below
set -u
unset LD_PRELOAD    # RCAC XALT tracker; can kill subshells on some nodes

# ---- mosdepth: module if present, else the cached static binary, else fetch it
if command -v mosdepth >/dev/null 2>&1; then
    MOSDEPTH="$(command -v mosdepth)"
else
    MOSDEPTH="${MOSDEPTH_DIR}/mosdepth"
    if [[ ! -x "$MOSDEPTH" ]]; then
        echo ">>> mosdepth not found — downloading static v${MOSDEPTH_VERSION} to ${MOSDEPTH_DIR}"
        mkdir -p "$MOSDEPTH_DIR"
        if command -v wget >/dev/null 2>&1; then
            wget -q -O "${MOSDEPTH}.part" "$MOSDEPTH_URL" || rm -f "${MOSDEPTH}.part"
        else
            curl -fsSL -o "${MOSDEPTH}.part" "$MOSDEPTH_URL" || rm -f "${MOSDEPTH}.part"
        fi
        if [[ -s "${MOSDEPTH}.part" ]]; then
            chmod +x "${MOSDEPTH}.part" && mv "${MOSDEPTH}.part" "$MOSDEPTH"
        fi
    fi
    if [[ ! -x "$MOSDEPTH" ]] || ! "$MOSDEPTH" --version >/dev/null 2>&1; then
        echo "ERROR: could not obtain a working mosdepth at ${MOSDEPTH}"
        echo "       On a login node run:  INSTALL_ONLY=true bash 02_assembly_qc_array.sh"
        echo "       or download ${MOSDEPTH_URL} there yourself and chmod +x it."
        exit 1
    fi
fi
# ---- meryl: module if present, else the cached release binary, else fetch it
if command -v meryl >/dev/null 2>&1; then
    MERYL="$(command -v meryl)"
else
    MERYL="${MERYL_DIR}/bin/meryl"
    if [[ ! -x "$MERYL" ]]; then
        echo ">>> meryl not found — downloading v${MERYL_VERSION} to ${MERYL_DIR}"
        mkdir -p "$MERYL_DIR"
        TARBALL="${MERYL_DIR}/meryl.tar.xz"
        if command -v wget >/dev/null 2>&1; then
            wget -q -O "${TARBALL}.part" "$MERYL_URL" || rm -f "${TARBALL}.part"
        else
            curl -fsSL -o "${TARBALL}.part" "$MERYL_URL" || rm -f "${TARBALL}.part"
        fi
        if [[ -s "${TARBALL}.part" ]]; then
            # the archive holds meryl-<ver>/bin/... ; flatten into MERYL_DIR
            tar -xJf "${TARBALL}.part" -C "$MERYL_DIR" --strip-components=1 && rm -f "${TARBALL}.part"
        fi
    fi
    if [[ ! -x "$MERYL" ]] || ! "$MERYL" --version >/dev/null 2>&1; then
        echo "ERROR: could not obtain a working meryl at ${MERYL}"
        echo "       On a login node run:  INSTALL_ONLY=true bash 02_assembly_qc_array.sh"
        echo "       or unpack ${MERYL_URL} into ${MERYL_DIR} yourself."
        exit 1
    fi
    # merqury.sh calls meryl itself; make this one visible to it too
    export PATH="${MERYL_DIR}/bin:${PATH}"
fi

if [[ "${INSTALL_ONLY:-false}" == true ]]; then
    echo ">>> mosdepth ready: ${MOSDEPTH} ($("$MOSDEPTH" --version 2>&1 | head -1))"
    echo ">>> meryl ready   : ${MERYL} ($("$MERYL" --version 2>&1 | head -1))"
    echo ">>> merqury.sh    : $(command -v merqury.sh || echo 'NOT FOUND — check: ml spider merqury')"
    exit 0
fi

MISSING=()
for B in samtools minimap2 merqury.sh; do
    command -v "$B" >/dev/null 2>&1 || MISSING+=("$B")
done
if (( ${#MISSING[@]} )); then
    echo "ERROR: not on PATH after module load: ${MISSING[*]}"
    echo "       check names with:  ml spider merqury ; ml spider minimap2"
    exit 1
fi
echo ">>> samtools $(command -v samtools)"
echo ">>> minimap2 $(minimap2 --version 2>/dev/null || true)"
echo ">>> meryl    ${MERYL} ($("$MERYL" --version 2>&1 | head -1 || true))"
echo ">>> merqury  $(command -v merqury.sh)"
echo ">>> mosdepth ${MOSDEPTH} ($("$MOSDEPTH" --version 2>&1 | head -1 || true))"

# =============================================================================
# SAMPLE FOR THIS TASK
# =============================================================================
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: submit with sbatch --array=..."; exit 1; }
[[ -s "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST"; exit 1; }

mapfile -t ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2 | grep -v '^[[:space:]]*$')
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"
[[ -n "$LINE" ]] || { echo "ERROR: no manifest row ${SLURM_ARRAY_TASK_ID} (${#ROWS[@]} rows)"; exit 1; }
IFS=$'\t' read -r SAMPLE SPECIES HIFI_BAMS _R1 _R2 _ONT <<< "$LINE"

H1=$(ls "${IN_DIR}"/*_"${SAMPLE}"_hap1.pg_input.fasta 2>/dev/null | head -1 || true)
H2=$(ls "${IN_DIR}"/*_"${SAMPLE}"_hap2.pg_input.fasta 2>/dev/null | head -1 || true)
for F in "$H1" "$H2"; do
    [[ -n "$F" && -s "$F" ]] || { echo "ERROR: missing pg_input FASTA for ${SAMPLE} in ${IN_DIR} — run step 01"; exit 1; }
done
P1=$(basename "$H1" .pg_input.fasta)
P2=$(basename "$H2" .pg_input.fasta)

WORK="${QC_DIR}/${SAMPLE}"
mkdir -p "$WORK" logs
echo ""
echo ">>> task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE} (${SPECIES})  threads=${THREADS}"
echo ">>> ${P1}  ${P2}"
echo ">>> work: ${WORK}"

# =============================================================================
# [1] READS
# =============================================================================
READS=()
if [[ -n "$HIFI_READS_GLOB" ]]; then
    PATTERN=$(eval echo "$HIFI_READS_GLOB")
    # shellcheck disable=SC2206
    READS=( $PATTERN )
    [[ -e "${READS[0]:-}" ]] || READS=()
fi
CONVERTED="${WORK}/${SAMPLE}.hifi.fa"
if (( ${#READS[@]} )); then
    echo ">>> [1] using filtered reads: ${READS[*]}"
else
    if [[ -s "$CONVERTED" && ! -e "${CONVERTED}.part" ]]; then
        echo ">>> [1] reads FASTA already present"
    else
        IFS=', ' read -r -a BAMS <<< "$HIFI_BAMS"
        (( ${#BAMS[@]} )) || { echo "ERROR: no HiFi BAMs in manifest for ${SAMPLE}"; exit 1; }
        for B in "${BAMS[@]}"; do [[ -s "$B" ]] || { echo "ERROR: HiFi BAM not found: $B"; exit 1; }; done
        echo ">>> [1] converting ${#BAMS[@]} HiFi BAM(s) to FASTA"
        : > "${CONVERTED}.part"
        for B in "${BAMS[@]}"; do
            samtools fasta -@ "$THREADS" "$B" >> "${CONVERTED}.part"
        done
        mv "${CONVERTED}.part" "$CONVERTED"
    fi
    READS=("$CONVERTED")
fi

# =============================================================================
# [2] MERYL
# =============================================================================
DB="${WORK}/${SAMPLE}.k${K}.meryl"
if [[ -d "$DB" && -e "${DB}/merylIndex" ]]; then
    echo ">>> [2] meryl db present — skipping"
else
    echo ">>> [2] meryl count k=${K}"
    rm -rf "$DB"
    "$MERYL" count k="$K" threads="$THREADS" memory="$MERYL_MEM" "${READS[@]}" output "$DB"
fi

# =============================================================================
# [3] MERQURY (run inside its own directory; it writes to the CWD)
# =============================================================================
MQ="${WORK}/merqury"
mkdir -p "$MQ"
if [[ -s "${MQ}/${SAMPLE}.qv" && -s "${MQ}/${SAMPLE}.completeness.stats" ]]; then
    echo ">>> [3] merqury done — skipping"
else
    echo ">>> [3] merqury"
    ln -sf "$H1" "${MQ}/${P1}.fasta"
    ln -sf "$H2" "${MQ}/${P2}.fasta"
    ( cd "$MQ" && OMP_NUM_THREADS="$THREADS" merqury.sh "$DB" "${P1}.fasta" "${P2}.fasta" "$SAMPLE" )
    [[ -s "${MQ}/${SAMPLE}.qv" ]] || { echo "ERROR: merqury produced no ${SAMPLE}.qv — see ${MQ}/logs"; exit 1; }
fi
echo "    QV:";            sed 's/^/      /' "${MQ}/${SAMPLE}.qv"
echo "    completeness:";  sed 's/^/      /' "${MQ}/${SAMPLE}.completeness.stats"

# =============================================================================
# [4] HIFI DEPTH against both haplotypes at once
# =============================================================================
DUAL="${WORK}/${SAMPLE}.dual.fa"
BAM="${WORK}/${SAMPLE}.hifi_vs_dual.bam"
MD="${WORK}/${SAMPLE}.depth"
if [[ -s "${MD}.regions.bed.gz" ]]; then
    echo ">>> [4] depth done — skipping"
else
    echo ">>> [4] HiFi depth (minimap2 map-hifi -> mosdepth ${WINDOW} bp windows)"
    # prefix names with the haplotype so both copies coexist in one reference
    { sed "s/^>/>h1|/" "$H1"; sed "s/^>/>h2|/" "$H2"; } > "$DUAL"
    if [[ ! -s "$BAM" || ! -s "${BAM}.bai" ]]; then
        minimap2 -ax map-hifi -t "$THREADS" --secondary=no "$DUAL" "${READS[@]}" \
            | samtools sort -@ 8 -m 2G -o "${BAM}.part" -
        mv "${BAM}.part" "$BAM"
        samtools index -@ 8 "$BAM"
    fi
    # -n: no per-base output; -x: fast mode; MAPQ 0 kept on purpose — reads
    # that cannot choose between two copies are exactly the false-duplication
    # signal, and dropping them would hide it.
    "$MOSDEPTH" -n -x -t 4 --by "$WINDOW" "$MD" "$BAM"
fi

FLAGS="${WORK}/${SAMPLE}.depth_flags.bed"
DSUM="${WORK}/${SAMPLE}.depth_summary.tsv"
WIN="${WORK}/${SAMPLE}.windows.tsv"
# hap  seq  beg  end  depth   (full-size windows only)
zcat "${MD}.regions.bed.gz" \
  | awk -F'\t' -v OFS='\t' -v minw="$MIN_WIN_BP" '
      $3 - $2 >= minw { i = index($1, "|"); print substr($1, 1, i - 1), substr($1, i + 1), $2, $3, $4 }' \
  > "$WIN"
printf 'haplotype\twindows\tmedian_depth\tlow_windows\tlow_bp_pct\thigh_windows\thigh_bp_pct\n' > "$DSUM"
: > "$FLAGS"
for HP in h1 h2; do
    NAME=$([[ "$HP" == h1 ]] && echo "$P1" || echo "$P2")
    MED=$(awk -F'\t' -v h="$HP" '$1 == h {print $5}' "$WIN" | sort -g \
          | awk '{v[NR] = $1} END {if (NR) print v[int((NR + 1) / 2)]; else print 0}')
    awk -F'\t' -v OFS='\t' -v h="$HP" -v name="$NAME" -v med="$MED" \
        -v lo="$LOW_FRAC" -v hi="$HIGH_FRAC" -v flags="$FLAGS" -v dsum="$DSUM" '
        $1 == h { n++; L = $4 - $3; tot += L
                  if ($5 < lo * med)      { nl++; lbp += L; print $2, $3, $4, name ":low:" $5  >> flags }
                  else if ($5 > hi * med) { nh++; hbp += L; print $2, $3, $4, name ":high:" $5 >> flags } }
        END { if (n) printf "%s\t%d\t%.1f\t%d\t%.2f\t%d\t%.2f\n", name, n, med, nl,
                            100 * lbp / tot, nh, 100 * hbp / tot >> dsum }' "$WIN"
done
rm -f "$WIN"
echo "    depth:"; show_tsv "$DSUM" | sed 's/^/      /'

# =============================================================================
# [5] OPTIONAL: FCS-adaptor
# =============================================================================
if [[ -n "$FCS_ADAPTOR_SCRIPT" && -n "$FCS_ADAPTOR_SIF" ]]; then
    for F in "$H1" "$H2"; do
        P=$(basename "$F" .pg_input.fasta)
        OUT="${WORK}/fcs_adaptor/${P}"
        if [[ -s "${OUT}/fcs_adaptor_report.txt" ]]; then
            echo ">>> [5] FCS-adaptor ${P} done — skipping"; continue
        fi
        mkdir -p "$OUT"
        echo ">>> [5] FCS-adaptor ${P}"
        bash "$FCS_ADAPTOR_SCRIPT" --fasta-input "$F" --output-dir "$OUT" --euk \
            --container-engine singularity --image "$FCS_ADAPTOR_SIF"
        echo "      $(grep -vc '^#' "${OUT}/fcs_adaptor_report.txt" || true) adaptor hit(s)"
    done
else
    echo ">>> [5] FCS-adaptor skipped (set FCS_ADAPTOR_SCRIPT and FCS_ADAPTOR_SIF to run it)"
fi

# =============================================================================
# CLEAN UP
# =============================================================================
rm -f "$DUAL"
[[ "$KEEP_BAM"   == true ]] || rm -f "$BAM" "${BAM}.bai"
[[ "$KEEP_READS" == true ]] || rm -f "$CONVERTED"

echo ""
echo ">>> ${SAMPLE} complete: $(date)"
echo "    merqury : ${MQ}/   (spectra-cn plots: ${SAMPLE}.*.spectra-cn.*.png)"
echo "    depth   : ${DSUM}"
echo "    flags   : ${FLAGS}"
