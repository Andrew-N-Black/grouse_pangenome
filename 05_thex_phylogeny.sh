#!/bin/bash
# =============================================================================
# STEP 05: WINDOWED PHYLOGENIES ACROSS THE 23 GROUSE GENOMES FOR TREE HOUSE EXPLORER
#
# Builds a Tree Viewer input for THEx (Tree House Explorer, Harris et al. 2022 MBE,
# https://github.com/harris-2374/THEx) from the 46 phased haplotype assemblies
# (23 birds x hap1/hap2) plus the rock ptarmigan as outgroup, and turns the
# Liftoff annotation into THEx's gene track. Runs on NEGISHI; the assemblies are
# pulled over from Gautschi first.
#
#   Gautschi  final_ptarmigan/<SP>_<S>_hap{1,2}.pseudo_chr.fasta(.fai) + .liftoff.gff3
#             ref/GCF_023343835.1_*_genomic.{fna,gff}  (rock ptarmigan)
#        |  rsync (stage: transfer, Negishi login node)
#        v
#   Negishi   1. pick-ref   one haplotype becomes the coordinate system (auto, see below)
#             2. project    array: minimap2 every haplotype + ptarmigan to it; write each
#                           one's bases at the reference coordinates (reference-projected
#                           alignment, insertions dropped, unaligned/duplicated = '-')
#             3. prep       per-chromosome MSAs -> thexb --minifastas (100 kb windows)
#                           -> thexb --trimal -> thexb --pw_filter (chunked in parallel)
#             4. iqtree     array: one IQ-TREE 2 ML tree + UFBoot per window
#             5. summarize  Tree Viewer .xlsx: trees rooted on ptarmigan, TopologyID =
#                           species-level topology, + per-window gene counts from the
#                           reference haplotype's Liftoff GFF3
#
# WHY A REFERENCE-PROJECTED ALIGNMENT
#   THExBuilder starts from one multiple alignment per chromosome in which every
#   sample has the same length. A whole-genome MSA of 47 genomes (Cactus, or the
#   step-03 PGGB graph) would also do, but projection onto one haplotype is what the
#   THEx authors used (Foley et al. 2021), needs nothing beyond minimap2, and is
#   insensitive to the RagTag misjoins in final_ptarmigan: a query's scaffolding does
#   not matter, only where its sequence aligns. Only the REFERENCE haplotype's own
#   chromosome structure matters, which is why pick-ref chooses it carefully.
#
# HOW THE REFERENCE HAPLOTYPE IS CHOSEN (override with REF_ID=...)
#   Among haplotypes with every chr_*, a full-length chr_Z and no chr_W, the one whose
#   chromosome lengths deviate least from the panel median (worst chromosome). A
#   misjoin such as the chr_23+chr_24 merger makes that haplotype's worst deviation
#   ~100%, so misjoined haplotypes are never picked. Table: work/ref_choice.tsv.
#   Windows are then in that haplotype's coordinates and named after the ptarmigan
#   chromosome numbering (chr_6 -> "chr-6" in THEx; THExBuilder swaps '_' for '-').
#
# WHAT THE TopologyID COLUMN MEANS
#   With 47 tips almost every window has its own exact topology, so THExBuilder's
#   topobinner (RF = 0 bins, all-against-all) would return thousands of singleton
#   bins and run for days. Instead each rooted tree is labelled by the species-level
#   relationship it shows: "((LEPC,GRPC),STGR)" when all three species are
#   reciprocally monophyletic, "... LEPC/GRPC unsorted" when the pair forms a clade but
#   its two species are mixed inside it, "... STGR paraphyletic", or
#   "No species pair forms a clade". THEx colours windows by this column. The full
#   47-tip trees are kept in NewickTree. Extra numeric columns THEx can plot:
#   CladeSupport (UFBoot of the clade that defines the label), Taxa, AlignedSites,
#   Genes (Liftoff genes overlapping the window).
#
# SEX CHROMOSOMES
#   chr_W and chr_MT are excluded. chr_Z is kept; a haplotype is left out of a
#   chromosome's alignment when < MIN_CHROM_COV of that chromosome is covered, which
#   removes the W-bearing haplotypes of the six ZW females from chr_Z only.
#
# A THExBuilder BUG THIS SCRIPT WORKS AROUND
#   thexb's pairwise filter computes "countmismatch + countmissing / denominator", so
#   its p-distance is a raw mismatch COUNT per sub-window, not a proportion, and the
#   default 0.15 cutoff means "any mismatch". The prep stage checks the installed
#   source and, if the bug is there, passes PW_MAX_PDIST x sub-window length (0.15 x
#   100 = 15) so the filter does what its documentation says. The outgroup is also in
#   the exclusion list: otherwise every ptarmigan sub-window is an outlier against the
#   ingroup median and gets masked, which drops nearly every window at the coverage
#   check.
#
# USAGE (Negishi, from the repo root so assembly_manifest.tsv is found)
#   bash 05_thex_phylogeny.sh transfer     # login node: pull assemblies from Gautschi
#   bash 05_thex_phylogeny.sh setup        # login node: THEx conda env (thexb, iqtree2, trimal)
#   bash 05_thex_phylogeny.sh pick-ref     # login node, seconds: shows the choice
#   bash 05_thex_phylogeny.sh submit       # submits stages 2-5 as one dependency chain
#   bash 05_thex_phylogeny.sh status       # progress counts
#   Any stage can be (re)submitted alone, e.g.
#     sbatch --export=ALL -c 128 --mem=240G -t 1-00:00:00 05_thex_phylogeny.sh prep
#   Every stage is resume-safe: finished haplotypes, chromosomes and windows are skipped.
#   Then copy work/thex_inputs/ to the laptop and run `thex` (see the end of this header).
#
# TRANSFER: needs passwordless ssh from Negishi to Gautschi. RCAC home directories are
#   shared between clusters, so once is enough:
#     ssh-keygen -t ed25519 && cat ~/.ssh/id_ed25519.pub >> ~/.ssh/authorized_keys
#   If that is not allowed, move final_ptarmigan/ with Globus (Purdue endpoints) or via
#   /depot/<lab>, into ${IN_DIR}, then run `bash 05_thex_phylogeny.sh transfer` anyway:
#   with the files present it only verifies and indexes them.
#
# VIEWING (laptop): conda create -n thex_env -c ajharris_2374 -c bioconda -c conda-forge thex
#   conda activate thex_env && thex   -> open http://127.0.0.1:8050 -> Tree Viewer, and load
#   grouse_TreeViewer_input.xlsx, chromosome_lengths.bed and reference_genes.gff.
# =============================================================================
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
# Literal address: SBATCH lines are not shell-expanded.
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail
# RCAC injects XALT through LD_PRELOAD; it fails on some nodes and kills subshells.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES (all can be overridden from the environment)
# =============================================================================
REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
MANIFEST="${MANIFEST:-${REPO_DIR}/assembly_manifest.tsv}"
SCRIPT="${REPO_DIR}/05_thex_phylogeny.sh"
PYHELP="${REPO_DIR}/05_thex_phylogeny.py"

# ---- source (Gautschi) ----
GAUTSCHI_HOST="${GAUTSCHI_HOST:-gautschi.rcac.purdue.edu}"
GAUTSCHI_ASM="${GAUTSCHI_ASM:-/scratch/gautschi/blackan/GROUSE/grouse_asm}"
REF_TAG="${REF_TAG:-ptarmigan}"                      # must match the assembly pipeline
OUTGROUP_ACC="${OUTGROUP_ACC:-GCF_023343835.1}"      # rock ptarmigan, Lagopus muta

# ---- destination (Negishi) ----
PROJECT_DIR="${PROJECT_DIR:-${CLUSTER_SCRATCH}/GROUSE/grouse_thex}"
IN_DIR="${PROJECT_DIR}/inputs/final_${REF_TAG}"     # copied assemblies + Liftoff GFF3
OUTGROUP_DIR="${PROJECT_DIR}/inputs/outgroup"        # copied ptarmigan FASTA + GFF
WORK="${PROJECT_DIR}/work"
ENV_PREFIX="${ENV_PREFIX:-${PROJECT_DIR}/envs/thex_env}"

# ---- SLURM (Negishi: 128 cores, 256 GB per cpu node; ~2 GB per core) ----
ACCOUNT="${ACCOUNT:-fnrdewoody}"
PARTITION="${PARTITION:-cpu}"
QOS="${QOS:-normal}"                                 # standby = 4 h limit, too short for prep/iqtree

# ---- alignment ----
REF_ID="${REF_ID:-}"                 # e.g. GRPC_F5596_hap1; empty = automatic (pick-ref)
OUTGROUP_NAME="${OUTGROUP_NAME:-Lmuta}"
PRESET_INGROUP="${PRESET_INGROUP:-asm5}"     # Tympanuchus vs Tympanuchus, < 1% divergence
PRESET_OUTGROUP="${PRESET_OUTGROUP:-asm20}"  # Lagopus vs Tympanuchus, a few % divergence
MIN_MAPQ="${MIN_MAPQ:-5}"
MIN_ALN="${MIN_ALN:-5000}"                   # minimum alignment block length (PAF col 11)
MIN_ALN_OUTGROUP="${MIN_ALN_OUTGROUP:-1000}"
EXCLUDE_CHROMS="${EXCLUDE_CHROMS:-chr_W,chr_MT}"
MIN_CHROM_COV="${MIN_CHROM_COV:-0.5}"        # drop a haplotype from a chromosome below this
PROJECT_CPUS="${PROJECT_CPUS:-32}"

# ---- THExBuilder ----
WINDOW="${WINDOW:-100kb}"
TRIMAL_GT="${TRIMAL_GT:-0.9}"               # trimal -gt: keep columns with <= 10% gaps
TRIMAL_MINLEN="${TRIMAL_MINLEN:-1kb}"        # drop windows shorter than this after trimal
PW_SUBWIN="${PW_SUBWIN:-100bp}"
PW_STEP="${PW_STEP:-10bp}"
PW_MAX_PDIST="${PW_MAX_PDIST:-0.15}"         # proportion; rescaled if the count bug is present
PW_ZSCORE="${PW_ZSCORE:-2}"                  # thexb's argparse default (0.15) is a typo for 2
PW_COV="${PW_COV:-0.9}"                      # every taxon needs >= 90% called bases, or the window is dropped
PW_FILTER="${PW_FILTER:-true}"               # false = skip it (trimal output goes straight to IQ-TREE)
PW_CHUNK="${PW_CHUNK:-40}"                   # windows per parallel thexb --pw_filter call
PREP_CPUS="${PREP_CPUS:-128}"

# ---- IQ-TREE ----
IQT_MODEL="${IQT_MODEL:-GTR+F+G4}"           # THEx's default GTR*H4 is several times slower
IQT_BOOT="${IQT_BOOT:-1000}"                 # UFBoot; THEx reads the .contree
IQT_THREADS="${IQT_THREADS:-2}"              # threads per IQ-TREE run
IQ_ARRAY="${IQ_ARRAY:-4}"                    # nodes (array tasks) sharing the windows
IQ_CPUS="${IQ_CPUS:-128}"                    # cores per array task
IQ_SEED="${IQ_SEED:-12345}"

SPECIES="${SPECIES:-STGR,LEPC,GRPC}"

# ---- derived paths ----
PAF_DIR="${WORK}/paf"
PROJ_DIR="${WORK}/projected"
MSA_DIR="${WORK}/msa"
THEXB_DIR="${WORK}/thexb"          # thexb working dir: windowed_fastas/, trimal_filtered_windows/, ...
PW_CHUNK_DIR="${WORK}/pw_chunks"
IQ_DIR="${WORK}/iqtree"
OUT_DIR="${WORK}/thex_inputs"
IDS_FILE="${WORK}/haplotype_ids.txt"      # 46 haplotypes, manifest order
NAMES_FILE="${WORK}/alignment_names.txt"  # 46 haplotypes + outgroup
CHROMS_FILE="${WORK}/chromosomes.txt"
REF_FILE="${WORK}/ref_id.txt"
MODE="${1:-}"

die()  { echo "ERROR: $*" >&2; exit 1; }
log()  { echo "[$(date '+%F %T')] $*"; }
sb()   { sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" --export=ALL "$@"; }
use_env() {
    [[ -x "${ENV_PREFIX}/bin/thexb" ]] || die "THEx env not found at ${ENV_PREFIX}; run: bash $0 setup"
    export PATH="${ENV_PREFIX}/bin:${PATH}"
    IQTREE_BIN="$(command -v iqtree2 || command -v iqtree3 || command -v iqtree || true)"
    [[ -n "$IQTREE_BIN" ]] || die "no iqtree2/iqtree3/iqtree in ${ENV_PREFIX}"
}

# Haplotype IDs from the manifest: <SPECIES>_<SAMPLE>_hap{1,2}, the assembly file prefix.
write_ids() {
    mkdir -p "$WORK"
    [[ -s "$MANIFEST" ]] || die "manifest not found: $MANIFEST (run from the repo root)"
    awk -F'\t' 'NR>1 && $0 !~ /^#/ && NF>=2 {
        gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)
        print $2"_"$1"_hap1"; print $2"_"$1"_hap2" }' "$MANIFEST" > "${IDS_FILE}.tmp"
    mv "${IDS_FILE}.tmp" "$IDS_FILE"
    { cat "$IDS_FILE"; echo "$OUTGROUP_NAME"; } > "$NAMES_FILE"
}

outgroup_fasta() { ls "${OUTGROUP_DIR}/${OUTGROUP_ACC}"_*_genomic.fna 2>/dev/null | head -1; }

resolve_ref() {
    if [[ -z "$REF_ID" ]]; then
        [[ -s "$REF_FILE" ]] || die "no reference chosen yet: run 'bash $0 pick-ref' (or set REF_ID)"
        REF_ID="$(cat "$REF_FILE")"
    fi
    REF_FA="${IN_DIR}/${REF_ID}.pseudo_chr.fasta"
    REF_GFF="${IN_DIR}/${REF_ID}.liftoff.gff3"
    [[ -s "$REF_FA" ]] || die "reference FASTA missing: $REF_FA"
}

# =============================================================================
# transfer: pull assemblies, Liftoff GFF3s and the ptarmigan reference from Gautschi
# =============================================================================
stage_transfer() {
    write_ids
    mkdir -p "$IN_DIR" "$OUTGROUP_DIR"
    local list="${WORK}/transfer_files.txt"
    sed 's/$/.pseudo_chr.fasta/' "$IDS_FILE" >  "$list"
    sed 's/$/.pseudo_chr.fasta.fai/' "$IDS_FILE" >> "$list"
    sed 's/$/.liftoff.gff3/' "$IDS_FILE" >> "$list"

    local need=0 f
    while read -r f; do [[ -s "${IN_DIR}/$f" ]] || need=1; done < "$list"
    [[ -n "$(outgroup_fasta)" ]] || need=1
    if [[ $need -eq 1 || "${FORCE_SYNC:-false}" == true ]]; then
        log "rsync ${GAUTSCHI_HOST}:${GAUTSCHI_ASM}/final_${REF_TAG}/ -> ${IN_DIR}/  ($(wc -l < "$list") files)"
        rsync -av --partial --files-from="$list" \
            "${GAUTSCHI_HOST}:${GAUTSCHI_ASM}/final_${REF_TAG}/" "${IN_DIR}/"
        log "rsync ptarmigan reference (${OUTGROUP_ACC}) -> ${OUTGROUP_DIR}/"
        rsync -av --partial --include="${OUTGROUP_ACC}_*_genomic.fna*" \
            --include="${OUTGROUP_ACC}_*_genomic.gff" --exclude='*' \
            "${GAUTSCHI_HOST}:${GAUTSCHI_ASM}/ref/" "${OUTGROUP_DIR}/"
    else
        log "all files already present in ${IN_DIR} (FORCE_SYNC=true to re-sync)"
    fi

    # verify: every file present, every FASTA indexed and agreeing with its index
    module load biocontainers samtools >/dev/null 2>&1 || true
    local missing=0
    while read -r f; do [[ -s "${IN_DIR}/$f" ]] || { echo "  MISSING ${IN_DIR}/$f"; missing=$((missing+1)); }; done < "$list"
    [[ $missing -eq 0 ]] || die "$missing files missing after transfer"
    local OG; OG="$(outgroup_fasta)"
    [[ -n "$OG" ]] || die "ptarmigan FASTA ${OUTGROUP_ACC}_*_genomic.fna not found in ${GAUTSCHI_ASM}/ref/"
    for f in "${IN_DIR}"/*.pseudo_chr.fasta "$OG"; do
        if [[ ! -s "${f}.fai" || "${f}.fai" -ot "$f" ]]; then
            command -v samtools >/dev/null || die "samtools needed to index $f"
            samtools faidx "$f"
        fi
    done
    log "transfer OK: $(wc -l < "$IDS_FILE") haplotypes + $(basename "$OG")"
    awk '{s+=$2} END {printf "  outgroup total length: %.2f Gb\n", s/1e9}' "${OG}.fai"
}

# =============================================================================
# setup: THEx conda environment on scratch (login node; needs network)
# =============================================================================
stage_setup() {
    module load conda >/dev/null 2>&1 || module load anaconda >/dev/null 2>&1 || true
    command -v conda >/dev/null || die "no conda: try 'module avail conda anaconda'"
    export CONDA_PKGS_DIRS="${PROJECT_DIR}/envs/pkgs"   # keep the package cache off the home quota
    mkdir -p "$CONDA_PKGS_DIRS"
    if [[ ! -x "${ENV_PREFIX}/bin/thexb" ]]; then
        local solver; solver="$(command -v mamba || command -v conda)"
        log "creating ${ENV_PREFIX} with $(basename "$solver")"
        "$solver" create -y -p "$ENV_PREFIX" -c ajharris_2374 -c bioconda -c conda-forge \
            thex iqtree trimal
    fi
    use_env
    log "thexb: $(thexb --version 2>&1 | tail -1)"
    log "iqtree: $("$IQTREE_BIN" --version 2>&1 | grep -m1 -i version)"
    log "trimal: $(trimal --version 2>&1 | head -1)"
    python -c "import ete3, pandas, openpyxl, pyfaidx" || die "python deps missing from env"
}

# =============================================================================
# pick-ref: choose the coordinate haplotype; write the chromosome list
# =============================================================================
stage_pick_ref() {
    write_ids
    if [[ -z "$REF_ID" ]]; then
        REF_ID="$(python3 "$PYHELP" pick-ref --ids "$IDS_FILE" --fai-dir "$IN_DIR" \
                  --exclude "$EXCLUDE_CHROMS" --out-table "${WORK}/ref_choice.tsv")"
        log "automatic reference choice (top of ${WORK}/ref_choice.tsv):"
        head -6 "${WORK}/ref_choice.tsv"
    fi
    grep -qx "$REF_ID" "$IDS_FILE" || die "REF_ID=$REF_ID is not one of the manifest haplotypes"
    if [[ -s "$REF_FILE" && "$(cat "$REF_FILE")" != "$REF_ID" && -d "$PROJ_DIR" ]]; then
        die "work/ already holds projections onto $(cat "$REF_FILE"); move ${WORK} aside before switching to $REF_ID"
    fi
    echo "$REF_ID" > "$REF_FILE"
    resolve_ref
    awk -v ex=",${EXCLUDE_CHROMS}," '$1 ~ /^chr_/ && index(ex, ","$1",")==0 {print $1}' \
        "${REF_FA}.fai" > "$CHROMS_FILE"
    log "reference haplotype: $REF_ID  ($(wc -l < "$CHROMS_FILE") chromosomes: $(paste -sd' ' "$CHROMS_FILE"))"
    grep -qx chr_Z "$CHROMS_FILE" || echo "WARNING: $REF_ID has no chr_Z, so the Z will have no windows"
}

# =============================================================================
# project (array, one task per haplotype + outgroup): align, then project
# =============================================================================
stage_project() {
    resolve_ref
    local i="${SLURM_ARRAY_TASK_ID:?run as an array job}"
    local NAME; NAME="$(sed -n "$((i+1))p" "$NAMES_FILE")"
    [[ -n "$NAME" ]] || die "no alignment name for task $i"
    local done_flag="${PROJ_DIR}/done/${NAME}.ok"
    [[ -e "$done_flag" ]] && { log "$NAME already projected"; return 0; }
    mkdir -p "$PAF_DIR" "${PROJ_DIR}/done"
    use_env    # python with nothing but the standard library is fine too; this one is known

    local Q PRESET MINALN
    if [[ "$NAME" == "$OUTGROUP_NAME" ]]; then
        Q="$(outgroup_fasta)"; PRESET="$PRESET_OUTGROUP"; MINALN="$MIN_ALN_OUTGROUP"
    else
        Q="${IN_DIR}/${NAME}.pseudo_chr.fasta"; PRESET="$PRESET_INGROUP"; MINALN="$MIN_ALN"
    fi
    [[ -s "$Q" ]] || die "query FASTA missing: $Q"

    if [[ "$NAME" == "$REF_ID" ]]; then
        python "$PYHELP" project --self --ref "$REF_FA" --name "$NAME" \
            --chroms "$CHROMS_FILE" --outdir "$PROJ_DIR"
    else
        local PAF="${PAF_DIR}/${NAME}.paf"
        if [[ ! -s "${PAF}.ok" ]]; then
            module load biocontainers minimap2 >/dev/null 2>&1 || true
            command -v minimap2 >/dev/null || die "minimap2 not found (module load biocontainers minimap2)"
            log "minimap2 -x $PRESET: $NAME -> $REF_ID"
            # -c --cs: base-level alignment with the short cs tag the projection reads
            minimap2 -c --cs -x "$PRESET" --secondary=no -t "${SLURM_CPUS_PER_TASK:-$PROJECT_CPUS}" \
                "$REF_FA" "$Q" > "${PAF}.part"
            mv "${PAF}.part" "$PAF"
            minimap2 --version > "${PAF}.ok"
        fi
        python "$PYHELP" project --ref "$REF_FA" --paf "$PAF" --name "$NAME" \
            --chroms "$CHROMS_FILE" --outdir "$PROJ_DIR" --min-mapq "$MIN_MAPQ" --min-aln "$MINALN"
    fi
    touch "$done_flag"
}

# =============================================================================
# prep (one node): MSAs, annotation track, thexb windowing + trimal + pairwise filter
# =============================================================================
stage_prep() {
    resolve_ref
    use_env
    local CPUS="${SLURM_CPUS_PER_TASK:-$PREP_CPUS}"
    local n_done n_need
    n_done=$(ls "${PROJ_DIR}/done/" 2>/dev/null | wc -l); n_need=$(wc -l < "$NAMES_FILE")
    [[ "$n_done" -eq "$n_need" ]] || die "only ${n_done}/${n_need} projections finished"

    # ---- 1. one alignment per chromosome ----
    if [[ ! -e "${MSA_DIR}/.done" ]]; then
        log "building per-chromosome alignments"
        python "$PYHELP" msa --projdir "$PROJ_DIR" --chroms "$CHROMS_FILE" --names "$NAMES_FILE" \
            --outdir "$MSA_DIR" --min-cov "$MIN_CHROM_COV" --always "${REF_ID},${OUTGROUP_NAME}" \
            --log "${WORK}/samples_per_chromosome.tsv"
        touch "${MSA_DIR}/.done"
    fi
    awk -F'\t' '$4=="no"' "${WORK}/samples_per_chromosome.tsv" | sed 's/^/  excluded: /' | head -40

    # ---- 2. THEx annotation track + chromosome lengths ----
    mkdir -p "$OUT_DIR"
    if [[ -s "$REF_GFF" ]]; then
        python "$PYHELP" annot --gff "$REF_GFF" --fai "${REF_FA}.fai" --chroms "$CHROMS_FILE" \
            --out-gff "${OUT_DIR}/reference_genes.gff" --out-bed "${OUT_DIR}/chromosome_lengths.bed"
    else
        echo "WARNING: no Liftoff GFF3 for $REF_ID; no gene track"
    fi

    # thexb writes into its -o dir and needs -i to exist; its stages read fixed subdirs.
    mkdir -p "$THEXB_DIR"
    cd "$THEXB_DIR"

    # ---- 3. windows ----
    if [[ ! -e windowed_fastas/.done ]]; then
        rm -rf windowed_fastas
        log "thexb --minifastas ($WINDOW)"
        thexb --minifastas -i "$MSA_DIR" -o "$THEXB_DIR" -w "$WINDOW" --cpu "$CPUS"
        local nchr; nchr=$(find windowed_fastas -mindepth 1 -maxdepth 1 -type d | wc -l)
        [[ "$nchr" -eq "$(wc -l < "$CHROMS_FILE")" ]] || die "minifastas wrote $nchr of $(wc -l < "$CHROMS_FILE") chromosomes; see ${THEXB_DIR}/logs/minifastas.log"
        touch windowed_fastas/.done
    fi
    log "$(find windowed_fastas -name '*.fasta' | wc -l) windows"

    # ---- 4. trimal ----
    if [[ ! -e trimal_filtered_windows/.done ]]; then
        rm -rf trimal_filtered_windows
        log "thexb --trimal (gap threshold $TRIMAL_GT, min length $TRIMAL_MINLEN)"
        thexb --trimal -i "$MSA_DIR" -o "$THEXB_DIR" --trimal_gap_threshold "$TRIMAL_GT" \
            --trimal_min_seq_len "$TRIMAL_MINLEN" --cpu "$CPUS"
        touch trimal_filtered_windows/.done
    fi
    log "$(find trimal_filtered_windows -name '*.fasta' ! -name '*-DROPPED*' | wc -l) windows pass trimal"

    # ---- 5. pairwise distance + coverage filter, chunked so it uses every core ----
    # (thexb parallelises this stage only across chromosomes, and it is pure Python.)
    if [[ "$PW_FILTER" != true ]]; then
        log "PW_FILTER=false: trimal output goes to IQ-TREE unfiltered"
        rm -rf pairwise_filtered_windows
        mkdir -p pairwise_filtered_windows
        cp -rs "${THEXB_DIR}/trimal_filtered_windows/"* pairwise_filtered_windows/
        find pairwise_filtered_windows -name '*.fai' -delete
        touch pairwise_filtered_windows/.done
    elif [[ ! -e pairwise_filtered_windows/.done ]]; then
        local maxp="$PW_MAX_PDIST" subwin_bp
        subwin_bp=$(python -c "
s='${PW_SUBWIN}'.lower(); m={'bp':1,'kb':1000,'mb':1000000}
u=[k for k in ('kb','mb','bp') if s.endswith(k)][0]; print(int(float(s[:-2])*m[u]))")
        if python - <<'EOF'
import inspect, sys, thexb.STAGE_pairwise_filter as m
sys.exit(0 if "countmismatch + countmissing / denominator" in inspect.getsource(m) else 1)
EOF
        then
            maxp=$(python -c "print(${PW_MAX_PDIST} * ${subwin_bp})")
            log "installed thexb has the mismatch-count bug: --pw_max_pdist ${maxp} (= ${PW_MAX_PDIST} x ${subwin_bp} bp)"
        fi
        rm -rf "$PW_CHUNK_DIR" pairwise_filtered_windows
        mkdir -p "$PW_CHUNK_DIR" pairwise_filtered_windows
        # chunk dirs: <chunk>/trimal_filtered_windows/<chrom>/ holding symlinks
        local chromdir chrom cdir k=0
        for chromdir in trimal_filtered_windows/*/; do
            chrom="$(basename "$chromdir")"
            mkdir -p "pairwise_filtered_windows/${chrom}"
            find "$chromdir" -maxdepth 1 -name '*.fasta' | sort | split -l "$PW_CHUNK" -d -a 5 - "${PW_CHUNK_DIR}/${chrom}.list."
        done
        for l in "${PW_CHUNK_DIR}"/*.list.*; do
            chrom="$(basename "$l" | sed 's/\.list\.[0-9]*$//')"
            cdir="${PW_CHUNK_DIR}/c$(printf '%05d' "$k")"
            mkdir -p "${cdir}/trimal_filtered_windows/${chrom}"
            while read -r f; do ln -s "${THEXB_DIR}/${f}" "${cdir}/trimal_filtered_windows/${chrom}/"; done < "$l"
            echo "$cdir" >> "${PW_CHUNK_DIR}/chunks.txt"
            k=$((k+1))
        done
        log "thexb --pw_filter on $k chunks, ${CPUS} at a time (sub-window $PW_SUBWIN, step $PW_STEP, z $PW_ZSCORE, coverage $PW_COV)"
        export REF_ID OUTGROUP_NAME PW_SUBWIN PW_STEP PW_ZSCORE PW_COV maxp
        xargs -P "$CPUS" -I{} bash -c '
            thexb --pw_filter -i "$1" -o "$1" -r "$REF_ID" --pw_exclude "${REF_ID},${OUTGROUP_NAME}" \
                --pw_subwindow_size "$PW_SUBWIN" --pw_step "$PW_STEP" --pw_max_pdist "$maxp" \
                --pw_zscore_cutoff "$PW_ZSCORE" --pw_seq_coverage "$PW_COV" --pw_missing_char N \
                --cpu 1 > "$1/pw.log" 2>&1 || echo "FAILED $1"' _ {} < "${PW_CHUNK_DIR}/chunks.txt" \
            | tee "${PW_CHUNK_DIR}/failed.txt"
        [[ ! -s "${PW_CHUNK_DIR}/failed.txt" ]] || die "$(wc -l < "${PW_CHUNK_DIR}/failed.txt") pw_filter chunks failed; see their pw.log"
        # gather
        while read -r cdir; do
            for d in "${cdir}/pairwise_filtered_windows/"*/; do
                mv "$d"*.fasta "pairwise_filtered_windows/$(basename "$d")/" 2>/dev/null || true
            done
        done < "${PW_CHUNK_DIR}/chunks.txt"
        local n_in n_out
        n_in=$(find trimal_filtered_windows -name '*.fasta' | wc -l)
        n_out=$(find pairwise_filtered_windows -name '*.fasta' | wc -l)
        [[ "$n_in" -eq "$n_out" ]] || die "pw_filter returned $n_out of $n_in windows"
        cat "${PW_CHUNK_DIR}"/c*/logs/pairwise_filter.log > logs/pairwise_filter_chunks.log 2>/dev/null || true
        rm -rf "$PW_CHUNK_DIR"
        touch pairwise_filtered_windows/.done
    fi
    find pairwise_filtered_windows -name '*.fasta' ! -name '*-DROPPED*' | sort > "${WORK}/iqtree_windows.txt"
    log "$(wc -l < "${WORK}/iqtree_windows.txt") windows go to IQ-TREE"
}

# =============================================================================
# iqtree (array): one ML tree + UFBoot per window; task t takes every IQ_ARRAY-th window
# =============================================================================
run_one_window() {   # $1 = window FASTA path relative to THEXB_DIR
    local f="$1" chrom stem pre
    chrom="$(basename "$(dirname "$f")")"; stem="$(basename "$f" .fasta)"
    pre="${IQ_DIR}/${chrom}/${stem}"
    [[ -s "${pre}.contree" ]] && return 0
    mkdir -p "${IQ_DIR}/${chrom}"
    if "$IQTREE_BIN" -s "${THEXB_DIR}/${f}" -m "$IQT_MODEL" -B "$IQT_BOOT" -T "$IQT_THREADS" \
        --prefix "$pre" --seed "$IQ_SEED" --quiet > /dev/null 2>&1 && [[ -s "${pre}.contree" ]]; then
        rm -f "${pre}".{ckp.gz,bionj,mldist,uniqueseq.phy,model.gz,splits.nex}
    else
        echo "FAILED ${f} (see ${pre}.log)"
    fi
}
stage_iqtree() {
    use_env
    local t="${SLURM_ARRAY_TASK_ID:-0}" n="${SLURM_ARRAY_TASK_COUNT:-$IQ_ARRAY}"
    local CPUS="${SLURM_CPUS_PER_TASK:-$IQ_CPUS}" P
    P=$(( CPUS / IQT_THREADS )); [[ $P -ge 1 ]] || P=1
    [[ -s "${WORK}/iqtree_windows.txt" ]] || die "no window list; run the prep stage first"
    mkdir -p "$IQ_DIR"
    local mine="${IQ_DIR}/task_${t}.windows"
    awk -v t="$t" -v n="$n" '(NR-1) % n == t' "${WORK}/iqtree_windows.txt" > "$mine"
    log "task $t/$n: $(wc -l < "$mine") windows, $P IQ-TREE runs x $IQT_THREADS threads ($IQT_MODEL, UFBoot $IQT_BOOT)"
    export -f run_one_window
    export IQ_DIR THEXB_DIR IQTREE_BIN IQT_MODEL IQT_BOOT IQT_THREADS IQ_SEED
    xargs -P "$P" -I{} bash -c 'run_one_window "$1"' _ {} < "$mine" | tee "${IQ_DIR}/task_${t}.failed"
    local nf; nf=$(grep -c FAILED "${IQ_DIR}/task_${t}.failed" || true)
    log "task $t done; $nf windows failed"
    [[ "$nf" -eq 0 ]] || echo "Failed windows become NoTree in THEx; resubmit this stage to retry them."
}

# =============================================================================
# summarize: Tree Viewer input
# =============================================================================
stage_summarize() {
    resolve_ref
    use_env
    mkdir -p "$OUT_DIR"
    python "$PYHELP" summarize \
        --windowed-dir "${THEXB_DIR}/windowed_fastas" \
        --trimal-dir "${THEXB_DIR}/trimal_filtered_windows" \
        --pw-dir "${THEXB_DIR}/pairwise_filtered_windows" \
        --iqtree-dir "$IQ_DIR" --names "$NAMES_FILE" --chroms "$CHROMS_FILE" \
        --gff "${OUT_DIR}/reference_genes.gff" --outgroup "$OUTGROUP_NAME" --species "$SPECIES" \
        --out-xlsx "${OUT_DIR}/grouse_TreeViewer_input.xlsx" --out-prefix "${OUT_DIR}/grouse_thex"
    cp "${WORK}/samples_per_chromosome.tsv" "${WORK}/ref_choice.tsv" "$OUT_DIR"/ 2>/dev/null || true
    cat > "${OUT_DIR}/README.txt" <<EOF
THEx Tree Viewer inputs, prairie grouse, $(date +%F)
Coordinates: ${REF_ID} (chromosome names as in THEx: '_' -> '-')
Tips: 46 haplotypes (<SPECIES>_<SAMPLE>_hap<N>) + ${OUTGROUP_NAME} (rock ptarmigan, outgroup; trees rooted on it)
Windows: ${WINDOW}; trimal -gt ${TRIMAL_GT}; pairwise filter ${PW_FILTER} (coverage ${PW_COV}); IQ-TREE ${IQT_MODEL}, UFBoot ${IQT_BOOT} (.contree)

grouse_TreeViewer_input.xlsx   Tree Viewer file. TopologyID = species-level topology (see 05_thex_phylogeny.sh header)
chromosome_lengths.bed         load as the chromosome-length file
reference_genes.gff            Liftoff genes on ${REF_ID}, load as the GFF/GTF annotation
grouse_thex.topology_summary.tsv, .topology_by_chromosome.tsv, .window_status.tsv
EOF
    log "Tree Viewer inputs in ${OUT_DIR}/"
}

# =============================================================================
# submit: the whole chain on Negishi
# =============================================================================
stage_submit() {
    mkdir -p "${REPO_DIR}/logs"
    [[ -s "$PYHELP" ]] || die "missing $PYHELP"
    [[ -x "${ENV_PREFIX}/bin/thexb" ]] || die "run 'bash $0 setup' first"
    stage_pick_ref
    local N; N=$(wc -l < "$NAMES_FILE")
    local J1 J2 J3 J4
    J1=$(sb --array=0-$((N-1)) -c "$PROJECT_CPUS" --mem=$(( PROJECT_CPUS * 2 ))G -t 12:00:00 \
            -J thex_project -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" project)
    J2=$(sb --dependency=afterok:"$J1" -N 1 -c "$PREP_CPUS" --mem=240G -t 2-00:00:00 \
            -J thex_prep -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" prep)
    J3=$(sb --dependency=afterok:"$J2" --array=0-$((IQ_ARRAY-1)) -N 1 -c "$IQ_CPUS" --mem=240G \
            -t 4-00:00:00 -J thex_iqtree -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" iqtree)
    J4=$(sb --dependency=afterany:"$J3" -c 8 --mem=32G -t 12:00:00 \
            -J thex_summary -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" summarize)
    log "submitted: project $J1 ($N tasks) -> prep $J2 -> iqtree $J3 ($IQ_ARRAY nodes) -> summarize $J4"
}

stage_status() {
    resolve_ref
    echo "reference:   $REF_ID"
    echo "projected:   $(ls "${PROJ_DIR}/done" 2>/dev/null | wc -l) / $(wc -l < "$NAMES_FILE")"
    [[ -d "${THEXB_DIR}/windowed_fastas" ]] && echo "windows:     $(find "${THEXB_DIR}/windowed_fastas" -name '*.fasta' | wc -l)"
    [[ -d "${THEXB_DIR}/trimal_filtered_windows" ]] && echo "pass trimal: $(find "${THEXB_DIR}/trimal_filtered_windows" -name '*.fasta' ! -name '*-DROPPED*' | wc -l)"
    [[ -s "${WORK}/iqtree_windows.txt" ]] && echo "to IQ-TREE:  $(wc -l < "${WORK}/iqtree_windows.txt")"
    [[ -d "$IQ_DIR" ]] && echo "trees done:  $(find "$IQ_DIR" -name '*.contree' | wc -l)"
    return 0
}

case "$MODE" in
    transfer)  stage_transfer ;;
    setup)     stage_setup ;;
    pick-ref)  stage_pick_ref ;;
    submit)    stage_submit ;;
    project)   stage_project ;;
    prep)      stage_prep ;;
    iqtree)    stage_iqtree ;;
    summarize) stage_summarize ;;
    status)    stage_status ;;
    *) sed -n '2,91p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
