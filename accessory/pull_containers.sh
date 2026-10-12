#!/bin/bash
# =============================================================================
# accessory/pull_containers.sh -- fetch the Apptainer images steps 03-04 use
# for tools with no Gautschi module (pggb + odgi, vg, RepeatMasker).
# Run ONCE on a LOGIN node (compute nodes may have no internet):
#   bash accessory/pull_containers.sh
# Images go to ${CLUSTER_SCRATCH}/GROUSE/grouse_asm/containers/, the default
# location read by 03_pangenome_build.sh (PGGB_SIF, TETOOLS_SIF) and
# 04_shortread_introgression.sh (VG_SIF). Versions actually pulled are written
# to containers/VERSIONS.txt for the report's methods.
# =============================================================================
set -euo pipefail
DEST="${CLUSTER_SCRATCH}/GROUSE/grouse_asm/containers"
export APPTAINER_CACHEDIR="${CLUSTER_SCRATCH}/.apptainer_cache"   # keep $HOME quota free
mkdir -p "$DEST" "$APPTAINER_CACHEDIR"
AP="$(command -v apptainer || command -v singularity)" || { echo "ERROR: no apptainer/singularity"; exit 1; }
export APPTAINERENV_LD_PRELOAD="" SINGULARITYENV_LD_PRELOAD=""

# Edit tags here to pin other versions.
declare -A IMAGES=(
  [pggb]="docker://ghcr.io/pangenome/pggb:latest"     # pggb + wfmash/seqwish/smoothxg/odgi
  [vg]="docker://quay.io/vgteam/vg:v1.59.0"           # vg giraffe/call/autoindex
  [tetools]="docker://dfam/tetools:1.88"              # RepeatModeler2 + RepeatMasker (only used if
                                                      # the repeatmodeler module lacks RepeatMasker)
)
for name in pggb vg tetools; do
  sif="${DEST}/${name}.sif"
  if [[ -s "$sif" ]]; then echo "== ${name}: ${sif} exists, skipping"; continue; fi
  echo "== ${name}: pulling ${IMAGES[$name]}"
  "$AP" pull "$sif" "${IMAGES[$name]}"
done

{
  echo "# pulled $(date)"
  for name in pggb vg tetools; do echo "${name}	${IMAGES[$name]}"; done
  echo "pggb version: $("$AP" exec "${DEST}/pggb.sif" pggb --version 2>&1 | head -1)"
  echo "odgi version: $("$AP" exec "${DEST}/pggb.sif" odgi version 2>&1 | head -1)"
  echo "vg version:   $("$AP" exec "${DEST}/vg.sif" vg version 2>&1 | head -1)"
  echo "RepeatMasker: $("$AP" exec "${DEST}/tetools.sif" RepeatMasker -v 2>&1 | head -1)"
} | tee "${DEST}/VERSIONS.txt"

echo
echo "Does the repeatmodeler module already provide RepeatMasker? (if yes, tetools.sif is unused)"
( module load biocontainers repeatmodeler/2.0.3 2>/dev/null && command -v RepeatMasker ) || echo "  no -> tetools.sif will be used"
