#!/bin/bash
# =============================================================================
# PANGENOME INPUT PREP 02b — collate step 02 into one table and flag outliers
#
#   bash 02b_collate_qc.sh          (login node, seconds)
#
# OUTPUT ${PROJECT_DIR}/pangenome_qc/qc_summary.tsv, one row per haplotype:
#   QV, k-mer completeness (Merqury), median HiFi depth, % of the haplotype in
#   low/high-depth windows, number of joins cut in step 01, and REVIEW reasons.
#
# The REVIEW thresholds are triage cut-offs for deciding where to look first,
# not pass/fail criteria; change them below to suit what the panel looks like.
# =============================================================================
set -euo pipefail

show_tsv() { if command -v column >/dev/null 2>&1; then column -t -s$'\t' "$1"; else cat "$1"; fi; }

PG_WORK="${PG_WORK:-${CLUSTER_SCRATCH}/GROUSE/grouse_asm}"
QC_DIR="${PG_WORK}/pangenome_qc"
IN_SUMMARY="${PG_WORK}/pangenome_inputs/summary.tsv"
OUT="${QC_DIR}/qc_summary.tsv"

MIN_QV=40           # below -> review
COMP_DROP=5         # flag a haplotype whose k-mer completeness is this many
                    # points below the panel median. Per-haplotype completeness
                    # depends on heterozygosity, so it is judged against the
                    # other haplotypes rather than against a fixed number.
MAX_LOW_PCT=5       # % of haplotype in < 0.5x-median windows
MAX_HIGH_PCT=5      # % of haplotype in > 1.5x-median windows

[[ -d "$QC_DIR" ]] || { echo "ERROR: $QC_DIR not found — run step 02 first"; exit 1; }

TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT

# haplotype -> qv / completeness / depth fields, from every sample directory
for D in "$QC_DIR"/*/; do
    S=$(basename "$D")
    MQ="${D}merqury"
    [[ -s "${MQ}/${S}.qv" ]] && awk -F'\t' '$1 != "Both" {print $1"\tqv\t"$4}' "${MQ}/${S}.qv" >> "$TMP"
    [[ -s "${MQ}/${S}.completeness.stats" ]] && \
        awk -F'\t' '$1 != "Both" {print $1"\tcomp\t"$5}' "${MQ}/${S}.completeness.stats" >> "$TMP"
    [[ -s "${D}${S}.depth_summary.tsv" ]] && \
        awk -F'\t' 'NR > 1 {print $1"\tmed\t"$3; print $1"\tlow\t"$5; print $1"\thigh\t"$7}' \
            "${D}${S}.depth_summary.tsv" >> "$TMP"
done
[[ -s "$IN_SUMMARY" ]] && awk -F'\t' 'NR > 1 {print $1"\tcuts\t"$5}' "$IN_SUMMARY" >> "$TMP"

[[ -s "$TMP" ]] || { echo "ERROR: no step 02 results found under $QC_DIR"; exit 1; }

MED_COMP=$(awk -F'\t' '$2 == "comp" {print $3}' "$TMP" | sort -g \
           | awk '{v[NR] = $1} END {if (NR) print v[int((NR + 1) / 2)]; else print 0}')
MIN_COMPLETE=$(awk -v m="$MED_COMP" -v d="$COMP_DROP" 'BEGIN {print m - d}')
echo ">>> panel median completeness ${MED_COMP}% -> flag below ${MIN_COMPLETE}%"
echo ""

awk -F'\t' -v OFS='\t' -v qv="$MIN_QV" -v cp="$MIN_COMPLETE" \
    -v lo="$MAX_LOW_PCT" -v hi="$MAX_HIGH_PCT" '
    { h[$1] = 1; v[$1, $2] = $3 }
    END {
      print "haplotype", "QV", "completeness_pct", "median_depth", "low_depth_pct",
            "high_depth_pct", "joins_cut", "REVIEW"
      for (k in h) {
        r = ""
        if ((k, "qv") in v && v[k, "qv"] + 0 < qv)          r = r "low_QV;"
        if ((k, "comp") in v && v[k, "comp"] + 0 < cp)      r = r "low_completeness;"
        if ((k, "low") in v && v[k, "low"] + 0 > lo)        r = r "low_depth_windows;"
        if ((k, "high") in v && v[k, "high"] + 0 > hi)      r = r "high_depth_windows;"
        if (!((k, "qv") in v) || !((k, "med") in v))        r = r "qc_incomplete;"
        print k, ((k, "qv") in v ? v[k, "qv"] : "NA"), ((k, "comp") in v ? v[k, "comp"] : "NA"),
              ((k, "med") in v ? v[k, "med"] : "NA"), ((k, "low") in v ? v[k, "low"] : "NA"),
              ((k, "high") in v ? v[k, "high"] : "NA"), ((k, "cuts") in v ? v[k, "cuts"] : "NA"),
              (r == "" ? "-" : r)
      }
    }' "$TMP" | (read -r HEADER; echo "$HEADER"; sort) > "$OUT"

show_tsv "$OUT"
echo ""
echo ">>> written: $OUT"
echo ">>> $(awk -F'\t' 'NR > 1 && $8 != "-"' "$OUT" | wc -l) haplotype(s) flagged for review"
echo "    Sex-chromosome-depleted female haplotypes are EXPECTED to show lower completeness."
