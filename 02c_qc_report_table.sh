#!/bin/bash
# =============================================================================
# 02c_qc_report_table.sh -- collect the per-sample pangenome-input QC
# (Merqury QV + k-mer completeness, HiFi depth) into one table and flag
# outliers. Runs in seconds on a login node; no modules needed.
#
#   bash 02c_qc_report_table.sh               # default QC_DIR below (after step 02)
#   bash 02c_qc_report_table.sh /path/to/pangenome_qc
#
# Complements 02b (per-haplotype triage with REVIEW reasons): this adds the
# combined hap1 + hap2 Merqury row per bird (the headline completeness/QV for
# a diploid assembly), panel ranges, and gathers every spectra-cn plot into
# one folder. Its table is the input to 02d_plot_qc_report.R (report Fig S8).
#
# Reads, per sample directory <QC_DIR>/<SAMPLE>/:
#   merqury/<SAMPLE>.qv                  asm  asm_only_kmers  total_kmers  QV  error
#   merqury/<SAMPLE>.completeness.stats  asm  all  found  total  completeness(%)
#   <SAMPLE>.depth_summary.tsv           haplotype windows median_depth low_windows
#                                        low_bp_pct high_windows high_bp_pct
# Writes to <QC_DIR>/summary/:
#   pangenome_qc_summary.tsv   one row per haplotype + one "both" row per sample
#   pangenome_qc_flags.txt     anything outside the thresholds below
#   spectra_cn/                copies of every spectra-cn PNG, for quick review
# =============================================================================
set -euo pipefail
QC_DIR="${1:-/scratch/gautschi/blackan/GROUSE/grouse_asm/pangenome_qc}"
OUT="${QC_DIR}/summary"
mkdir -p "${OUT}/spectra_cn"

# Review thresholds (flags only; nothing is excluded)
MIN_QV=50            # per haplotype; Q50 = 1 error per 100 kb
MIN_COMP_BOTH=98     # % of read solid k-mers found in hap1 + hap2 together
MIN_COMP_HAP=85      # % per haplotype (expected ~88-95 in a diploid; lower
                     # in the W-bearing haplotype of a female, which lacks the Z)
MAX_LOW_BP=5         # % of a haplotype in low-depth windows (possible false dup / errors)
MAX_HIGH_BP=1        # % in high-depth windows (possible collapsed repeats)

TAB="${OUT}/pangenome_qc_summary.tsv"
printf "sample\tspecies\tassembly\tlevel\tQV\terror_rate\tasm_only_kmers\tcompleteness_pct\tfound_kmers\tread_kmers\tmedian_depth\tlow_bp_pct\thigh_bp_pct\n" > "$TAB"
: > "${OUT}/pangenome_qc_flags.txt"

n_ok=0; missing=()
for d in "${QC_DIR}"/*/; do
    s=$(basename "$d"); [[ "$s" == summary ]] && continue
    qv="${d}merqury/${s}.qv"; cs="${d}merqury/${s}.completeness.stats"; dp="${d}${s}.depth_summary.tsv"
    if [[ ! -s "$qv" || ! -s "$cs" ]]; then missing+=("$s"); continue; fi
    n_ok=$((n_ok+1))
    cp -f "${d}"merqury/*.spectra-cn.*png "${OUT}/spectra_cn/" 2>/dev/null || true
    # Merqury files are tab-delimited; assembly names can contain a space
    # (manifest "LEPC " -> "LEPC _F5540_hap1"), so never split on spaces.
    awk -F'\t' -v OFS='\t' -v s="$s" -v dp="$dp" '
        FNR==1 { f++ }
        NF==1 { n=split($0, t, / +/); for (i=1;i<=n;i++) $i=t[i] }   # space-delimited fallback
        f==1 { qv[$1]=$4; err[$1]=$5; only[$1]=$2; next }                 # .qv
        f==2 { a=($1=="both")?"Both":$1; comp[a]=$5; fnd[a]=$3; tot[a]=$4; next }  # completeness
        f==3 && FNR>1 { med[$1]=$3; low[$1]=$5; high[$1]=$7; next }       # depth
        END {
            for (a in qv) {
                lvl = (a=="Both") ? "both" : "haplotype"
                asm = (a=="Both") ? s "_both" : a
                sp = a; sub(/_F[0-9]+_hap[12].*$/, "", sp); gsub(/ /, "", sp)
                if (a=="Both") sp="."
                print s, sp, asm, lvl, qv[a], err[a], only[a], comp[a], fnd[a], tot[a],
                      (a in med ? med[a] : "NA"), (a in low ? low[a] : "NA"), (a in high ? high[a] : "NA")
            }
        }' "$qv" "$cs" $( [[ -s "$dp" ]] && echo "$dp" || echo /dev/null ) >> "$TAB"
done

# fill species on the "both" rows from that sample's haplotype rows
awk -F'\t' -v OFS='\t' 'NR==FNR { if (NR>1 && $4=="haplotype") sp[$1]=$2; next }
    FNR==1 || $2!="." { print; next } { $2=sp[$1]; print }' "$TAB" "$TAB" > "${TAB}.tmp" && mv "${TAB}.tmp" "$TAB"
# sort: species, sample, level
{ head -1 "$TAB"; tail -n +2 "$TAB" | sort -t$'\t' -k2,2 -k1,1 -k4,4r -k3,3; } > "${TAB}.tmp" && mv "${TAB}.tmp" "$TAB"

awk -F'\t' -v q="$MIN_QV" -v cb="$MIN_COMP_BOTH" -v ch="$MIN_COMP_HAP" -v lo="$MAX_LOW_BP" -v hi="$MAX_HIGH_BP" '
    NR==1 { next }
    $4=="haplotype" && $5+0 < q   { print "QV < " q ":\t" $3 "\tQV=" $5 }
    $4=="both"      && $8+0 < cb  { print "completeness (both) < " cb "%:\t" $1 "\t" $8 }
    $4=="haplotype" && $8+0 < ch  { print "completeness (hap) < " ch "%:\t" $3 "\t" $8 }
    $4=="haplotype" && $12!="NA" && $12+0 > lo { print "low-depth > " lo "% of bp:\t" $3 "\t" $12 }
    $4=="haplotype" && $13!="NA" && $13+0 > hi { print "high-depth > " hi "% of bp:\t" $3 "\t" $13 }
' "$TAB" > "${OUT}/pangenome_qc_flags.txt"

echo "Samples summarized: ${n_ok}"
if (( ${#missing[@]} )); then echo "MISSING Merqury output: ${missing[*]}"; fi
echo; echo "== Per-haplotype and combined metrics (${TAB}) =="
if command -v column >/dev/null; then column -t -s$'\t' "$TAB"; else cat "$TAB"; fi
echo; echo "== Range by level =="
awk -F'\t' 'NR>1 { k=$4; n[k]++; if(!(k in qmin)||$5<qmin[k])qmin[k]=$5; if($5>qmax[k])qmax[k]=$5;
                   if(!(k in cmin)||$8<cmin[k])cmin[k]=$8; if($8>cmax[k])cmax[k]=$8 }
    END { for (k in n) printf "  %-9s n=%-3d QV %.1f-%.1f   completeness %.2f-%.2f%%\n", k, n[k], qmin[k], qmax[k], cmin[k], cmax[k] }' "$TAB"
echo; echo "== Flags (${OUT}/pangenome_qc_flags.txt) =="
if [[ -s "${OUT}/pangenome_qc_flags.txt" ]]; then cat "${OUT}/pangenome_qc_flags.txt"; else echo "  none"; fi
echo; echo "spectra-cn plots copied to ${OUT}/spectra_cn/ ($(ls "${OUT}/spectra_cn" | wc -l) files)"
