# After the misjoin fix: what's left to do

Written 2026-10-13. State when written:
- 24 RagTag misjoins (chr_20 in F5544 hap2, chr_23-29 in 18 haplotypes) were cut in the
  **pangenome inputs** (`pangenome_inputs/`); originals are in `pangenome_inputs/pre_misjoin_fix/`.
- Round-2 chrom_qc is clean on chr_20 and chr_23-29.
- 9 per-chromosome pggb rebuilds were submitted (jobs 17002918-17002926: chr_20, chr_23 ... chr_30).
- The **deliverable assemblies** (`final_ptarmigan/`) and their Liftoff GFF3s still contain the misjoins.

Do NOT re-run grouse_pangenome step 01 (`01_break_unsupported_joins.sh`): it rewrites
`pangenome_inputs/` (undoing the fix), and after Part B its AGP check will no longer match
`final_ptarmigan/`.

Set these on every Gautschi login (edit the first two):

```bash
REPO=/path/to/grouse_pangenome
PG=/path/to/prairie_grouse
ANCHOR=/full/path/to/the/anchor.fa
ASM=$CLUSTER_SCRATCH/GROUSE/grouse_asm
QC=$REPO/scripts/chrom_qc
```

Part B does not depend on the graphs and can start any time.

---

## Part A: when the 9 graph jobs finish

### A1. Did they all finish?

```bash
sacct -j 17002918,17002919,17002920,17002921,17002922,17002923,17002924,17002925,17002926 \
      --format=JobID,JobName%12,State,Elapsed,MaxRSS | grep -v '\.'
ls $ASM/pangenome/chr_{20,23,24,25,26,27,28,29,30}/*.final.gfa
```

All 9 should show COMPLETED and a `*.final.gfa`. If one failed on time, resubmit just that one:

```bash
cd $REPO
CHROM=chr_NN PGGB_EXTRA="--resume" sbatch --export=ALL -A dewoody -t 4-00:00:00 --mem=128G -J pg_chr_NN 03_pangenome_build.sh
```

### A2. Check chr_20 is fixed

```bash
grep -c 'chr_26_p1' $ASM/pangenome/chr_26/combined_haplotypes.fa.gz.fai    # expect 1
ls $ASM/pangenome/chr_20/*.png                                             # pggb's viz + layout images
```

Copy the chr_20 `*.viz*.png` and `*.lay*.png` to the laptop and compare with the old ones:
- one open line instead of a loop;
- no long links spanning the graph;
- no block owned by only LEPC F5544#2;
- the ~0.4 Mb presence/absence insertion still visible as a bubble.

If STGR F5457#1 still sits apart from the others, that is a graph problem, not an assembly
problem: chrom_qc found its chr_20 normal.

### A2b. Leave chr_MT out of the merged graph (2026-10-14)

The hifiasm chr_MT sequences are multi-copy (20-50 kb). The mitogenomes were
re-assembled separately with MitoHiFi (prairie_grouse `accessory_scripts/mitohifi_array.sh`
and `mito_qc.sh`) and are kept out of the pangenome. The pangenome input FASTAs are not
edited for this.

Only once NO 03_pangenome_build array task is still pending or running, because the
array reads `chromosomes.txt` by line number and removing a line shifts every later
task onto the wrong chromosome:

```bash
squeue -u $USER -n grouse_pangenome,pg_chr_20,pg_chr_23,pg_chr_24,pg_chr_25,pg_chr_26,pg_chr_27,pg_chr_28,pg_chr_29,pg_chr_30   # must be empty
cp $ASM/pangenome/chromosomes.txt $ASM/pangenome/chromosomes.with_MT.txt
sed -i '/^chr_MT$/d' $ASM/pangenome/chromosomes.txt
grep -c chr_MT $ASM/pangenome/chromosomes.txt      # expect 0
```

### A3. Merge all per-chromosome graphs

Wait until every chromosome graph is done (including the untouched ones from the original
arrays), and until Part B (B5) is finished, because the summary's gene tables read
`final_ptarmigan` + Liftoff.

```bash
cd $REPO
STAGE=summary sbatch --export=ALL -t 4-00:00:00 03_pangenome_build.sh
```

---

## Part B: apply the same fix to the deliverable assemblies (prairie_grouse)

The pangenome-input coordinates differ from `final_ptarmigan` (step 01 removed other gaps), so
`cuts.tsv` is not reused directly. The same two scripts are re-run on the final FASTAs instead.

Affected samples (13; manifest array indices in brackets):
F5457 [0], F5462 [1], F5468 [3], F5478 [5], F5502 [8], F5503 [9], F5541 [11], F5542 [12],
F5543 [13], F5544 [14], F5545 [15], F5598 [20], F5599 [21].

### B1. Map the 18 affected final assemblies to the anchor

```bash
mkdir -p $QC/final_fix && cd $QC/final_fix
ls $ASM/final_ptarmigan/*_hap[12].pseudo_chr.fasta > haps_final_all.txt     # all 46: medians for length checks
wc -l < haps_final_all.txt                                                  # expect 46
sed "s#^#$ASM/final_ptarmigan/#; s#\$#.pseudo_chr.fasta#" $QC/fixed_haps.txt > haps_final_changed.txt
wc -l < haps_final_changed.txt                                              # expect 18
sbatch -A dewoody -p cpu -q standby -t 04:00:00 -n 16 --mem=32G \
  --array=1-18 $QC/run_chrom_qc.sh map $ANCHOR haps_final_changed.txt
```

### B2. Dry run: it should find the same 24 joins

```bash
cd $QC/final_fix
python3 $QC/fix_misjoins.py --haps haps_final_all.txt --pafdir chromqc_paf --outdir fixed --dry-run
cut -f1,2,7,8 fixed/cuts.tsv
cut -f1,2,7,8 $QC/fixed/cuts.tsv
```

Expect 24 rows naming the same haplotypes, chromosomes and pieces as the pangenome run (the
coordinates will differ). The 28 unchanged haplotypes print "no PAF ... skipped"; that is
expected, as they are only there for the length medians. Stop and check if the lists differ.

### B3. Real run, then swap into final_ptarmigan (originals backed up)

```bash
cd $QC/final_fix
sbatch -A dewoody -p cpu -q standby -t 04:00:00 -n 4 --mem=16G -J fixfinal -o fixfinal.log \
  --wrap "python3 $QC/fix_misjoins.py --haps haps_final_all.txt --pafdir chromqc_paf --outdir fixed"
```

When that job is done:

```bash
module load biocontainers samtools
cd $ASM/final_ptarmigan
mkdir -p pre_misjoin_fix
while read h; do
  cp -p ${h}.pseudo_chr.fasta* ${h}.liftoff.gff3 pre_misjoin_fix/ 2>/dev/null
  cp $QC/final_fix/fixed/${h}.pseudo_chr.fixed.fa ${h}.pseudo_chr.fasta.tmp
  mv ${h}.pseudo_chr.fasta.tmp ${h}.pseudo_chr.fasta
  samtools faidx ${h}.pseudo_chr.fasta
done < $QC/fixed_haps.txt
for h in $(cat $QC/fixed_haps.txt); do
  echo -e "$h\t$(( $(awk '{s+=$2}END{print s}' pre_misjoin_fix/${h}.pseudo_chr.fasta.fai) - $(awk '{s+=$2}END{print s}' ${h}.pseudo_chr.fasta.fai) ))"
done
```

The last loop should print 100 for every haplotype except STGR_F5462_hap2 (200), as before.

### B4. Regenerate QUAST, Liftoff and dot plots (prairie_grouse step 05)

Step 05 is timestamp-driven:
- the RagTag and renaming stages see that the final FASTA is newer and skip, so the fixed
  FASTA is not overwritten;
- QUAST, Liftoff and the dot plot are older than the new FASTA, so they rerun;
- unchanged haplotypes skip everything.

Hi-C and depth are off here to keep it fast (see B7).

```bash
cd $PG && git pull
sbatch --export=ALL,DO_HIC=false,DO_DEPTH=false,DO_TIDK=false \
  --array=0,1,3,5,8,9,11,12,13,14,15,20,21 05_ragtag_liftoff_array.sh
```

In each task log, check that the rename stage says "already renamed — skipping", and that
QUAST and Liftoff actually ran. Then:

```bash
bash 05b_relabel_male_W.sh     # safe to rerun; should report nothing to change
```

### B5. Regenerate BUSCO and the BUSCO figures (steps 07, 09)

Step 08 (ptarmigan) and step 10 (sex check) are unaffected; Z and W were not touched.

```bash
cd $PG
J7=$(sbatch --parsable --array=0,1,3,5,8,9,11,12,13,14,15,20,21 07_busco_array.sh)
sbatch --dependency=afterok:$J7 09_busco_plots.sh
```

BUSCO completeness will not change, because the sequence is identical. What changes is
which chromosome each ortholog sits on, so Figures 13-14 (which include STGR F5457 hap1)
must be redrawn.

### B6. Update the contiguity figure (Fig. S5)

`accessory_scripts/plot_assembly_contiguity.R` holds the values inline. Pull the new numbers
for the 18 haplotypes, then edit those rows and re-plot locally:

```bash
for h in $(cat $QC/fixed_haps.txt); do
  s=$(echo $h | cut -d_ -f2); p=$(echo $h | cut -d_ -f3)
  f=$ASM/qc_ptarmigan/quast/${s}.${p}.post_ragtag/report.tsv
  echo "== $h"; grep -E '^(# contigs|Total length|N50)' $f
done
```

Expect +1 or 2 scaffolds per haplotype, total length -100 bp, and N50 unchanged.

### B7. Optional, before NCBI submission

- **Contact maps and join scores:** rerun step 05 with `DO_HIC=true` for the same indices,
  then step 11 (`--array=` same indices). This gives Pretext maps of the corrected
  assemblies and join-support scores that match them.
- **AGP:** the RagTag AGP in `ragtag_ptarmigan/` no longer matches the 18 corrected FASTAs,
  and NCBI needs an AGP that does. Regenerate it from the corrected FASTA (split at N runs)
  at submission time.

---

## Report / manuscript updates

- Figures 13-14: redraw from B5 (STGR F5457 hap1 carried a chr_23+chr_24 misjoin).
- Fig. S5 (contiguity): update from B6.
- Fig. S6 (BUSCO completeness): values unchanged; regenerate only for consistency.
- Text: the "apparent chr_23-chr_24 merger ... most likely an assembly artifact" is now
  confirmed as RagTag scaffolding misjoins at gaps, and wider than first described:
  - 24 joins in 18 haplotypes;
  - mostly between ptarmigan chr_23-chr_27, plus chr_26 attached to chr_20 in LEPC F5544
    hap2 and chr_29 attached to chr_28 in STGR F5462 hap1;
  - all detected by alignment to an anchor haplotype and corrected before graph
    construction.
