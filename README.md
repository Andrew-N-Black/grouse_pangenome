# Prairie Grouse Pangenome and Introgression

SLURM/bash pipeline (Purdue RCAC, Gautschi) building **one pangenome graph
across all prairie grouse haplotype assemblies**, then using it to compare the
three species and to quantify introgression from population-level short-read
data.

| Species | Code | Scientific name |
|---|---|---|
| Lesser Prairie-Chicken | LEPC | *Tympanuchus pallidicinctus* |
| Greater Prairie-Chicken | GRPC | *T. cupido* |
| Sharp-tailed Grouse | STGR | *T. phasianellus* |

## Relationship to the assembly pipeline

**This repo starts where [`prairie_grouse`](https://github.com/Andrew-N-Black/prairie_grouse)
finishes.** That repo assembles and QCs the 23 individuals (46 haplotypes);
this one consumes its output. Nothing here re-runs an assembler.

Specifically, steps 01–03 read, for every sample in `assembly_manifest.tsv`:

```
$CLUSTER_SCRATCH/GROUSE/grouse_asm/final_ptarmigan/<SPECIES>_<SAMPLE>_hap{1,2}.pseudo_chr.fasta
$CLUSTER_SCRATCH/GROUSE/grouse_asm/final_ptarmigan/<SPECIES>_<SAMPLE>_hap{1,2}.liftoff.gff3
$CLUSTER_SCRATCH/GROUSE/grouse_asm/final_ptarmigan/W_relabel_log.tsv
$CLUSTER_SCRATCH/GROUSE/grouse_asm/ragtag_ptarmigan/<SAMPLE>.<HAP>/ragtag.scaffold.agp
$CLUSTER_SCRATCH/GROUSE/grouse_asm/ragtag_ptarmigan/<SAMPLE>.<HAP>/<SAMPLE>.<HAP>.rename_map.tsv
$CLUSTER_SCRATCH/GROUSE/grouse_asm/qc_ptarmigan/join_support/<SPECIES>_<SAMPLE>_<HAP>.joins.tsv
$CLUSTER_SCRATCH/GROUSE/grouse_asm/liftoff_ptarmigan/<SAMPLE>.<HAP>.unmapped_features.txt
```

so the assembly pipeline must have completed through its Hi-C join-support
step (`13_join_support_ptarmigan.sh`), its chromosome renaming (step 14) and
`relabel_male_W.sh` first. Both
repos share the same `PROJECT_DIR`, and the `_ptarmigan` suffix comes from
`REF_TAG`, which must match between them.

### Chromosome names

The assemblies are ordered against the rock ptarmigan reference (*Lagopus
muta*, bLagMut1, `GCF_023343835.1`), which numbers its chromosomes **by
descending length, not by homology** — `chr_6` here is not chicken chromosome 6.
Nothing in this pipeline depends on that, because every assembly in the graph
uses the same convention. But any result reported *per chromosome* has to be
translated through the ptarmigan↔chicken correspondence table that the assembly
repo's step 04 produces before it can be set against the galliform literature.
Three chromosomes genuinely differ between chicken and these birds (chicken 6+8
are one chromosome here; GGA4p is free-standing here; chicken 2 corresponds to
two chromosomes), so this is not a formality.

---

## Pipeline

| Step | Script | What it does | Shape |
|---|---|---|---|
| 01 | [`01_break_unsupported_joins.sh`](01_break_unsupported_joins.sh) (+ [`.py`](01_break_unsupported_joins.py)) | Cuts every RagTag join that the assembly repo's Hi-C join-support step called `unsupported`, removing its 100 bp gap, and writes one graph-input FASTA per haplotype. Before cutting it proves the AGP coordinates still match each FASTA (object length; the cut gap must be all N) and finds the joins table's position column from the data; any mismatch stops the run. A cut chromosome becomes `<chrom>_p1`, `<chrom>_p2`, …; `breaks.tsv`, `pieces.tsv` and `summary.tsv` record every cut and the bp accounting. Joins RagTag made on the reference's *unplaced* scaffolds (objects named `scaffold_N`) were never scored by the join-support step and are cut too (`CUT_UNPLACED=false` keeps them). `BREAK=all` returns every chromosome to its yahs scaffolds. | single job, < 1 h |
| 02 | [`02_assembly_qc_array.sh`](02_assembly_qc_array.sh) | Per sample, on step 01's FASTAs: Merqury (k = 21 HiFi k-mers: QV, completeness, copy-number spectra); HiFi depth with reads mapped to hap1 + hap2 together, flagging 100 kb windows below 0.5× or above 1.5× the haplotype median to a BED; optional NCBI FCS-adaptor. | array, one task per sample |
| 02b | [`02b_collate_qc.sh`](02b_collate_qc.sh) | One table, one row per haplotype, with a `REVIEW` column (QV < 40; completeness > 5 points below the panel median; > 5% of the haplotype in low- or high-depth windows). Triage thresholds, not pass/fail. | login node, seconds |
| 03 | [`03_pangenome_build.sh`](03_pangenome_build.sh) | Renames every haplotype to PanSN (`species.sample#hap#contig`), concatenates them, and builds one PGGB graph. Then: `odgi` graph stats and per-haplotype path lengths, `panacus` core/accessory growth by species, de novo RepeatModeler library per species → RepeatMasker on all 46 haplotypes, and a Python pass that joins genome length + Liftoff gene counts + TE content into per-haplotype and per-species tables. | single job, 64 cpu, 500 G, 14 d |
| 04 | [`04_shortread_introgression.sh`](04_shortread_introgression.sh) | Converts the graph to GBZ + giraffe indexes, then per short-read sample: `vg giraffe` → `vg pack` → `vg call` per reference contig → concatenate. Merges to a joint VCF with `bcftools`, then `Dsuite Dtrios` for genome-wide ABBA-BABA D-statistics and `Dsuite Dinvestigate` for windowed D/f_d on the trio the data supports. | single job, 64 cpu, 250 G, 10 d |

Steps 03 and 04 are single jobs rather than arrays, because a pangenome is
inherently a joint analysis. The embarrassingly parallel sub-steps (RepeatMasker per genome,
mapping per sample) are parallelised internally with `xargs -P`.

```bash
# run everything from the repo root, so assembly_manifest.tsv is found
DRY_RUN=true bash 01_break_unsupported_joins.sh     # plan only: check the cuts first
sbatch 01_break_unsupported_joins.sh
N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
sbatch --array=0-$((N-1))%6 02_assembly_qc_array.sh
bash 02b_collate_qc.sh
# manual review of flagged regions, then:
sbatch 03_pangenome_build.sh
# then, once the graph exists:
sbatch 04_shortread_introgression.sh
```

### Why the joins are cut before the graph (step 01)

RagTag joined each haplotype's yahs scaffolds in the order the rock ptarmigan
reference implied. A graph builder places sequences by alignment and gains
nothing from that order, but an incorrect join becomes a false structural
variant. So only joins with Hi-C support survive into the graph. The TE and
gene comparisons in step 03 still read the uncut `final_ptarmigan/` assemblies:
cutting a join removes only a 100 bp gap, and the Liftoff GFF3 seqids name the
uncut chromosomes.

### Notes on step 02

HiFi reads come from the manifest's BAMs unless `HIFI_READS_GLOB` points to
adapter-filtered reads. FCS-GX (contamination) is deliberately not included:
its ~470 GB database needs a large-memory node, and NCBI runs it at submission.
Sex-chromosome-depleted female haplotypes are expected to show lower k-mer
completeness. Steps 01 and 02 are resume-safe.

### Honest scope note on step 03

Whole-genome PGGB across ~46 haplotypes of ~1.1 Gb bird genomes is a heavy
alignment job — potentially multi-day to multi-week. The 14-day/500 G request is
a starting point, not a guarantee. If it does not finish in walltime, use pggb's
`partition-before-pggb` helper to run per chromosome; because every assembly is
already RagTag-oriented against the same reference, chromosome-level
partitioning should partition cleanly.

Step 03 is resumable at its two expensive stages: an existing `*.final.gfa`
skips the pggb run, and an existing `.tbl` skips RepeatMasker for that
haplotype. `panacus` failure is caught and warned about rather than fatal, since
the TE and gene comparisons do not depend on it.

---

## Manifests

### `assembly_manifest.tsv`
Copied from the assembly repo and **must stay in sync with it** — step 03 builds
assembly filenames from its `sample_id` and `species` fields verbatim.

```
sample_id  species  hifi_bams  hic_r1  hic_r2  ont_ul
```

Steps 03–04 use only the first two columns. Step 02 also reads `hifi_bams`
(comma-separated if a sample has several) for Merqury and the depth check.

> **Fixed:** the seven LEPC rows previously had a trailing space in the
> `species` field (`"LEPC "`). The assembly files on disk are named without it
> (`LEPC_F5540_hap1…`), so the space would have made step 03 miss every LEPC
> haplotype. It has been trimmed; step 03 still warns if whitespace reappears.

### `shortread_manifest.tsv`
One row per population-resequencing individual — there can be many more of these
than the 23 assembled genomes.

```
sample_id  species  fastq_r1  fastq_r2
```

`species` must be one of `LEPC`, `GRPC`, `STGR`, or the literal `Outgroup`.

**At least one `Outgroup` sample is required and step 04 refuses to run without
one.** D-statistics need a fourth, more distantly related taxon to polarise
ancestral versus derived alleles; without it, "quantify introgression" is not a
defined computation. Chicken resequencing data is the obvious choice.

---

## Two deliberate non-decisions

**No hardcoded species tree.** With three ingroup species there is exactly one
possible trio, and we do not have independent knowledge of these species'
phylogeny or of which pairs are sympatric. Rather than invent either, step 04
uses Dsuite's untreed mode: it infers the best-supported P1/P2/P3 arrangement
from allele-sharing patterns in the data (the BBAA pattern), and that
arrangement's P2/P3 pair is the one showing excess sharing. Cross-check the
result against known species ranges — a statistically supported pair that is not
geographically sympatric is a reason to doubt the signal, not to report it.

**Reference haplotype for VCF coordinates is not hand-picked.** Step 04 anchors
calling on `hap1` of the first sample in `assembly_manifest.tsv`, so it always
points at a path that really exists in the graph. Override `REF_HAP` /
`REF_SAMPLE_TAG` near the top if you would rather anchor on a specific
higher-quality individual.

---

## Output layout

Under `$CLUSTER_SCRATCH/GROUSE/grouse_asm` (shared with the assembly pipeline):

```
pangenome_inputs/       <SP>_<S>_<H>.pg_input.fasta (+.fai), breaks.tsv,
                        pieces.tsv, summary.tsv                                 (01)
pangenome_qc/<SAMPLE>/  Merqury (QV, completeness, spectra-cn), depth summary,
                        depth_flags.bed, optional fcs_adaptor/                  (02)
pangenome_qc/qc_summary.tsv   one row per haplotype, with REVIEW reasons       (02b)
pangenome/              combined_haplotypes.fa.gz, PGGB graph (*.final.gfa/.og)  (03)
te_annotation/          RepeatModeler libraries, per-haplotype RepeatMasker      (03)
comparison/             graph summary, path lengths, per-haplotype and
                        per-species tables, panacus growth curves               (03)
panacus-0.5.1/          cached precompiled binary                               (03)
graph_index/            GBZ + giraffe indexes, reference contig list            (04)
mapped/                 per-sample .gam and .pack                               (04)
called/                 per-sample VCFs (per contig, then concatenated)         (04)
joint_vcf/joint.vcf.gz  merged multi-sample VCF                                 (04)
introgression/          SETS.tsv, Dtrios BBAA/Dmin, windowed localFstats         (04)
```

The headline tables are
`comparison/species_comparison_summary.tsv` (length, gene content and TE content
per species) and `introgression/dtrios_BBAA.txt` (genome-wide D).

---

## Tooling caveats carried over

These were noted when the scripts were written and have not been re-verified
against the current module stack:

- `merqury` and `minimap2` (step 02) are Gautschi biocontainer modules; step 02
  stops at the start and names any that are missing. Gautschi has **no mosdepth
  module** and its merqury module does not expose `meryl`, so step 02 uses the
  official release binaries (mosdepth v0.3.14, meryl v1.4.1), cached under
  `PROJECT_DIR/tools/`. Fetch it once on a login node, where there is
  network access: `INSTALL_ONLY=true bash 02_assembly_qc_array.sh`.
- `pggb`, `odgi`, `vg` and `dsuite` module availability varies by cluster. All
  four are also available via bioconda or as Singularity containers; `pggb`
  especially is often easier as a container.
- `panacus` has no module and the cluster has no Rust toolchain, so step 03
  downloads the precompiled musl-static binary from GitHub releases and caches it
  under `PROJECT_DIR`. It downloads once.
- `vg gbwt -G` is used to make the GBZ rather than `vg autoindex` reading the
  PGGB GFA directly: PGGB writes all paths as P-lines and autoindex's GFA loader
  expects W-lines for per-haplotype paths
  ([vg#4302](https://github.com/vgteam/vg/issues/4302)).
- `vg pack` and `vg paths` flags come from long-stable usage rather than a
  freshly checked `--help`. If either errors, check `--help` first.
