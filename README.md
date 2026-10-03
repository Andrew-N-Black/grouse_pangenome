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

Specifically, step 01 reads, for every sample in `assembly_manifest.tsv`:

```
$CLUSTER_SCRATCH/GROUSE/grouse_asm/final_ptarmigan/<SPECIES>_<SAMPLE>_hap{1,2}.pseudo_chr.fasta
$CLUSTER_SCRATCH/GROUSE/grouse_asm/final_ptarmigan/<SPECIES>_<SAMPLE>_hap{1,2}.liftoff.gff3
$CLUSTER_SCRATCH/GROUSE/grouse_asm/liftoff_ptarmigan/<SAMPLE>.<HAP>.unmapped_features.txt
```

so the assembly pipeline must have completed through **its step 05** first. Both
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
| 01 | [`01_pangenome_build.sh`](01_pangenome_build.sh) | Renames every haplotype to PanSN (`species.sample#hap#contig`), concatenates them, and builds one PGGB graph. Then: `odgi` graph stats and per-haplotype path lengths, `panacus` core/accessory growth by species, de novo RepeatModeler library per species → RepeatMasker on all 46 haplotypes, and a Python pass that joins genome length + Liftoff gene counts + TE content into per-haplotype and per-species tables. | single job, 64 cpu, 500 G, 14 d |
| 02 | [`02_shortread_introgression.sh`](02_shortread_introgression.sh) | Converts the graph to GBZ + giraffe indexes, then per short-read sample: `vg giraffe` → `vg pack` → `vg call` per reference contig → concatenate. Merges to a joint VCF with `bcftools`, then `Dsuite Dtrios` for genome-wide ABBA-BABA D-statistics and `Dsuite Dinvestigate` for windowed D/f_d on the trio the data supports. | single job, 64 cpu, 250 G, 10 d |

Both are single jobs rather than arrays, because a pangenome is inherently a
joint analysis. The embarrassingly parallel sub-steps (RepeatMasker per genome,
mapping per sample) are parallelised internally with `xargs -P`.

```bash
sbatch 01_pangenome_build.sh
# then, once the graph exists:
sbatch 02_shortread_introgression.sh
```

### Honest scope note on step 01

Whole-genome PGGB across ~46 haplotypes of ~1.1 Gb bird genomes is a heavy
alignment job — potentially multi-day to multi-week. The 14-day/500 G request is
a starting point, not a guarantee. If it does not finish in walltime, use pggb's
`partition-before-pggb` helper to run per chromosome; because every assembly is
already RagTag-oriented against the same reference, chromosome-level
partitioning should partition cleanly.

Step 01 is resumable at its two expensive stages: an existing `*.final.gfa`
skips the pggb run, and an existing `.tbl` skips RepeatMasker for that
haplotype. `panacus` failure is caught and warned about rather than fatal, since
the TE and gene comparisons do not depend on it.

---

## Manifests

### `assembly_manifest.tsv`
Copied from the assembly repo and **must stay in sync with it** — step 01 builds
assembly filenames from its `sample_id` and `species` fields verbatim.

```
sample_id  species  hifi_bams  hic_r1  hic_r2  ont_ul
```

Only the first two columns are used here; the read paths are vestigial but kept
so the file can be diffed against the assembly repo's copy.

> **Known issue in the current file:** the seven LEPC rows have a **trailing
> space** in the `species` field (`"LEPC "`, not `"LEPC"`). Step 01 therefore
> looks for `LEPC _F5540_hap1.pseudo_chr.fasta`, with the space. That works only
> if the assembly pipeline ran with the same whitespace — which it did, since it
> reads the same manifest — but it is fragile and it makes the PanSN sample name
> `LEPC_.F5540`. Step 01 prints a warning listing the affected rows and does
> **not** silently trim them, because trimming would break paths that currently
> resolve. Fix it in both repos together, or leave it deliberately.

### `shortread_manifest.tsv`
One row per population-resequencing individual — there can be many more of these
than the 23 assembled genomes.

```
sample_id  species  fastq_r1  fastq_r2
```

`species` must be one of `LEPC`, `GRPC`, `STGR`, or the literal `Outgroup`.

**At least one `Outgroup` sample is required and step 02 refuses to run without
one.** D-statistics need a fourth, more distantly related taxon to polarise
ancestral versus derived alleles; without it, "quantify introgression" is not a
defined computation. Chicken resequencing data is the obvious choice.

---

## Two deliberate non-decisions

**No hardcoded species tree.** With three ingroup species there is exactly one
possible trio, and we do not have independent knowledge of these species'
phylogeny or of which pairs are sympatric. Rather than invent either, step 02
uses Dsuite's untreed mode: it infers the best-supported P1/P2/P3 arrangement
from allele-sharing patterns in the data (the BBAA pattern), and that
arrangement's P2/P3 pair is the one showing excess sharing. Cross-check the
result against known species ranges — a statistically supported pair that is not
geographically sympatric is a reason to doubt the signal, not to report it.

**Reference haplotype for VCF coordinates is not hand-picked.** Step 02 anchors
calling on `hap1` of the first sample in `assembly_manifest.tsv`, so it always
points at a path that really exists in the graph. Override `REF_HAP` /
`REF_SAMPLE_TAG` near the top if you would rather anchor on a specific
higher-quality individual.

---

## Output layout

Under `$CLUSTER_SCRATCH/GROUSE/grouse_asm` (shared with the assembly pipeline):

```
pangenome/              combined_haplotypes.fa.gz, PGGB graph (*.final.gfa/.og)  (01)
te_annotation/          RepeatModeler libraries, per-haplotype RepeatMasker      (01)
comparison/             graph summary, path lengths, per-haplotype and
                        per-species tables, panacus growth curves               (01)
panacus-0.5.1/          cached precompiled binary                               (01)
graph_index/            GBZ + giraffe indexes, reference contig list            (02)
mapped/                 per-sample .gam and .pack                               (02)
called/                 per-sample VCFs (per contig, then concatenated)         (02)
joint_vcf/joint.vcf.gz  merged multi-sample VCF                                 (02)
introgression/          SETS.tsv, Dtrios BBAA/Dmin, windowed localFstats         (02)
```

The headline tables are
`comparison/species_comparison_summary.tsv` (length, gene content and TE content
per species) and `introgression/dtrios_BBAA.txt` (genome-wide D).

---

## Tooling caveats carried over

These were noted when the scripts were written and have not been re-verified
against the current module stack:

- `pggb`, `odgi`, `vg` and `dsuite` module availability varies by cluster. All
  four are also available via bioconda or as Singularity containers; `pggb`
  especially is often easier as a container.
- `panacus` has no module and the cluster has no Rust toolchain, so step 01
  downloads the precompiled musl-static binary from GitHub releases and caches it
  under `PROJECT_DIR`. It downloads once.
- `vg gbwt -G` is used to make the GBZ rather than `vg autoindex` reading the
  PGGB GFA directly: PGGB writes all paths as P-lines and autoindex's GFA loader
  expects W-lines for per-haplotype paths
  ([vg#4302](https://github.com/vgteam/vg/issues/4302)).
- `vg pack` and `vg paths` flags come from long-stable usage rather than a
  freshly checked `--help`. If either errors, check `--help` first.
