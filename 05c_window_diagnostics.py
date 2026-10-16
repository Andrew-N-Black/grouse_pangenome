#!/usr/bin/env python3
"""
Step 05c: why does a stretch of a chromosome have no THEx trees? Per-window diagnostics.

    python3 05c_window_diagnostics.py <PROJECT_DIR>/work chr_Z [window_bp]

Standard library only (runs with the system python3). For every window of one chromosome it
reports, from files the pipeline already wrote:

    ref_N          fraction of the ptarmigan reference that is N (assembly gap)
    ref_masked     fraction soft-masked (lower case) in the RefSeq FASTA = repeats
    ref_GC         GC fraction of the ptarmigan sequence
    hap_cov_mean   mean fraction of the window aligned, across the grouse haplotypes
    haps_ge90      grouse haplotypes with >= 90% of the window aligned (the pairwise filter
                   needs every taxon at 90% after trimming)
    genes          RefSeq genes overlapping the window
    status         tree / dropped_trimal / dropped_pairwise_filter (from summarize)

and then prints the runs of consecutive windows without a tree with their averages, so a gap
can be read as an assembly gap (ref_N high), a repeat block / satellite (ref_masked high,
hap_cov low), a region the grouse lack or do not align (hap_cov low, ref_masked normal), or a
filter casualty (hap_cov high but dropped). Writes <work>/thex_inputs/<chrom>.window_diagnostics.tsv.
"""
import os
import sys

work, chrom = sys.argv[1], sys.argv[2]
W = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
og = open(os.path.join(work, "ref_id.txt")).read().strip()      # Lmuta = ptarmigan reference
proj = os.path.dirname(os.path.abspath(work))
ref_fa = os.path.join(proj, "inputs", "outgroup", og + ".chr.fasta")
if not os.path.exists(ref_fa):
    sys.exit("ERROR: %s not found (this script expects the ptarmigan reference, REF_ID=Lmuta)" % ref_fa)
msa = os.path.join(work, "msa", chrom + ".fasta")
if not os.path.exists(msa):
    msa = os.path.join(work, "msa_excluded", chrom + ".fasta")
gff = os.path.join(work, "thex_inputs", "reference_genes.gff")
status_f = os.path.join(work, "thex_inputs", "grouse_thex.window_status.tsv")
tchrom = chrom.replace("_", "-")


def read_record(path, name=None):
    """yield (name, seq) one record at a time; if name is given, return only that one."""
    cur, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and (name is None or cur == name):
                    yield cur, "".join(buf)
                    if name is not None:
                        return
                cur, buf = line[1:].split()[0], []
            elif name is None or cur == name:
                buf.append(line.strip())
    if cur is not None and (name is None or cur == name):
        yield cur, "".join(buf)


_, ref = next(read_record(ref_fa, chrom))
L = len(ref)
nw = (L + W - 1) // W
rows = []
for i in range(nw):
    s = ref[i * W:(i + 1) * W]
    n = len(s)
    up = s.upper()
    nN = up.count("N")
    acgt = n - nN
    rows.append({
        "start": i * W + 1, "end": min((i + 1) * W, L),
        "ref_N": nN / n,
        "ref_masked": sum(1 for c in s if c.islower()) / n,
        "ref_GC": (up.count("G") + up.count("C")) / acgt if acgt else 0.0,
        "covs": [],
    })
del ref

n_hap = 0
for name, seq in read_record(msa):
    if name == og:
        continue
    n_hap += 1
    for i, r in enumerate(rows):
        s = seq[i * W:(i + 1) * W]
        r["covs"].append(1 - s.count("-") / max(len(s), 1))

genes = [0] * nw
if os.path.exists(gff):
    with open(gff) as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] == tchrom:
                a, b = int(f[3]), int(f[4])
                for i in range((a - 1) // W, min((b - 1) // W, nw - 1) + 1):
                    genes[i] += 1
status = {}
if os.path.exists(status_f):
    with open(status_f) as fh:
        next(fh)
        for line in fh:
            c, s, e, st = line.rstrip("\n").split("\t")
            if c == tchrom:
                status[int(s)] = st

outp = os.path.join(work, "thex_inputs", chrom + ".window_diagnostics.tsv")
with open(outp, "w") as fh:
    fh.write("chrom\tstart\tend\tref_N\tref_masked\tref_GC\thap_cov_mean\thaps_ge90\tgenes\tstatus\n")
    for i, r in enumerate(rows):
        c = r["covs"]
        r["hap_cov_mean"] = sum(c) / len(c) if c else 0.0
        r["haps_ge90"] = sum(1 for x in c if x >= 0.9)
        r["genes"] = genes[i]
        r["status"] = status.get(r["start"], "no_window")
        fh.write("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%s\n" % (
            chrom, r["start"], r["end"], r["ref_N"], r["ref_masked"], r["ref_GC"],
            r["hap_cov_mean"], r["haps_ge90"], r["genes"], r["status"]))


def avg(rs, k):
    return sum(r[k] for r in rs) / len(rs)


tree = [r for r in rows if r["status"] == "tree"]
print("%s: %.1f Mb, %d windows, %d grouse haplotypes in the alignment" % (chrom, L / 1e6, nw, n_hap))
if tree:
    print("  windows WITH trees (baseline): ref_N %.3f  masked %.3f  GC %.3f  hap_cov %.3f  haps>=90%% %.1f  genes %.2f"
          % (avg(tree, "ref_N"), avg(tree, "ref_masked"), avg(tree, "ref_GC"),
             avg(tree, "hap_cov_mean"), avg(tree, "haps_ge90"), avg(tree, "genes")))
print("  runs of >= 3 consecutive windows without a tree:")
print("    %-24s %6s %6s %7s %6s %8s %9s %6s  %s" % ("region (Mb)", "n_win", "ref_N", "masked", "GC",
                                                  "hap_cov", "haps>=90", "genes", "statuses"))
run = []
for r in rows + [{"status": "tree"}]:
    if r["status"] != "tree":
        run.append(r)
        continue
    if len(run) >= 3:
        sts = {}
        for x in run:
            sts[x["status"]] = sts.get(x["status"], 0) + 1
        print("    %-24s %6d %6.3f %7.3f %6.3f %8.3f %9.1f %6.2f  %s" % (
            "%.2f-%.2f" % ((run[0]["start"] - 1) / 1e6, run[-1]["end"] / 1e6), len(run),
            avg(run, "ref_N"), avg(run, "ref_masked"), avg(run, "ref_GC"), avg(run, "hap_cov_mean"),
            avg(run, "haps_ge90"), avg(run, "genes"),
            ", ".join("%s:%d" % kv for kv in sorted(sts.items()))))
    run = []
print("  table: %s" % outp)
