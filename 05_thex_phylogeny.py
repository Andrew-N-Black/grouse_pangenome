#!/usr/bin/env python3
"""
Helper for 05_thex_phylogeny.sh: builds THEx (Tree House Explorer) inputs from the
prairie grouse haplotype assemblies. Subcommands:

  outgroup-ref  rock ptarmigan as the coordinate reference, renamed to the grouse chr_ names (default)
  pick-ref   alternatively, choose a grouse haplotype as the reference from the .fai files
  project    project one assembly onto the reference's coordinates (from a minimap2 PAF with cs tags)
  msa        gather the projections into one multiple alignment per reference chromosome
  annot      THEx annotation track (genes only) + chromosome-lengths BED, in THEx chromosome names
  summarize  IQ-TREE window trees -> rooted Tree Viewer file with a species-level TopologyID

The alignment is reference-projected (as in Foley et al. 2021, which THExBuilder's filters
come from): every haplotype is aligned to ONE reference haplotype, and its bases are written
at the reference coordinate they align to. Insertions relative to the reference are dropped,
deletions and unaligned reference positions become '-', and reference positions covered by
more than one alignment of the same query (collapsed duplicates, paralogs) are set to '-'
rather than guessed. Every sequence in a chromosome's alignment therefore has exactly the
reference chromosome's length, which is what THExBuilder's windower assumes.

Standard library only, except `summarize` (pandas, openpyxl, ete3: all in the THEx env).
"""

import argparse
import bisect
import glob
import os
import re
import shutil
import statistics
import sys
from collections import defaultdict

# ----------------------------------------------------------------------- utils --
# A, C, G, T (either case) -> upper case; anything else (N, IUPAC, '-') -> '-'
_TT = bytearray(b"-" * 256)
for _c in b"ACGT":
    _TT[_c] = _c
    _TT[_c + 32] = _c
TRANSLATE = bytes(_TT)

CS_RE = re.compile(rb"(:[0-9]+|\*[a-zA-Z][a-zA-Z]|[+\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2})")


def thex_name(chrom):
    """THExBuilder's windower turns '_' into '-' in chromosome names (chr_6 -> chr-6)."""
    return chrom.replace("_", "-")


def read_fai(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


def read_fasta(path, keep=None):
    """name -> bytes (raw case). keep: optional set of names to load."""
    seqs, name, buf = {}, None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None and (keep is None or name in keep):
                    seqs[name] = b"".join(buf)
                name = line[1:].split()[0].decode()
                buf = []
            elif keep is None or name in keep:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        seqs[name] = b"".join(buf)
    return seqs


def read_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


# -------------------------------------------------------------------- pick-ref --
def cmd_pick_ref(a):
    """Pick the haplotype whose chromosome lengths are most typical of the panel.

    A RagTag misjoin makes one chromosome much longer than the panel median (e.g. the
    chr_23 + chr_24 merger) and leaves another short or missing, so the haplotype with the
    smallest worst-case deviation from the per-chromosome median is the safest coordinate
    system. Haplotypes missing chr_Z (the W-bearing haplotype of a ZW female) are not
    eligible, otherwise the Z would have no coordinates at all.
    """
    ids = read_list(a.ids)
    fais = {}
    for i in ids:
        p = os.path.join(a.fai_dir, i + ".pseudo_chr.fasta.fai")
        if not os.path.exists(p):
            sys.exit("ERROR: missing %s (run the transfer stage first)" % p)
        fais[i] = read_fai(p)
    excl = set(a.exclude.split(","))
    # The comparison set is every chr_* carried by at least half the haplotypes, so one
    # haplotype's odd extra sequence (or a W-bearing haplotype's chr_W) cannot make every
    # other haplotype "incomplete".
    counts = defaultdict(int)
    for f in fais.values():
        for c in f:
            if c.startswith("chr_") and c not in excl:
                counts[c] += 1
    chroms = sorted(c for c, n in counts.items() if n >= len(fais) / 2.0)
    med = {c: statistics.median([f[c] for f in fais.values() if c in f]) for c in chroms}
    zmax = max([f.get("chr_Z", 0) for f in fais.values()] or [0])
    rows = []
    for i, f in fais.items():
        missing = [c for c in chroms if c not in f]
        devs = {c: abs(f[c] / med[c] - 1) for c in chroms if c in f}
        worst_c = max(devs, key=devs.get)
        has_z = f.get("chr_Z", 0) >= 0.8 * zmax if zmax else True
        eligible = has_z and not missing and "chr_W" not in f
        rows.append((not eligible, devs[worst_c], i, worst_c, len(missing), has_z, "chr_W" in f))
    rows.sort()
    with open(a.out_table, "w") as fh:
        fh.write("haplotype\teligible\tmax_rel_dev\tworst_chrom\tn_chrom_missing\thas_full_Z\thas_W\n")
        for r in rows:
            fh.write("%s\t%s\t%.4f\t%s\t%d\t%s\t%s\n" % (r[2], not r[0], r[1], r[3], r[4], r[5], r[6]))
    best = rows[0]
    if best[0]:
        sys.exit("ERROR: no haplotype has every common chromosome, a full-length chr_Z and no chr_W "
                 "(see %s); set REF_ID by hand" % a.out_table)
    print(best[2])


# ------------------------------------------------------------- outgroup-ref --
def cmd_outgroup_ref(a):
    """Rock ptarmigan as the coordinate reference, with its chromosomes renamed exactly the
    way the assembly pipeline named the grouse chromosomes (chr_ + the assembly report's
    name, from <ASM>.chr_map.tsv: accession -> name). Only assembled molecules are kept.
    Writes FASTA (80 columns) + .fai, and the RefSeq GFF restricted to those sequences."""
    amap = {}
    with open(a.chr_map) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#"):
                amap[f[0]] = "chr_" + f[1]
    if not amap:
        sys.exit("ERROR: empty chromosome map %s" % a.chr_map)
    tmp = a.out_fasta + ".part"
    fai = []
    offset = 0

    def flush(name, chunks, oh):
        nonlocal offset
        seq = b"".join(chunks)
        hdr = b">" + name.encode() + b"\n"
        oh.write(hdr)
        offset += len(hdr)
        fai.append((name, len(seq), offset))
        for i in range(0, len(seq), 80):
            oh.write(seq[i:i + 80] + b"\n")
        offset += len(seq) + (len(seq) + 79) // 80

    with open(a.fasta, "rb") as ih, open(tmp, "wb") as oh:
        name, chunks = None, []
        for line in ih:
            if line.startswith(b">"):
                if name is not None:
                    flush(name, chunks, oh)
                acc = line[1:].split()[0].decode()
                name, chunks = amap.get(acc), []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            flush(name, chunks, oh)
    if len(fai) != len(amap):
        sys.exit("ERROR: %d of %d chromosome-map accessions found in %s" % (len(fai), len(amap), a.fasta))
    os.replace(tmp, a.out_fasta)
    with open(a.out_fasta + ".fai", "w") as fh:
        for n, L, off in fai:
            fh.write("%s\t%d\t%d\t80\t81\n" % (n, L, off))
    n = 0
    if a.gff:
        with open(a.gff) as ih, open(a.out_gff, "w") as oh:
            for line in ih:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) == 9 and f[0] in amap:
                    f[0] = amap[f[0]]
                    oh.write("\t".join(f))
                    n += 1
    sys.stderr.write("[outgroup-ref] %d chromosomes (%s ... %s), %d GFF lines\n"
                     % (len(fai), fai[0][0], fai[-1][0], n))


# --------------------------------------------------------------------- project --
def mask_multi_covered(out, ivs):
    """Set to '-' every position covered by >= 2 of the given half-open intervals."""
    ev = []
    for s, e in ivs:
        ev.append((s, 1))
        ev.append((e, -1))
    ev.sort(key=lambda x: (x[0], x[1]))  # ends (-1) before starts at the same position
    depth, prev, masked = 0, 0, 0
    for pos, d in ev:
        if depth >= 2 and pos > prev:
            out[prev:pos] = b"-" * (pos - prev)
            masked += pos - prev
        depth += d
        prev = pos
    return masked


def walk_cs(out, refseq, ts, te, cs):
    """Write the query's bases into `out` at reference coordinates, following a short cs string.
    Bases are on the reference strand (minimap2 reverse-complements the query for '-' hits)."""
    pos = ts
    for m in CS_RE.finditer(cs):
        tok = m.group(0)
        op = tok[0:1]
        if op == b":":
            n = int(tok[1:])
            out[pos:pos + n] = refseq[pos:pos + n]
            pos += n
        elif op == b"*":
            out[pos] = tok[2]
            pos += 1
        elif op == b"-":
            n = len(tok) - 1
            out[pos:pos + n] = b"-" * n
            pos += n
        elif op == b"+":
            pass  # insertion in the query: no reference column for it
        elif op == b"~":
            pos += int(re.search(rb"[0-9]+", tok).group(0))
    return pos == te


def cmd_project(a):
    chroms = read_list(a.chroms)
    keep = set(chroms)
    sys.stderr.write("[project] %s: loading reference\n" % a.name)
    ref = {k: v.upper() for k, v in read_fasta(a.ref, keep).items()}
    missing = [c for c in chroms if c not in ref]
    if missing:
        sys.exit("ERROR: chromosomes not in reference: %s" % ",".join(missing))

    alns = defaultdict(list)
    n_tot = n_kept = n_bad = 0
    if not a.self:
        with open(a.paf, "rb") as fh:
            for line in fh:
                f = line.rstrip(b"\n").split(b"\t")
                n_tot += 1
                t = f[5].decode()
                if t not in keep:
                    continue
                if int(f[11]) < a.min_mapq or int(f[10]) < a.min_aln:
                    continue
                tags = {x[:2]: x[5:] for x in f[12:]}
                if tags.get(b"tp", b"P") != b"P" or b"cs" not in tags:
                    continue
                alns[t].append((int(f[7]), int(f[8]), tags[b"cs"]))
                n_kept += 1

    os.makedirs(os.path.join(a.outdir, "coverage"), exist_ok=True)
    cov_rows = []
    for c in chroms:
        refseq = ref[c]
        L = len(refseq)
        if a.self:
            out = bytearray(refseq)
            masked = 0
        else:
            out = bytearray(b"-" * L)
            for ts, te, cs in alns[c]:
                if not walk_cs(out, refseq, ts, te, cs):
                    n_bad += 1
            masked = mask_multi_covered(out, [(s, e) for s, e, _ in alns[c]])
        out = out.translate(TRANSLATE)
        covered = L - out.count(b"-")
        d = os.path.join(a.outdir, c)
        os.makedirs(d, exist_ok=True)
        tmp = os.path.join(d, a.name + ".fa.part")
        with open(tmp, "wb") as fh:
            fh.write(b">" + a.name.encode() + b"\n")
            for i in range(0, L, 80):
                fh.write(out[i:i + 80] + b"\n")
        os.replace(tmp, os.path.join(d, a.name + ".fa"))
        cov_rows.append((c, L, covered, masked, len(alns[c])))
    with open(os.path.join(a.outdir, "coverage", a.name + ".tsv"), "w") as fh:
        fh.write("chrom\tlength\tcovered_bp\tfrac_covered\tmulti_covered_masked_bp\tn_alignments\n")
        for c, L, cv, mk, n in cov_rows:
            fh.write("%s\t%d\t%d\t%.4f\t%d\t%d\n" % (c, L, cv, cv / L if L else 0, mk, n))
    tot_L = sum(r[1] for r in cov_rows)
    tot_c = sum(r[2] for r in cov_rows)
    sys.stderr.write("[project] %s: %d/%d PAF lines kept, %d cs walks ended off-coordinate; "
                     "%.1f%% of reference chromosomes covered\n"
                     % (a.name, n_kept, n_tot, n_bad, 100.0 * tot_c / max(tot_L, 1)))
    if n_bad:
        sys.exit("ERROR: %d alignments' cs strings did not reach their end coordinate "
                 "(was minimap2 run with -c --cs?)" % n_bad)


# ------------------------------------------------------------------------- msa --
def cmd_msa(a):
    chroms = read_list(a.chroms)
    names = read_list(a.names)
    os.makedirs(a.outdir, exist_ok=True)
    cov = {}
    for n in names:
        p = os.path.join(a.projdir, "coverage", n + ".tsv")
        with open(p) as fh:
            next(fh)
            for line in fh:
                f = line.split("\t")
                cov[(n, f[0])] = float(f[3])
    log = open(a.log, "w")
    log.write("chrom\tsample\tfrac_covered\tincluded\n")
    for c in chroms:
        out = os.path.join(a.outdir, c + ".fasta")
        tmp = out + ".part"
        n_in = 0
        with open(tmp, "wb") as oh:
            for n in names:
                fc = cov[(n, c)]
                inc = (n in a.always.split(",")) or fc >= a.min_cov
                log.write("%s\t%s\t%.4f\t%s\n" % (c, n, fc, "yes" if inc else "no"))
                if inc:
                    with open(os.path.join(a.projdir, c, n + ".fa"), "rb") as ih:
                        shutil.copyfileobj(ih, oh, 16 << 20)
                    n_in += 1
        os.replace(tmp, out)
        sys.stderr.write("[msa] %s: %d of %d sequences\n" % (c, n_in, len(names)))
    log.close()


# ----------------------------------------------------------------------- annot --
def cmd_annot(a):
    chroms = read_list(a.chroms)
    keep = set(chroms)
    lens = read_fai(a.fai)
    with open(a.out_bed, "w") as fh:
        fh.write("Chromosome\tStart\tEnd\n")
        for c in chroms:
            fh.write("%s\t0\t%d\n" % (thex_name(c), lens[c]))
    n = 0
    with open(a.gff) as ih, open(a.out_gff, "w") as oh:
        for line in ih:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != a.feature or f[0] not in keep:
                continue
            f[0] = thex_name(f[0])
            oh.write("\t".join(f) + "\n")
            n += 1
    sys.stderr.write("[annot] %d %s features written\n" % (n, a.feature))


# ------------------------------------------------------------------- summarize --
def window_files(windowed_dir):
    """(thex_chrom, start, end, stem) for every window THExBuilder's windower wrote."""
    out = []
    for d in sorted(glob.glob(os.path.join(windowed_dir, "*"))):
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if not f.endswith(".fasta"):
                continue
            stem = f[:-len(".fasta")]
            p = stem.split("_")
            out.append((p[0], int(p[1]), int(p[2]), stem, os.path.basename(d)))
    return out


def alignment_stats(path):
    """(n_taxa, aligned_length) of a FASTA alignment."""
    n, L, first = 0, 0, True
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
                if n == 2:
                    first = False
            elif first:
                L += len(line.strip())
    return n, L


def classify(tree, sp_of, outgroup, species):
    """Species-level topology of one rooted tree, plus the support of the defining clade."""
    from itertools import combinations
    leaves = tree.get_leaf_names()
    present = {s: [l for l in leaves if sp_of.get(l) == s] for s in species}
    if any(len(v) == 0 for v in present.values()):
        return "Species missing", None
    mono = {s: tree.check_monophyly(values=present[s], target_attr="name")[0] for s in species}
    for A, B in combinations(species, 2):
        C = [s for s in species if s not in (A, B)][0]
        ab = present[A] + present[B]
        if tree.check_monophyly(values=ab, target_attr="name")[0]:
            node = tree.get_common_ancestor(ab)
            sup = node.support if node is not tree else None
            label = "((%s,%s),%s)" % (A, B, C)
            if not mono[C]:
                label += " %s paraphyletic" % C
            elif not (mono[A] and mono[B]):
                label += " %s/%s unsorted" % (A, B)
            return label, sup
    return "No species pair forms a clade", None


def cmd_summarize(a):
    import pandas as pd
    from ete3 import Tree

    species = a.species.split(",")
    names = read_list(a.names)
    sp_of = {n: n.split("_")[0] for n in names if n.split("_")[0] in species}

    # genes per window, from the THEx GFF (already in THEx chromosome names)
    starts, ends = defaultdict(list), defaultdict(list)
    if a.gff and os.path.exists(a.gff):
        with open(a.gff) as fh:
            for line in fh:
                f = line.split("\t")
                starts[f[0]].append(int(f[3]))
                ends[f[0]].append(int(f[4]))
        for c in starts:
            starts[c].sort()
            ends[c].sort()

    wins = window_files(a.windowed_dir)
    rows, status = [], []
    for chrom, ws, we, stem, cdir in wins:
        pw = os.path.join(a.pw_dir, cdir, stem + ".fasta")
        tr = os.path.join(a.trimal_dir, cdir, stem + ".fasta")
        ct = os.path.join(a.iqtree_dir, cdir, stem + ".contree")
        genes = (bisect.bisect_right(starts[chrom], we) - bisect.bisect_left(ends[chrom], ws)) if chrom in starts else 0
        rec = {"Chromosome": chrom, "Window": we, "NewickTree": "NoTree", "TopologyID": None,
               "CladeSupport": None, "Taxa": None, "AlignedSites": None, "Genes": genes}
        if os.path.exists(ct) and os.path.getsize(ct) > 0:
            t = Tree(open(ct).readline().strip())
            if a.outgroup in t.get_leaf_names():
                t.set_outgroup(t & a.outgroup)
                t.ladderize()
                rec["TopologyID"], rec["CladeSupport"] = classify(t, sp_of, a.outgroup, species)
                st = "tree"
            else:
                rec["TopologyID"] = "Outgroup missing"
                st = "tree_no_outgroup"
            rec["NewickTree"] = t.write(format=0)
            rec["Taxa"], rec["AlignedSites"] = alignment_stats(pw)
        elif os.path.exists(pw):
            st = "tree_missing"
        elif os.path.exists(os.path.join(a.pw_dir, cdir, stem + "-DROPPED.fasta")):
            st = "dropped_pairwise_filter" if os.path.exists(tr) else "dropped_trimal"
        else:
            st = "dropped_trimal"
        rows.append(rec)
        status.append((chrom, ws, we, st))

    df = pd.DataFrame(rows, columns=["Chromosome", "Window", "NewickTree", "TopologyID",
                                     "CladeSupport", "Taxa", "AlignedSites", "Genes"])
    order = {thex_name(c): i for i, c in enumerate(read_list(a.chroms))}
    df["_o"] = df["Chromosome"].map(order)
    df = df.sort_values(["_o", "Window"]).drop(columns="_o")
    df.to_excel(a.out_xlsx, index=False)

    st = pd.DataFrame(status, columns=["Chromosome", "Start", "End", "Status"])
    st.to_csv(a.out_prefix + ".window_status.tsv", sep="\t", index=False)
    trees = df[df["NewickTree"] != "NoTree"]
    is_z = trees["Chromosome"] == "chr-Z"
    summ = pd.DataFrame({
        "all_windows": trees["TopologyID"].value_counts(),
        "autosomes": trees.loc[~is_z, "TopologyID"].value_counts(),
        "chr_Z": trees.loc[is_z, "TopologyID"].value_counts(),
    }).fillna(0).astype(int)
    summ["frac_all"] = (summ["all_windows"] / summ["all_windows"].sum()).round(4)
    summ.sort_values("all_windows", ascending=False).to_csv(
        a.out_prefix + ".topology_summary.tsv", sep="\t", index_label="TopologyID")
    per_chrom = trees.groupby(["Chromosome", "TopologyID"]).size().unstack(fill_value=0)
    per_chrom = per_chrom.loc[sorted(per_chrom.index, key=lambda c: order.get(c, 1e9))]
    per_chrom.to_csv(a.out_prefix + ".topology_by_chromosome.tsv", sep="\t")
    sys.stderr.write("[summarize] %d windows, %d with trees\n%s\n%s\n"
                     % (len(df), len(trees), st["Status"].value_counts().to_string(),
                        summ.to_string()))


# ------------------------------------------------------------------------ main --
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd")  # (required= needs Python 3.7; the login node's python3 may be 3.6)

    s = sp.add_parser("pick-ref")
    s.add_argument("--ids", required=True)
    s.add_argument("--fai-dir", required=True)
    s.add_argument("--exclude", default="chr_W,chr_MT")
    s.add_argument("--out-table", required=True)

    s = sp.add_parser("outgroup-ref")
    s.add_argument("--fasta", required=True)
    s.add_argument("--gff")
    s.add_argument("--chr-map", required=True)
    s.add_argument("--out-fasta", required=True)
    s.add_argument("--out-gff")

    s = sp.add_parser("project")
    s.add_argument("--ref", required=True)
    s.add_argument("--paf")
    s.add_argument("--self", action="store_true", help="the query IS the reference")
    s.add_argument("--name", required=True)
    s.add_argument("--chroms", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--min-mapq", type=int, default=5)
    s.add_argument("--min-aln", type=int, default=5000)

    s = sp.add_parser("msa")
    s.add_argument("--projdir", required=True)
    s.add_argument("--chroms", required=True)
    s.add_argument("--names", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--min-cov", type=float, default=0.5)
    s.add_argument("--always", default="", help="comma-separated names included regardless of coverage")
    s.add_argument("--log", required=True, help="per chromosome x sample coverage table (keep it out of --outdir)")

    s = sp.add_parser("annot")
    s.add_argument("--gff", required=True)
    s.add_argument("--fai", required=True)
    s.add_argument("--chroms", required=True)
    s.add_argument("--feature", default="gene")
    s.add_argument("--out-gff", required=True)
    s.add_argument("--out-bed", required=True)

    s = sp.add_parser("summarize")
    s.add_argument("--windowed-dir", required=True)
    s.add_argument("--trimal-dir", required=True)
    s.add_argument("--pw-dir", required=True)
    s.add_argument("--iqtree-dir", required=True)
    s.add_argument("--names", required=True)
    s.add_argument("--chroms", required=True)
    s.add_argument("--gff")
    s.add_argument("--outgroup", required=True)
    s.add_argument("--species", default="STGR,LEPC,GRPC")
    s.add_argument("--out-xlsx", required=True)
    s.add_argument("--out-prefix", required=True)

    a = p.parse_args()
    if not a.cmd:
        p.error("choose a subcommand")
    {"pick-ref": cmd_pick_ref, "outgroup-ref": cmd_outgroup_ref, "project": cmd_project, "msa": cmd_msa,
     "annot": cmd_annot, "summarize": cmd_summarize}[a.cmd](a)


if __name__ == "__main__":
    main()
