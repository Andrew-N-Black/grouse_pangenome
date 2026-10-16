#!/usr/bin/env python3
"""
Helper for 05d_thex_rescue.sh.

  prune  For every window that got no tree in step 05, write a copy without the haplotypes
         that barely cover it (< --min-cov aligned), so the window can be retried with the
         remaining taxa. The outgroup and reference are never removed; a window is not
         retried if more than --max-drop haplotypes would go or if any species would be
         left with fewer than --min-per-species haplotypes.
  merge  Move the retried windows that passed trimal + the pairwise filter into step 05's
         pairwise_filtered_windows/ (replacing the -DROPPED placeholder), so step 05's
         iqtree and summarize stages pick them up unchanged.
"""
import argparse
import os
import shutil
import sys


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def cmd_prune(a):
    keep_always = set(x for x in a.keep.split(",") if x)
    species = a.species.split(",")
    tried = skipped = 0
    log = open(a.log, "w")
    log.write("chrom_dir\twindow\tstatus\tn_taxa_before\tn_taxa_after\tdropped\n")
    for cdir in sorted(os.listdir(a.windowed)):
        wdir = os.path.join(a.windowed, cdir)
        if not os.path.isdir(wdir):
            continue
        for f in sorted(os.listdir(wdir)):
            if not f.endswith(".fasta"):
                continue
            stem = f[:-6]
            if os.path.exists(os.path.join(a.iqtree, cdir, stem + ".contree")):
                continue                                   # already has a tree
            if not os.path.exists(os.path.join(a.pw, cdir, stem + "-DROPPED.fasta")):
                continue                                   # not a dropped window
            recs = read_fasta(os.path.join(wdir, f))
            if not recs or len(recs[0][1]) < a.min_len:
                continue                                   # chromosome-end stub
            drop = []
            for n, s in recs:
                cov = 1 - (s.count("-") + s.count("N")) / float(len(s))
                if cov < a.min_cov and n not in keep_always:
                    drop.append((n, cov))
            kept = [n for n, _ in recs if n not in {d for d, _ in drop}]
            per_sp = {sp: sum(1 for n in kept if n.split("_")[0] == sp) for sp in species}
            reason = None
            if not drop:
                reason = "nothing_to_drop"
            elif len(drop) > a.max_drop:
                reason = "too_many_missing"
            elif min(per_sp.values()) < a.min_per_species:
                reason = "species_too_thin"
            dropped_txt = ",".join("%s:%.2f" % d for d in drop)
            if reason:
                skipped += 1
                log.write("%s\t%s\t%s\t%d\t%d\t%s\n" % (cdir, stem, reason, len(recs), len(kept), dropped_txt))
                continue
            od = os.path.join(a.out, cdir)
            os.makedirs(od, exist_ok=True)
            with open(os.path.join(od, f), "w") as oh:
                for n, s in recs:
                    if n in kept:
                        oh.write(">%s\n" % n)
                        for i in range(0, len(s), 80):
                            oh.write(s[i:i + 80] + "\n")
            tried += 1
            log.write("%s\t%s\tretry\t%d\t%d\t%s\n" % (cdir, stem, len(recs), len(kept), dropped_txt))
    log.close()
    sys.stderr.write("[prune] %d windows retried with fewer haplotypes, %d not retried (see %s)\n"
                     % (tried, skipped, a.log))


def cmd_merge(a):
    n = 0
    for cdir in sorted(os.listdir(a.rescued_pw)):
        d = os.path.join(a.rescued_pw, cdir)
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if not f.endswith(".fasta") or "-DROPPED" in f:
                continue
            stem = f[:-6]
            dst = os.path.join(a.pw, cdir)
            placeholder = os.path.join(dst, stem + "-DROPPED.fasta")
            if os.path.exists(placeholder):
                os.remove(placeholder)
            shutil.copy(os.path.join(d, f), os.path.join(dst, f))
            n += 1
    sys.stderr.write("[merge] %d rescued windows added to %s\n" % (n, a.pw))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd")
    s = sp.add_parser("prune")
    s.add_argument("--windowed", required=True)
    s.add_argument("--pw", required=True)
    s.add_argument("--iqtree", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--log", required=True)
    s.add_argument("--keep", default="")
    s.add_argument("--species", default="STGR,LEPC,GRPC")
    s.add_argument("--min-cov", type=float, default=0.5)
    s.add_argument("--max-drop", type=int, default=4)
    s.add_argument("--min-per-species", type=int, default=4)
    s.add_argument("--min-len", type=int, default=10000)
    s = sp.add_parser("merge")
    s.add_argument("--rescued-pw", required=True)
    s.add_argument("--pw", required=True)
    a = p.parse_args()
    if not a.cmd:
        p.error("choose a subcommand")
    {"prune": cmd_prune, "merge": cmd_merge}[a.cmd](a)


if __name__ == "__main__":
    main()
