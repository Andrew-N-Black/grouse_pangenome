#!/usr/bin/env python3
"""
Break reference-imposed RagTag joins that the Hi-C does not support, and write
pangenome-input FASTAs.

WHY
    RagTag ordered each haplotype's yahs scaffolds against the rock ptarmigan
    reference and joined them with 100 bp gaps. A pangenome builder (PGGB,
    Minigraph-Cactus) assigns sequences to chromosomes by alignment and gains
    nothing from that order, but an incorrect join becomes a false structural
    variant in the graph. So every join whose Hi-C support was judged
    'unsupported' by the assembly repo's step 13 (13_join_support_ptarmigan.sh)
    is cut at its gap,
    and the gap is removed. Nothing else changes: no bases are added, removed
    or reordered apart from the N runs of the broken gaps.

INPUTS (per haplotype, prefix <SP>_<SAMPLE>_<HAP>)
    <final-dir>/<P>.pseudo_chr.fasta (+ .fai)     current assembly
    <ragtag-dir>/<S>.<H>/ragtag.scaffold.agp       where every join is
    <ragtag-dir>/<S>.<H>/<S>.<H>.rename_map.tsv    AGP object -> current name
    <joins-dir>/<P>.joins.tsv                       step 13 verdicts
    <final-dir>/W_relabel_log.tsv (optional)        male chr_W -> scaffold_N

    The chromosome names were changed twice after RagTag (homology names, then
    native names in step 14; males' chr_W later relabelled). AGP coordinates are
    still valid because renaming never changed a base, and this script PROVES
    that per chromosome before cutting: the AGP object length must equal the
    FASTA length, and every gap it cuts must be all-N in the FASTA. Any mismatch
    stops the run.

    The position column of the joins table is not assumed: the script finds the
    numeric column whose values fall on AGP gaps, and refuses to proceed if no
    column does so for at least 95% of rows.

OUTPUTS (in --out-dir)
    <P>.pg_input.fasta    assembly with broken chromosomes replaced by pieces
                          <chrom>_p1, <chrom>_p2 ... in position order
    breaks.tsv            every join cut: where, why, which pieces it made
    pieces.tsv            every piece: parent chromosome and coordinates
    summary.tsv           per haplotype: sequences/bp in and out
    The unplaced_short FASTA from step 09 is not touched or copied.

JOINS ON THE REFERENCE'S UNPLACED SCAFFOLDS
    RagTag also concatenated grouse scaffolds that align to the same UNPLACED
    ptarmigan scaffold (NW_ accessions); step 09 named those objects
    scaffold_N. Step 13 scored joins on chr_* sequences only, so these joins
    were never tested — and they rest on the weakest evidence RagTag has (an
    unplaced, often repetitive reference scaffold). By default they are all
    cut, verdict 'untested_unplaced'. --keep-unplaced-joins leaves them.

MODES (--break)
    unsupported              default: cut joins step 13 called unsupported
    unsupported,ambiguous    also cut joins it could not call either way
    all                      cut EVERY RagTag join, i.e. return each chromosome
                             to its yahs scaffolds (no joins table needed)

USAGE
    python3 01_break_unsupported_joins.py --final-dir ... --ragtag-dir ... \
        --joins-dir ... --out-dir ... [--dry-run] [--only PREFIX]
"""

import argparse
import glob
import os
import re
import sys
from collections import OrderedDict, defaultdict

KNOWN_VERDICTS = {"supported", "unsupported", "ambiguous",
                  "not_scored_small", "not_scored", "uncalibrated"}
POS_MATCH_MIN = 0.95      # fraction of joins rows that must land on an AGP gap
POS_TOL = 1               # bp slack either side of a gap (0/1-based conventions)


# --------------------------------------------------------------------- readers
def read_fai(path):
    out = OrderedDict()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


def read_agp(path):
    """object -> list of rows; each row a dict. Gaps are type N or U."""
    objs = OrderedDict()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            row = {"beg": int(f[1]), "end": int(f[2]), "type": f[4]}
            if f[4] in ("N", "U"):
                row["gap"] = True
            else:
                row["gap"] = False
                row["comp"] = f[5]
            objs.setdefault(f[0], []).append(row)
    return objs


def read_two_col(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[1]:
                out[f[0]] = f[1]
    return out


def read_w_log(path):
    """haplotype prefix -> new name given to its chr_W (males only)."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1] == "chr_W":
                out[f[0]] = f[2]
    return out


def read_joins(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    return header, rows


def fasta_records(path):
    """Yield (name, header_line, sequence_string)."""
    name, head, chunks = None, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, head, "".join(chunks)
                head = line.rstrip("\n")
                name = head[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, head, "".join(chunks)


def write_seq(fh, name, seq, width=80):
    fh.write(">" + name + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


# ------------------------------------------------------------------ resolving
def parse_prefix(prefix):
    m = re.fullmatch(r"([A-Za-z]+)_([A-Za-z0-9.]+)_(hap\d+)", prefix)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3)


def resolve_chromosomes(prefix, fai, agp, rmap, w_new):
    """
    current FASTA name -> AGP object, for chromosomes RagTag built.

    Primary route: rename map (AGP object -> current name), with the male W
    relabel applied on top. Fallback for any chr_* left unmatched: a unique
    AGP object of exactly the same length. Every match is length-checked.
    """
    name_to_obj, notes = {}, []
    for obj, rows in agp.items():
        if not obj.endswith("_RagTag"):
            continue
        name = rmap.get(obj)
        if name is None:
            continue
        if name == "chr_W" and w_new:
            name = w_new
        if name not in fai:
            continue          # e.g. a placement that went to unplaced_short
        obj_len = max(r["end"] for r in rows)
        if fai[name] != obj_len:
            sys.exit("ERROR: %s: %s is %d bp in the FASTA but its AGP object %s is "
                     "%d bp. AGP coordinates do not transfer; stopping."
                     % (prefix, name, fai[name], obj, obj_len))
        name_to_obj[name] = obj

    used = set(name_to_obj.values())
    for name, L in fai.items():
        if not name.startswith("chr_") or name in name_to_obj:
            continue
        cands = [o for o, rows in agp.items()
                 if o.endswith("_RagTag") and o not in used
                 and max(r["end"] for r in rows) == L]
        if len(cands) == 1:
            name_to_obj[name] = cands[0]
            used.add(cands[0])
            notes.append("%s matched to %s by length" % (name, cands[0]))
    return name_to_obj, notes


def gaps_of(rows):
    """[(row_index, beg, end)] for gap rows."""
    return [(i, r["beg"], r["end"]) for i, r in enumerate(rows) if r["gap"]]


def match_gap(gaps, pos):
    hits = [g for g in gaps if g[1] - POS_TOL <= pos <= g[2] + POS_TOL]
    if not hits:
        return None
    return min(hits, key=lambda g: abs((g[1] + g[2]) / 2.0 - pos))


# ---------------------------------------------------------- column detection
def detect_columns(tables, chrom_gaps):
    """
    tables: prefix -> (header, rows). chrom_gaps: prefix -> {chrom: gaps}.
    Returns (chrom_col, verdict_col, pos_col, header) or exits with a report.
    """
    header = None
    all_rows = []
    for prefix, (h, rows) in tables.items():
        if header is None:
            header = h
        elif h != header:
            sys.exit("ERROR: joins tables do not share one header (%s differs)" % prefix)
        all_rows.extend((prefix, r) for r in rows)
    if not all_rows:
        sys.exit("ERROR: the joins tables contain no rows")
    ncol = max(len(r) for _, r in all_rows)

    def frac(pred, col):
        vals = [r[col] for _, r in all_rows if len(r) > col]
        return sum(1 for v in vals if pred(v)) / float(len(vals)) if vals else 0.0

    verdict_col = max(range(ncol), key=lambda c: frac(lambda v: v in KNOWN_VERDICTS, c))
    if frac(lambda v: v in KNOWN_VERDICTS, verdict_col) < 0.9:
        sys.exit("ERROR: no column of the joins table holds verdicts. Header: %s" % header)
    chrom_col = max(range(ncol), key=lambda c: frac(lambda v: v.startswith("chr_"), c))
    if frac(lambda v: v.startswith("chr_"), chrom_col) < 0.9:
        sys.exit("ERROR: no column of the joins table holds chromosome names. Header: %s"
                 % header)

    scores = {}
    for c in range(ncol):
        if c in (verdict_col, chrom_col):
            continue
        hit = tot = 0
        for prefix, r in all_rows:
            if len(r) <= c or not re.fullmatch(r"-?\d+", r[c] or ""):
                continue
            gaps = chrom_gaps.get(prefix, {}).get(r[chrom_col])
            if gaps is None:
                continue
            tot += 1
            if match_gap(gaps, int(r[c])) is not None:
                hit += 1
        if tot:
            scores[c] = (hit / float(tot), tot)
    if not scores:
        sys.exit("ERROR: no integer column in the joins table to locate joins by")
    pos_col = max(scores, key=lambda c: (scores[c][0], scores[c][1]))
    print(">>> joins table columns (0-based):")
    print("    chromosome : %d (%s)" % (chrom_col, header[chrom_col]))
    print("    verdict    : %d (%s)" % (verdict_col, header[verdict_col]))
    for c in sorted(scores):
        flag = "  <- position" if c == pos_col else ""
        print("    col %-3d %-22s lands on an AGP gap in %5.1f%% of %d rows%s"
              % (c, header[c][:22], 100 * scores[c][0], scores[c][1], flag))
    if scores[pos_col][0] < POS_MATCH_MIN:
        sys.exit("ERROR: no column lands on AGP gaps for >= %.0f%% of rows; the "
                 "joins table and the AGPs do not describe the same assemblies."
                 % (100 * POS_MATCH_MIN))
    return chrom_col, verdict_col, pos_col, header


# ------------------------------------------------------------------ cutting
def plan_pieces(rows, break_idx):
    """
    Cut an AGP object at the gap rows in break_idx. Returns [(beg, end)] in
    object coordinates (1-based, inclusive). Leading/trailing gaps of a piece
    are dropped; unbroken interior gaps are kept as they are.
    """
    pieces, start, end = [], None, None
    for i, r in enumerate(rows):
        if r["gap"]:
            if i in break_idx and start is not None:
                pieces.append((start, end))
                start = end = None
            continue
        if start is None:
            start = r["beg"]
        end = r["end"]
    if start is not None:
        pieces.append((start, end))
    return pieces


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--ragtag-dir", required=True)
    ap.add_argument("--joins-dir", help="step 13 output (not needed with --break all)")
    ap.add_argument("--w-log", help="W_relabel_log.tsv (default: <final-dir>/W_relabel_log.tsv)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--break", dest="brk", default="unsupported",
                    help="comma list of verdicts to cut, or 'all' (default: unsupported)")
    ap.add_argument("--only", default=None, help="process only this prefix, e.g. GRPC_F5595_hap1")
    ap.add_argument("--keep-unplaced-joins", action="store_true",
                    help="do NOT cut joins inside scaffold_N objects (RagTag placements on the "
                         "reference's unplaced scaffolds), which step 13 never scored")
    ap.add_argument("--dry-run", action="store_true", help="report the plan, write nothing")
    args = ap.parse_args()

    cut_all = args.brk.strip() == "all"
    cut_verdicts = set() if cut_all else {v.strip() for v in args.brk.split(",") if v.strip()}
    bad = cut_verdicts - KNOWN_VERDICTS
    if bad:
        sys.exit("ERROR: unknown verdict(s) in --break: %s" % ", ".join(sorted(bad)))
    if not cut_all and not args.joins_dir:
        sys.exit("ERROR: --joins-dir is required unless --break all")

    w_log = args.w_log or os.path.join(args.final_dir, "W_relabel_log.tsv")
    w_new = read_w_log(w_log)

    fastas = sorted(glob.glob(os.path.join(args.final_dir, "*_hap[0-9].pseudo_chr.fasta")))
    if args.only:
        fastas = [f for f in fastas if os.path.basename(f).startswith(args.only + ".")]
    if not fastas:
        sys.exit("ERROR: no *_hapN.pseudo_chr.fasta in %s" % args.final_dir)

    # ---- pass 1: resolve every haplotype's chromosomes to AGP objects -------
    hap = OrderedDict()
    for fa in fastas:
        prefix = os.path.basename(fa)[:-len(".pseudo_chr.fasta")]
        parsed = parse_prefix(prefix)
        if not parsed:
            sys.exit("ERROR: cannot parse <SP>_<SAMPLE>_<HAP> from %s" % prefix)
        _sp, sample, h = parsed
        rdir = os.path.join(args.ragtag_dir, "%s.%s" % (sample, h))
        agp_p = os.path.join(rdir, "ragtag.scaffold.agp")
        map_p = os.path.join(rdir, "%s.%s.rename_map.tsv" % (sample, h))
        fai_p = fa + ".fai"
        for p in (agp_p, map_p, fai_p):
            if not os.path.isfile(p):
                sys.exit("ERROR: %s: missing %s" % (prefix, p))
        fai, agp, rmap = read_fai(fai_p), read_agp(agp_p), read_two_col(map_p)
        name_to_obj, notes = resolve_chromosomes(prefix, fai, agp, rmap, w_new.get(prefix))
        for n in notes:
            print("  note: %s: %s" % (prefix, n))
        unresolved = [n for n in fai if n.startswith("chr_") and n not in name_to_obj]
        if unresolved:
            # a chr_ name with no AGP object cannot be cut; it can only be
            # carried through untouched, which is safe, so warn rather than stop
            print("  warn: %s: no AGP object for %s (carried through uncut)"
                  % (prefix, ", ".join(unresolved)))
        chrom_gaps = {n: gaps_of(agp[o]) for n, o in name_to_obj.items()}
        hap[prefix] = {"fa": fa, "fai": fai, "agp": agp, "name_to_obj": name_to_obj,
                       "chrom_gaps": chrom_gaps}
    n_joins = sum(len(g) for h in hap.values() for g in h["chrom_gaps"].values())
    print(">>> %d haplotypes, %d RagTag joins in total" % (len(hap), n_joins))

    # ---- pass 2: which gaps to cut ------------------------------------------
    cuts = defaultdict(dict)        # prefix -> chrom -> {row_idx: (verdict, row)}
    joins_header = []
    if cut_all:
        for prefix, h in hap.items():
            for chrom, gaps in h["chrom_gaps"].items():
                for gi, b, e in gaps:
                    cuts[prefix].setdefault(chrom, {})[gi] = ("all", [])
    else:
        tables = OrderedDict()
        for prefix in hap:
            p = os.path.join(args.joins_dir, prefix + ".joins.tsv")
            if not os.path.isfile(p):
                sys.exit("ERROR: no joins table for %s at %s" % (prefix, p))
            header, rows = read_joins(p)
            # males: the joins table still says chr_W for what is now scaffold_N
            wn = w_new.get(prefix)
            tables[prefix] = (header, rows if not wn else
                              [[(wn if v == "chr_W" else v) for v in r] for r in rows])
        chrom_col, verdict_col, pos_col, joins_header = detect_columns(
            tables, {p: h["chrom_gaps"] for p, h in hap.items()})
        unmatched = 0
        for prefix, (_h, rows) in tables.items():
            for r in rows:
                if len(r) <= max(chrom_col, verdict_col, pos_col):
                    continue
                if r[verdict_col] not in cut_verdicts:
                    continue
                chrom = r[chrom_col]
                gaps = hap[prefix]["chrom_gaps"].get(chrom)
                g = match_gap(gaps, int(r[pos_col])) if gaps else None
                if g is None:
                    unmatched += 1
                    print("  ERROR: %s %s %s=%s (%s) is not on any AGP gap"
                          % (prefix, chrom, joins_header[pos_col], r[pos_col], r[verdict_col]))
                    continue
                cuts[prefix].setdefault(chrom, {})[g[0]] = (r[verdict_col], r)
        if unmatched:
            sys.exit("ERROR: %d join(s) selected for cutting could not be placed; "
                     "nothing written." % unmatched)

        # joins step 13 never saw: inside objects RagTag built on the
        # reference's unplaced scaffolds (named scaffold_N, not chr_*)
        seen = {(p, r[chrom_col]) for p, (_h, rows) in tables.items()
                for r in rows if len(r) > chrom_col}
        n_unpl = 0
        for prefix, h in hap.items():
            for chrom, gaps in h["chrom_gaps"].items():
                if chrom.startswith("chr_") or (prefix, chrom) in seen or not gaps:
                    continue
                if args.keep_unplaced_joins:
                    n_unpl += len(gaps)
                    continue
                for gi, _b, _e in gaps:
                    cuts[prefix].setdefault(chrom, {})[gi] = ("untested_unplaced", [])
                    n_unpl += 1
        print(">>> joins inside reference-unplaced placements (never scored): %d — %s"
              % (n_unpl, "kept" if args.keep_unplaced_joins else "cut"))
        n_unscored_chr = sum(1 for prefix, h in hap.items()
                             for chrom, gaps in h["chrom_gaps"].items()
                             if chrom.startswith("chr_") and (prefix, chrom) not in seen
                             for _ in gaps)
        if n_unscored_chr:
            print("  warn: %d join(s) on chr_* sequences have no row in the joins tables "
                  "(not cut)" % n_unscored_chr)

    # ---- pass 3: cut and write ----------------------------------------------
    if not args.dry_run:
        os.makedirs(args.out_dir, exist_ok=True)
    brk_fh = pcs_fh = sum_fh = None
    if not args.dry_run:
        brk_fh = open(os.path.join(args.out_dir, "breaks.tsv.part"), "w")
        pcs_fh = open(os.path.join(args.out_dir, "pieces.tsv.part"), "w")
        sum_fh = open(os.path.join(args.out_dir, "summary.tsv.part"), "w")
        brk_fh.write("haplotype\tchrom\tgap_beg\tgap_end\tverdict\tleft_piece\tright_piece\t"
                     + "\t".join("joins:" + c for c in joins_header) + "\n")
        pcs_fh.write("haplotype\tpiece\tparent\tparent_beg\tparent_end\tlength\n")
        sum_fh.write("haplotype\tseqs_in\tbp_in\tchroms_cut\tcuts\tseqs_out\tbp_out\tgap_bp_removed\n")

    tot_cuts = 0
    for prefix, h in hap.items():
        hcuts = cuts.get(prefix, {})
        ncut = sum(len(v) for v in hcuts.values())
        tot_cuts += ncut
        lab_chr = ", ".join("%s x%d" % (c, len(v)) for c, v in sorted(hcuts.items())
                            if c.startswith("chr_"))
        n_up = sum(len(v) for c, v in hcuts.items() if not c.startswith("chr_"))
        label = (lab_chr or "none") + ("  (+%d unplaced)" % n_up if n_up else "")
        print("  %-22s cuts: %s" % (prefix, label))
        if args.dry_run:
            continue

        out_fa = os.path.join(args.out_dir, prefix + ".pg_input.fasta")
        seqs_in = bp_in = seqs_out = bp_out = gap_removed = 0
        with open(out_fa + ".part", "w") as out:
            for name, _head, seq in fasta_records(h["fa"]):
                seqs_in += 1
                bp_in += len(seq)
                if len(seq) != h["fai"].get(name, -1):
                    sys.exit("ERROR: %s: %s length differs from its .fai; reindex first"
                             % (prefix, name))
                if name not in hcuts:
                    write_seq(out, name, seq)
                    seqs_out += 1
                    bp_out += len(seq)
                    continue
                rows = h["agp"][h["name_to_obj"][name]]
                bidx = set(hcuts[name])
                for i in bidx:
                    gs = seq[rows[i]["beg"] - 1:rows[i]["end"]]
                    if gs.strip("Nn"):
                        sys.exit("ERROR: %s %s: AGP gap %d-%d is not all N in the FASTA. "
                                 "Coordinates do not line up; stopping."
                                 % (prefix, name, rows[i]["beg"], rows[i]["end"]))
                pieces = plan_pieces(rows, bidx)
                pnames = ["%s_p%d" % (name, k + 1) for k in range(len(pieces))]
                kept = sum(e - b + 1 for b, e in pieces)
                cut_bp = sum(rows[i]["end"] - rows[i]["beg"] + 1 for i in bidx)
                edge_gap = len(seq) - kept - cut_bp
                if edge_gap < 0:
                    sys.exit("ERROR: %s %s: pieces overlap; stopping" % (prefix, name))
                for pn, (b, e) in zip(pnames, pieces):
                    write_seq(out, pn, seq[b - 1:e])
                    pcs_fh.write("%s\t%s\t%s\t%d\t%d\t%d\n" % (prefix, pn, name, b, e, e - b + 1))
                    seqs_out += 1
                    bp_out += e - b + 1
                gap_removed += len(seq) - kept
                # record each cut with the pieces either side of it
                for i in sorted(bidx):
                    gb, ge = rows[i]["beg"], rows[i]["end"]
                    left = [pn for pn, (b, e) in zip(pnames, pieces) if e < gb]
                    right = [pn for pn, (b, e) in zip(pnames, pieces) if b > ge]
                    verdict, jrow = hcuts[name][i]
                    brk_fh.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (
                        prefix, name, gb, ge, verdict, left[-1] if left else "-",
                        right[0] if right else "-", "\t".join(jrow)))
        if seqs_out == 0:
            sys.exit("ERROR: %s: wrote no sequences" % prefix)
        os.replace(out_fa + ".part", out_fa)
        sum_fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
            prefix, seqs_in, bp_in, len(hcuts), ncut, seqs_out, bp_out, gap_removed))

    if not args.dry_run:
        for fh, nm in ((brk_fh, "breaks.tsv"), (pcs_fh, "pieces.tsv"), (sum_fh, "summary.tsv")):
            fh.close()
            os.replace(os.path.join(args.out_dir, nm + ".part"), os.path.join(args.out_dir, nm))
    print(">>> %s: %d cut(s) across %d haplotype(s)%s"
          % ("plan" if args.dry_run else "done", tot_cuts, len(hap),
             " — dry run, nothing written" if args.dry_run else ""))


if __name__ == "__main__":
    main()
