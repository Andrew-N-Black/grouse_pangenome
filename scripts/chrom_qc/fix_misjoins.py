#!/usr/bin/env python3
"""
fix_misjoins.py - split chromosome-scale misjoins found by chrom_qc.py at the nearest
scaffolding gap, and relabel each piece by the anchor chromosome it aligns to.

For every sequence, alignments to the anchor are collapsed into blocks along the query.
A boundary between an expected-chromosome block and an off-target block is CUT when the
off-target block is large (>= --min-off) or moderately large (>= --min-off-len) on a
chromosome that is also a length outlier. Each cut is snapped to the nearest run of Ns
(a RagTag/scaffolding gap) within --window bp; the Ns are removed.
If no gap is close enough, the join is inside a contig (an assembly chimera): it is NOT
cut and is written to manual_review.tsv with coordinates to check against Hi-C.

Pieces that belong to the expected chromosome are re-joined in their original order
(100-N gap) under the original name; pieces that belong elsewhere are named <true_chrom>_pN (first unused N in that haplotype), so
downstream steps that group chr_10, chr_10_p1, chr_10_p2 together pick them up.

Outputs (in --outdir):
  <hap>.fixed.fa      corrected FASTA, only for haplotypes that changed
  cuts.tsv            every cut made
  manual_review.tsv   joins with no nearby gap, plus smaller off-target blocks to eyeball
  haps.fixed.txt      haps.txt with changed haplotypes pointing to the fixed FASTA

Usage:
  python3 fix_misjoins.py --haps haps.txt --pafdir chromqc_paf --outdir fixed
"""
import argparse, collections, gzip, os, re, statistics, sys

SPLIT = re.compile(r'_p\d+$')
EXT = re.compile(r'(\.fai)?(\.(fa|fasta|fna))?(\.gz)?(\.paf)?$')
NRUN = re.compile(r'[Nn]{10,}')

def base_chrom(name): return SPLIT.sub('', name.split('#')[-1])

def hap_id(path):
    b = os.path.basename(path); prev = None
    while prev != b:
        prev, b = b, EXT.sub('', b)
    return b

def opener(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def read_fai(fa):
    with open(fa + '.fai') as fh:
        return {l.split('\t')[0]: int(l.split('\t')[1]) for l in fh}

def read_paf(paf, min_mapq, min_len, prefix):
    alns = collections.defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 12 or 'tp:A:S' in c[12:] or int(c[11]) < min_mapq or int(c[10]) < min_len:
                continue
            t = base_chrom(c[5])
            if t.startswith(prefix):
                alns[c[0]].append((int(c[2]), int(c[3]), t))
    return alns

def blocks(alist, min_block):
    out = []
    for qs, qe, t in sorted(alist):
        if out and out[-1][2] == t:
            out[-1][1] = max(out[-1][1], qe); out[-1][3] += qe - qs
        else:
            out.append([qs, qe, t, qe - qs])
    out = [b for b in out if b[3] >= min_block]
    merged = []
    for b in out:
        if merged and merged[-1][2] == b[2]:
            merged[-1][1] = max(merged[-1][1], b[1]); merged[-1][3] += b[3]
        else:
            merged.append(list(b))
    return merged

def read_seqs(fa, wanted):
    seqs, name, buf = {}, None, []
    with opener(fa) as fh:
        for line in fh:
            if line.startswith('>'):
                if name in wanted: seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted: seqs[name] = ''.join(buf)
    return seqs

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--haps', required=True, help='haps.txt used for the QC run')
    ap.add_argument('--pafdir', default='chromqc_paf')
    ap.add_argument('--outdir', default='fixed')
    ap.add_argument('--chrom-prefix', default='chr')
    ap.add_argument('--min-mapq', type=int, default=20)
    ap.add_argument('--min-aln', type=int, default=10000)
    ap.add_argument('--min-block', type=int, default=500000, help='blocks smaller than this are ignored when segmenting')
    ap.add_argument('--min-off', type=int, default=2500000, help='off-target block size cut unconditionally')
    ap.add_argument('--min-off-len', type=int, default=1000000, help='off-target block size cut if chromosome is a length outlier')
    ap.add_argument('--len-tol', type=float, default=0.15)
    ap.add_argument('--window', type=int, default=300000, help='max distance from breakpoint to a gap')
    ap.add_argument('--dry-run', action='store_true', help='report cuts but do not write FASTA')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    fastas = [l.strip() for l in open(a.haps) if l.strip()]
    lens = {f: read_fai(f) for f in fastas}

    tot = collections.Counter()
    for f, d in lens.items():
        for n, l in d.items():
            c = base_chrom(n)
            if c.startswith(a.chrom_prefix): tot[(f, c)] += l
    per_c = collections.defaultdict(list)
    for (f, c), l in tot.items(): per_c[c].append(l)
    med = {c: statistics.median(v) for c, v in per_c.items()}

    cutlog = open(os.path.join(a.outdir, 'cuts.tsv'), 'w')
    cutlog.write('haplotype\tsequence\tbreak_est\tgap_start\tgap_end\tdist_to_gap\tleft_piece\tright_piece\n')
    review = open(os.path.join(a.outdir, 'manual_review.tsv'), 'w')
    review.write('haplotype\tsequence\tbreak_est\toff_target\toff_target_mb\treason\n')
    newlist = []
    n_cut = n_rev = 0

    for f in fastas:
        h = hap_id(f)
        paf = os.path.join(a.pafdir, h + '.paf')
        if not os.path.exists(paf):
            print(f'WARNING: no PAF for {h}, skipped', file=sys.stderr); newlist.append(f); continue
        alns = read_paf(paf, a.min_mapq, a.min_aln, a.chrom_prefix)
        plan = {}                                   # seqname -> list of (break_est, off_target, off_bp)
        for n, L in lens[f].items():
            c = base_chrom(n)
            if not c.startswith(a.chrom_prefix) or n not in alns: continue
            bl = blocks(alns[n], a.min_block)
            outlier = abs(tot[(f, c)] / med[c] - 1) > a.len_tol if med.get(c) else False
            for i in range(len(bl) - 1):
                l, r = bl[i], bl[i + 1]
                if l[2] == r[2]: continue
                off = r if l[2] == c else l if r[2] == c else (l if l[3] < r[3] else r)
                brk = (l[1] + r[0]) // 2
                if off[3] >= a.min_off or (off[3] >= a.min_off_len and outlier):
                    plan.setdefault(n, []).append((brk, off[2], off[3]))
                elif off[3] >= a.min_off_len:
                    review.write(f'{h}\t{n}\t{brk}\t{off[2]}\t{off[3]/1e6:.2f}\tsmall_off_target_block\n'); n_rev += 1
        if not plan:
            newlist.append(f); continue

        seqs = read_seqs(f, set(plan))
        used = set(lens[f])
        pieces_for = {}
        for n, brks in plan.items():
            s = seqs[n]
            gaps = [(m.start(), m.end()) for m in NRUN.finditer(s)]
            cuts = []
            for brk, off_t, off_bp in brks:
                best = min(gaps, key=lambda g: abs((g[0] + g[1]) // 2 - brk), default=None)
                if best and abs((best[0] + best[1]) // 2 - brk) <= a.window:
                    cuts.append((best, brk))
                else:
                    review.write(f'{h}\t{n}\t{brk}\t{off_t}\t{off_bp/1e6:.2f}\tno_gap_within_window(contig_chimera?)\n'); n_rev += 1
            if not cuts: continue
            cuts.sort()
            bounds, start = [], 0
            for (gs, ge), _ in cuts:
                bounds.append((start, gs)); start = ge
            bounds.append((start, len(s)))
            # label each piece by majority anchor chromosome
            pieces = []
            for ps, pe in bounds:
                cnt = collections.Counter()
                for qs, qe, t in alns[n]:
                    ov = min(qe, pe) - max(qs, ps)
                    if ov > 0: cnt[t] += ov
                t = cnt.most_common(1)[0][0] if cnt else base_chrom(n)
                pieces.append([ps, pe, t])
            # expected-chromosome pieces are re-joined (in original order, 100-N gap) under the
            # original name, so an off-target insertion in the middle is lifted out cleanly
            names, groups = [], collections.OrderedDict()
            for ps, pe, t in pieces:
                if t == base_chrom(n):
                    nm = n
                else:
                    k = 1
                    while f'{t}_p{k}' in used: k += 1
                    nm = f'{t}_p{k}'; used.add(nm)
                names.append(nm)
                groups.setdefault(nm, []).append((ps, pe))
            pieces_for[n] = list(groups.items())
            for i, ((gs, ge), brk) in enumerate(cuts):
                cutlog.write(f'{h}\t{n}\t{brk}\t{gs}\t{ge}\t{abs((gs+ge)//2-brk)}\t{names[i]}\t{names[i+1]}\n'); n_cut += 1

        if not pieces_for or a.dry_run:
            newlist.append(f); continue
        out_fa = os.path.join(a.outdir, h + '.fixed.fa')
        with opener(f) as fh, open(out_fa, 'w') as out:
            name, keep = None, True
            def flush(nm, buf):
                if nm in pieces_for:
                    s = ''.join(buf)
                    for pn, ranges in pieces_for[nm]:
                        sub = ('N' * 100).join(s[ps:pe].strip('Nn') for ps, pe in ranges)
                        out.write(f'>{pn}\n')
                        for i in range(0, len(sub), 80): out.write(sub[i:i+80] + '\n')
            buf = []
            for line in fh:
                if line.startswith('>'):
                    flush(name, buf); buf = []
                    name = line[1:].split()[0]; keep = name not in pieces_for
                    if keep: out.write(line)
                elif keep:
                    out.write(line)
                else:
                    buf.append(line.strip())
            flush(name, buf)
        newlist.append(os.path.abspath(out_fa))

    with open(os.path.join(a.outdir, 'haps.fixed.txt'), 'w') as fh:
        fh.write('\n'.join(newlist) + '\n')
    cutlog.close(); review.close()
    print(f'{n_cut} cuts made; {n_rev} joins for manual review -> {a.outdir}/', file=sys.stderr)

if __name__ == '__main__':
    main()
