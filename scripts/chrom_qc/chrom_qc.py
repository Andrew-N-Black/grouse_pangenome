#!/usr/bin/env python3
"""
chrom_qc.py - flag mislabelled, chimeric and length-outlier chromosomes in a set of
RagTag-scaffolded haplotype assemblies BEFORE building per-chromosome pangenome graphs.

Inputs
  --paf   one PAF per haplotype: minimap2 -x asm10 --secondary=no ANCHOR.fa HAP.fa
  --fai   samtools faidx .fai for every haplotype FASTA (for the length screen)

Sequence names are PanSN (SPECIES.ID#hap#chr_N); the chromosome is the field after the
last '#', so query and anchor names are compared on that field only.

Flags
  LEN_OUTLIER   length differs from the cross-haplotype median by > --len-tol
  MISLABEL      most aligned bases go to a different anchor chromosome
  CHIMERA       >=2 blocks (each >= --min-block) align to different anchor chromosomes
  LOW_ALIGN     < --min-aln-frac of the sequence aligns at all
  MOSTLY_REV    most bases aligned to the expected chromosome are on the reverse strand
"""
import argparse, collections, statistics, sys

def chrom_of(name):
    return name.split('#')[-1]

def hap_of(name):
    return '#'.join(name.split('#')[:-1]) if '#' in name else name

def read_fai(files, chrom_prefix):
    lengths = {}
    for f in files:
        with open(f) as fh:
            for line in fh:
                name, ln = line.split('\t')[:2]
                if chrom_of(name).startswith(chrom_prefix):
                    lengths[name] = int(ln)
    return lengths

def read_paf(files, min_mapq, min_len):
    alns = collections.defaultdict(list)
    for f in files:
        with open(f) as fh:
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) < 12:
                    continue
                if 'tp:A:S' in c[12:]:
                    continue
                if int(c[11]) < min_mapq or int(c[10]) < min_len:
                    continue
                alns[c[0]].append((int(c[2]), int(c[3]), chrom_of(c[5]), c[4]))
    return alns

def blocks_for(alist, min_block):
    alist = sorted(alist)
    blocks = []
    for qs, qe, t, s in alist:
        if blocks and blocks[-1]['t'] == t:
            b = blocks[-1]; b['e'] = max(b['e'], qe); b['bp'] += qe - qs; b[s] += qe - qs
        else:
            blocks.append({'s': qs, 'e': qe, 't': t, 'bp': qe - qs, '+': 0, '-': 0})
            blocks[-1][s] += qe - qs
    big = [b for b in blocks if b['bp'] >= min_block]
    merged = []
    for b in big:
        if merged and merged[-1]['t'] == b['t']:
            m = merged[-1]; m['e'] = b['e']; m['bp'] += b['bp']; m['+'] += b['+']; m['-'] += b['-']
        else:
            merged.append(dict(b))
    return merged

def mb(x):
    return f"{x/1e6:.2f}"

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--paf', nargs='+', required=True)
    ap.add_argument('--fai', nargs='+', required=True)
    ap.add_argument('--out', default='chrom_qc.tsv')
    ap.add_argument('--chrom-prefix', default='chr')
    ap.add_argument('--min-mapq', type=int, default=20)
    ap.add_argument('--min-aln', type=int, default=10000)
    ap.add_argument('--min-block', type=int, default=500000)
    ap.add_argument('--len-tol', type=float, default=0.15)
    ap.add_argument('--min-aln-frac', type=float, default=0.5)
    a = ap.parse_args()

    lengths = read_fai(a.fai, a.chrom_prefix)
    alns = read_paf(a.paf, a.min_mapq, a.min_aln)

    by_chrom = collections.defaultdict(list)
    for n, l in lengths.items():
        by_chrom[chrom_of(n)].append(l)
    med = {c: statistics.median(v) for c, v in by_chrom.items()}

    cols = ['haplotype', 'chrom', 'length_mb', 'len_vs_median', 'aligned_frac', 'frac_expected',
            'top_target', 'frac_top', 'frac_rev_expected', 'blocks', 'flags']
    rows = []
    for name in sorted(lengths, key=lambda n: (chrom_of(n), hap_of(n))):
        c, L = chrom_of(name), lengths[name]
        al = alns.get(name, [])
        per_t = collections.Counter(); rev_exp = 0
        for qs, qe, t, s in al:
            per_t[t] += qe - qs
            if t == c and s == '-':
                rev_exp += qe - qs
        tot = sum(per_t.values())
        top, top_bp = (per_t.most_common(1)[0] if per_t else ('NA', 0))
        blocks = blocks_for(al, a.min_block)
        bstr = ';'.join(f"{mb(b['s'])}-{mb(b['e'])}:{b['t']}({'+' if b['+'] >= b['-'] else '-'})" for b in blocks)

        ratio = L / med[c] if med[c] else float('nan')
        flags = []
        if abs(ratio - 1) > a.len_tol: flags.append('LEN_OUTLIER')
        if tot / L < a.min_aln_frac: flags.append('LOW_ALIGN')
        if top != c and top_bp > 0.5 * tot and tot: flags.append('MISLABEL')
        if len({b['t'] for b in blocks}) > 1: flags.append('CHIMERA')
        if per_t[c] and rev_exp > 0.5 * per_t[c]: flags.append('MOSTLY_REV')

        rows.append([hap_of(name), c, mb(L), f"{ratio:.2f}", f"{tot/L:.2f}",
                     f"{per_t[c]/L:.2f}", top, f"{(top_bp/tot if tot else 0):.2f}",
                     f"{(rev_exp/per_t[c] if per_t[c] else 0):.2f}", bstr or 'NA',
                     ','.join(flags) or 'OK'])

    with open(a.out, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for r in rows:
            out.write('\t'.join(r) + '\n')

    flagged = [r for r in rows if r[-1] != 'OK']
    print(f"{len(rows)} haplotype-chromosomes checked; {len(flagged)} flagged -> {a.out}", file=sys.stderr)
    for r in flagged:
        print(f"  {r[0]}#{r[1]}\t{r[2]} Mb (x{r[3]})\t{r[-1]}\t{r[9]}", file=sys.stderr)

if __name__ == '__main__':
    main()
