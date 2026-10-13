#!/usr/bin/env python3
"""
chrom_qc.py - flag mislabelled, chimeric and length-outlier chromosomes in a set of
haplotype assemblies BEFORE building per-chromosome pangenome graphs.

Inputs
  --paf   one PAF per haplotype: minimap2 -x asm10 --secondary=no ANCHOR.fa HAP.fa
  --fai   samtools faidx .fai for every haplotype FASTA

The haplotype is identified from the FILE name (GRPC_F5595_hap1.pg_input.paf <->
GRPC_F5595_hap1.pg_input.fa.fai), so sequence names do not need to be PanSN.
Chromosome names are compared after dropping any PanSN prefix (up to the last '#')
and any split suffix (_p1, _p2, ...), so chr_10_p1 is checked against anchor chr_10.

Flags (one row per sequence)
  LEN_OUTLIER   summed length of this chromosome (all _pN pieces) differs from the
                cross-haplotype median by > --len-tol
  MISLABEL      most aligned bases go to a different anchor chromosome
  CHIMERA       >= --min-block aligned bp in one block to a different anchor chromosome
  LOW_ALIGN     < --min-aln-frac of the sequence aligns (sequences >= 200 kb only)
  MOSTLY_REV    most bases aligned to the expected chromosome are on the reverse strand
Anchor sequences not starting with --chrom-prefix (e.g. scaffold_12) are ignored.
"""
import argparse, collections, os, re, statistics, sys

SPLIT = re.compile(r'_p\d+$')
EXT = re.compile(r'(\.fai)?(\.(fa|fasta|fna))?(\.gz)?(\.paf)?$')

def base_chrom(name):
    return SPLIT.sub('', name.split('#')[-1])

def hap_id(path):
    b = os.path.basename(path)
    prev = None
    while prev != b:
        prev, b = b, EXT.sub('', b)
    return b

def read_fai(files, prefix, skip):
    seqs = {}                                       # (hap, seqname) -> length
    for f in files:
        h = hap_id(f)
        with open(f) as fh:
            for line in fh:
                name, ln = line.split('\t')[:2]
                c = base_chrom(name)
                if c.startswith(prefix) and c not in skip:
                    seqs[(h, name)] = int(ln)
    return seqs

def read_paf(files, min_mapq, min_len, prefix):
    alns = collections.defaultdict(list)            # (hap, qname) -> [(qs, qe, tchrom, strand)]
    for f in files:
        h = hap_id(f)
        with open(f) as fh:
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) < 12 or 'tp:A:S' in c[12:]:
                    continue
                if int(c[11]) < min_mapq or int(c[10]) < min_len:
                    continue
                t = base_chrom(c[5])
                if not t.startswith(prefix):
                    continue
                alns[(h, c[0])].append((int(c[2]), int(c[3]), t, c[4]))
    return alns

def blocks_for(alist):
    blocks = []
    for qs, qe, t, s in sorted(alist):
        if blocks and blocks[-1]['t'] == t:
            b = blocks[-1]; b['e'] = max(b['e'], qe); b['bp'] += qe - qs; b[s] += qe - qs
        else:
            blocks.append({'s': qs, 'e': qe, 't': t, 'bp': qe - qs, '+': 0, '-': 0})
            blocks[-1][s] += qe - qs
    return blocks

def mb(x):
    return f"{x/1e6:.2f}"

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--paf', nargs='+', required=True)
    ap.add_argument('--fai', nargs='+', required=True)
    ap.add_argument('--out', default='chrom_qc.tsv')
    ap.add_argument('--chrom-prefix', default='chr')
    ap.add_argument('--skip', default='chr_MT', help='comma-separated chromosomes to ignore')
    ap.add_argument('--min-mapq', type=int, default=20)
    ap.add_argument('--min-aln', type=int, default=10000)
    ap.add_argument('--min-block', type=int, default=1000000, help='aligned bp to another chromosome needed to call CHIMERA')
    ap.add_argument('--len-tol', type=float, default=0.15)
    ap.add_argument('--min-aln-frac', type=float, default=0.5)
    a = ap.parse_args()

    skip = set(a.skip.split(',')) if a.skip else set()
    seqs = read_fai(a.fai, a.chrom_prefix, skip)
    alns = read_paf(a.paf, a.min_mapq, a.min_aln, a.chrom_prefix)
    haps = sorted({h for h, _ in seqs})

    # length screen on whole chromosomes (all _pN pieces summed)
    tot_len = collections.Counter()
    for (h, n), l in seqs.items():
        tot_len[(h, base_chrom(n))] += l
    by_chrom = collections.defaultdict(list)
    for (h, c), l in tot_len.items():
        by_chrom[c].append(l)
    med = {c: statistics.median(v) for c, v in by_chrom.items()}

    cols = ['haplotype', 'sequence', 'chrom', 'seq_mb', 'chrom_total_mb', 'chrom_vs_median',
            'aligned_frac', 'frac_expected', 'top_target', 'frac_rev_expected', 'off_target_blocks', 'flags']
    rows = []
    for (h, n) in sorted(seqs, key=lambda k: (base_chrom(k[1]), k[0], k[1])):
        c, L = base_chrom(n), seqs[(h, n)]
        al = alns.get((h, n), [])
        per_t = collections.Counter(); rev_exp = 0
        for qs, qe, t, s in al:
            per_t[t] += qe - qs
            if t == c and s == '-':
                rev_exp += qe - qs
        tot = sum(per_t.values())
        top, top_bp = (per_t.most_common(1)[0] if per_t else ('NA', 0))
        off = [b for b in blocks_for(al) if b['t'] != c and b['bp'] >= a.min_block]
        ostr = ';'.join(f"{mb(b['s'])}-{mb(b['e'])}:{b['t']}({mb(b['bp'])}Mb aligned)" for b in off)
        ratio = tot_len[(h, c)] / med[c] if med[c] else float('nan')

        flags = []
        if abs(ratio - 1) > a.len_tol: flags.append('LEN_OUTLIER')
        if L >= 200000 and tot / L < a.min_aln_frac: flags.append('LOW_ALIGN')
        if tot and top != c and top_bp > 0.5 * tot: flags.append('MISLABEL')
        if off and 'MISLABEL' not in flags: flags.append('CHIMERA')
        if per_t[c] and rev_exp > 0.5 * per_t[c]: flags.append('MOSTLY_REV')

        rows.append([h, n, c, mb(L), mb(tot_len[(h, c)]), f"{ratio:.2f}", f"{tot/L:.2f}",
                     f"{per_t[c]/L:.2f}", top, f"{(rev_exp/per_t[c] if per_t[c] else 0):.2f}",
                     ostr or '-', ','.join(flags) or 'OK'])

    present = collections.defaultdict(set)
    for (h, c) in tot_len:
        present[c].add(h)
    for c in sorted(present):
        if len(present[c]) >= 0.5 * len(haps):
            for h in haps:
                if h not in present[c]:
                    rows.append([h, c, c, '0.00', '0.00', '0.00', '0.00', '0.00', 'NA', '0.00', '-', 'MISSING'])

    with open(a.out, 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for r in rows:
            out.write('\t'.join(r) + '\n')

    flagged = [r for r in rows if r[-1] != 'OK']
    print(f"{len(haps)} haplotypes, {len(rows)} sequences checked; {len(flagged)} flagged -> {a.out}", file=sys.stderr)
    for r in flagged:
        print(f"  {r[0]}  {r[1]}\t{r[3]} Mb (chrom total x{r[5]})\t{r[-1]}\t{r[10]}", file=sys.stderr)

if __name__ == '__main__':
    main()
