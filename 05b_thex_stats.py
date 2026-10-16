#!/usr/bin/env python3
"""
Step 05b: statistics for the THEx window-tree results (run after 05_thex_phylogeny.sh).

    python 05b_thex_stats.py grouse_TreeViewer_input.xlsx chromosome_lengths.bed [out_prefix]

Runs in the THEx env (pandas, scipy, statsmodels, scikit-posthocs). Writes <out_prefix>.*.tsv
and prints a summary.

WHY NOT THEx's "Exploratory Statistical Testing" TAB
    It treats every 100 kb window as an independent observation. Neighbouring windows share
    linkage, recombination rate and ancestry, so a few thousand windows behave like far fewer
    independent data points and its p-values (e.g. 1e-33) are badly overstated. It also only
    offers ANOVA on the extra columns, and CladeSupport is partly defined by the labels
    themselves (a "sorted" label requires a well-resolved clade), so that comparison is
    circular. This script uses tests whose unit of replication is defensible:

    1. Z vs autosomes      the Z is ONE chromosome, so it is compared with the distribution
                           of the 31 autosomes (empirical rank), not by pooling windows.
    2. Chromosome size     one value per autosome (n = 31): Spearman correlation between
                           chromosome length / gene density and the share of each label.
    3. Gene density        window-level, but with chromosome as a cluster: a logistic GEE
       (window level)      (exchangeable within-chromosome correlation) of "label == X" on
                           genes per window, adjusting for chromosome length. Kruskal-Wallis
                           + Dunn's (Bonferroni) are reported too, marked as descriptive.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

xlsx, bed = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else "thex_stats"

SORTED = "((LEPC,GRPC),STGR)"
MIXED = "((LEPC,GRPC),STGR) LEPC/GRPC unsorted"
NOPAIR = "No species pair forms a clade"
MAIN = [SORTED, MIXED, NOPAIR]
SHORT = {SORTED: "sorted", MIXED: "LEPC_GRPC_mixed", NOPAIR: "no_pair"}

d = pd.read_excel(xlsx)
d = d[d["NewickTree"] != "NoTree"].copy()
d["label"] = np.where(d["TopologyID"].isin(MAIN), d["TopologyID"], "other")
L = pd.read_csv(bed, sep="\t").set_index("Chromosome")["End"]
d["chr_len_mb"] = d["Chromosome"].map(L) / 1e6
d["is_Z"] = d["Chromosome"].eq("chr-Z")

# ---------------------------------------------------------------- per chromosome
pc = pd.crosstab(d["Chromosome"], d["label"])
for t in MAIN:
    if t not in pc:
        pc[t] = 0
pc["n_windows"] = pc.sum(axis=1)
for t in MAIN:
    pc["frac_" + SHORT[t]] = pc[t] / pc["n_windows"]
pc["length_mb"] = pc.index.map(L) / 1e6
pc["genes_per_mb"] = d.groupby("Chromosome")["Genes"].sum() / pc["length_mb"]
pc = pc.sort_values("length_mb", ascending=False)
pc.to_csv(out + ".per_chromosome.tsv", sep="\t")
auto = pc.drop(index="chr-Z", errors="ignore")

print("=" * 78)
print("1. Z vs AUTOSOMES (Z compared with the distribution of %d autosomes)" % len(auto))
rows = []
if "chr-Z" in pc.index:
    for t in MAIN:
        col = "frac_" + SHORT[t]
        z = pc.loc["chr-Z", col]
        a = auto[col]
        # empirical one-sided rank: how many autosomes are as extreme as the Z
        n_ge, n_le = (a >= z).sum(), (a <= z).sum()
        rows.append({"label": SHORT[t], "chr_Z": round(z, 4), "autosome_median": round(a.median(), 4),
                     "autosome_min": round(a.min(), 4), "autosome_max": round(a.max(), 4),
                     "autosomes_>=_Z": int(n_ge), "autosomes_<=_Z": int(n_le),
                     "empirical_p_one_sided": round((min(n_ge, n_le) + 1) / (len(a) + 1), 4)})
    r1 = pd.DataFrame(rows)
    print(r1.to_string(index=False))
    r1.to_csv(out + ".Z_vs_autosomes.tsv", sep="\t", index=False)
    print("  -> empirical p = (autosomes at least as extreme + 1) / (n autosomes + 1); with 31 autosomes the")
    print("     smallest possible value is 1/32 = 0.031. Report it as 'the Z exceeds every autosome'.")

print("=" * 78)
print("2. CHROMOSOME SIZE / GENE DENSITY (autosomes, one value per chromosome, Spearman)")
rows = []
for t in MAIN:
    col = "frac_" + SHORT[t]
    for x in ("length_mb", "genes_per_mb"):
        rho, p = stats.spearmanr(auto[x], auto[col])
        rows.append({"label_share": SHORT[t], "vs": x, "spearman_rho": round(rho, 3), "p": p, "n": len(auto)})
r2 = pd.DataFrame(rows)
print(r2.to_string(index=False))
r2.to_csv(out + ".chromosome_correlations.tsv", sep="\t", index=False)
rho, p = stats.spearmanr(auto["length_mb"], auto["genes_per_mb"])
print("  (length vs gene density across autosomes: rho = %.3f, p = %.3g -- the two predictors are"
      " confounded; see the GEE in 3)" % (rho, p))

print("=" * 78)
print("3. GENE DENSITY PER WINDOW (autosomal windows, chromosome as cluster)")
import statsmodels.api as sm
import statsmodels.formula.api as smf
a = d[~d["is_Z"]].copy()
a["log_len"] = np.log(a["chr_len_mb"])
rows = []
for t in MAIN:
    a["y"] = (a["label"] == t).astype(int)
    m = smf.gee("y ~ Genes + log_len", groups="Chromosome", data=a,
                family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
    for term in ("Genes", "log_len"):
        rows.append({"outcome": SHORT[t], "term": term, "odds_ratio": round(np.exp(m.params[term]), 3),
                     "OR_95CI": "%.3f-%.3f" % tuple(np.exp(m.conf_int().loc[term])),
                     "p": m.pvalues[term]})
r3 = pd.DataFrame(rows)
print(r3.to_string(index=False))
r3.to_csv(out + ".gee_window_level.tsv", sep="\t", index=False)
print("  -> odds ratio per extra gene in the window (Genes) and per e-fold larger chromosome (log_len),")
print("     with robust SEs clustered by chromosome. Genes > 1 and log_len < 1 for no_pair = mixing rises")
print("     with gene density and on smaller chromosomes, each after adjusting for the other.")

print("=" * 78)
print("   Descriptive only (windows treated as independent; p-values overstated):")
g = [a.loc[a["label"] == t, "Genes"] for t in MAIN]
print("   median genes/window: " + ", ".join("%s %.1f (mean %.2f, n=%d)" % (SHORT[t], x.median(), x.mean(), len(x))
                                           for t, x in zip(MAIN, g)))
h, p = stats.kruskal(*g)
print("   Kruskal-Wallis H = %.1f, p = %.3g" % (h, p))
try:
    import scikit_posthocs as sp
    dn = sp.posthoc_dunn(a[a["label"].isin(MAIN)], val_col="Genes", group_col="label", p_adjust="bonferroni")
    dn.index = dn.columns = [SHORT[c] for c in dn.columns]
    print("   Dunn's test (Bonferroni):\n" + dn.to_string().replace("\n", "\n   "))
except ImportError:
    pass
print("=" * 78)
print("tables: %s.{per_chromosome,Z_vs_autosomes,chromosome_correlations,gee_window_level}.tsv" % out)
