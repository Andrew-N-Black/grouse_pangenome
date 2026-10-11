# =============================================================================
# 02d_plot_qc_report.R -- Merqury k-mer completeness and QV (and HiFi depth) of
# the 46 haplotype assemblies used as pangenome input.
# USFWS report Figure S8. Input: pangenome_qc_summary.tsv from 02c_qc_report_table.sh (copy it
# locally). Style matches plot_assembly_contiguity.R (Figure S5): one row per
# bird, grouped by species; hap1 / hap2 / both; W-bearing haplotypes open.
#   install.packages(c("ggplot2", "patchwork"))
# Outputs: FigS_merqury_qc.pdf / .png, and a summary printed for the text.
# =============================================================================
library(ggplot2)
infile <- "pangenome_qc_summary.tsv"
d <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE)
d$species <- trimws(d$species)
d$hap <- ifelse(d$level == "both", "both", sub(".*_(hap[12]).*$", "\\1", d$assembly))
d$id  <- paste(d$sample, d$hap, sep = "_")

# W-bearing haplotypes of the six females (Table S3; Figure S7)
w_hap <- c("F5503_hap1", "F5595_hap2", "F5596_hap2", "F5598_hap1", "F5599_hap1", "F5600_hap1")
d$pt <- ifelse(d$hap == "both", "both", ifelse(d$id %in% w_hap, paste0(d$hap, "_W"), d$hap))

d$species <- factor(d$species, levels = c("LEPC", "GRPC", "STGR"))
ord <- unique(d$sample[order(d$species, d$sample)])
d$sample <- factor(d$sample, levels = rev(ord))
sp_cols <- c(LEPC = "#A52A2A", GRPC = "#DAA520", STGR = "#000000")
shapes  <- c(hap1 = 16, hap2 = 17, hap1_W = 1, hap2_W = 2, both = 18)

panel <- function(dd, var, xlab, show_y = TRUE, show_strip = FALSE) {
  ggplot(dd, aes(x = .data[[var]], y = sample)) +
    geom_line(aes(group = sample, colour = species), linewidth = 0.4, alpha = 0.5) +
    geom_point(aes(colour = species, shape = pt), size = 2.2, stroke = 0.9, show.legend = TRUE) +
    scale_shape_manual(values = shapes, limits = names(shapes), breaks = c("hap1", "hap2", "hap1_W", "hap2_W", "both"),
                       labels = c("hap1", "hap2", "hap1, W-bearing", "hap2, W-bearing", "hap1 + hap2"), name = NULL) +
    scale_colour_manual(values = sp_cols, guide = "none") +
    facet_grid(species ~ ., scales = "free_y", space = "free_y") +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = if (show_y) element_text(size = 8) else element_blank(),
          axis.ticks.y = if (show_y) element_line() else element_blank(),
          strip.text.y.right = if (show_strip) element_text(angle = 0, face = "bold") else element_blank(),
          strip.background = if (show_strip) element_rect(fill = "grey95") else element_blank(),
          legend.position = "bottom")
}

pA <- panel(d, "completeness_pct", "k-mer completeness (%)") + labs(tag = "A")
# Only panel A carries the shape legend. ggplot2 >= 3.5 drops legend keys for
# levels a panel does not plot (panel C has no "hap1 + hap2"), so per-panel
# legends differ and patchwork cannot merge them into one.
pB <- panel(d, "QV", "Consensus quality (QV)", show_y = FALSE) + labs(tag = "B") + guides(shape = "none")
h  <- d[d$level == "haplotype" & !is.na(d$median_depth), ]
has_depth <- nrow(h) > 0
if (has_depth) pC <- panel(h, "median_depth", "Median HiFi depth (x)", show_y = FALSE, show_strip = TRUE) + labs(tag = "C") + guides(shape = "none")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- (if (has_depth) (pA | pB | pC) else (pA | pB)) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  print(fig)
  ggsave("FigS_merqury_qc.pdf", fig, width = 9, height = 6.5, device = cairo_pdf)
  ggsave("FigS_merqury_qc.png", fig, width = 9, height = 6.5, dpi = 600)
} else {
  ggsave("FigS_merqury_completeness.png", pA, width = 4.5, height = 6.5, dpi = 600)
  ggsave("FigS_merqury_QV.png", pB, width = 4, height = 6.5, dpi = 600)
}

# ---- numbers for the text ----------------------------------------------------
f <- function(x) sprintf("%.2f-%.2f (median %.2f)", min(x), max(x), median(x))
hh <- d[d$level == "haplotype", ]; bb <- d[d$level == "both", ]
cat("Completeness, hap1 + hap2:", f(bb$completeness_pct), "%\n")
cat("Completeness, haplotypes: ", f(hh$completeness_pct), "%\n")
cat("  W-bearing haplotypes:   ", f(hh$completeness_pct[hh$id %in% w_hap]), "%\n")
cat("  other haplotypes:       ", f(hh$completeness_pct[!hh$id %in% w_hap]), "%\n")
cat("QV, haplotypes:           ", f(hh$QV), "\n")
cat("QV, hap1 + hap2:          ", f(bb$QV), "\n")
if (has_depth) {
  cat("Median HiFi depth/hap:    ", f(h$median_depth), "x\n")
  cat("Low-depth bp (%):         ", f(h$low_bp_pct), "\n")
  cat("High-depth bp (%):        ", f(h$high_bp_pct), "\n")
}
print(aggregate(cbind(completeness_pct, QV) ~ species + level, d, median))
