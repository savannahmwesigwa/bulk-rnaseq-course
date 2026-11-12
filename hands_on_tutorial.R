## ----setup------------------------------------------------------------------------------------------
#| include: false
knitr::opts_chunk$set(comment = "#>", collapse = FALSE)


## ----packages---------------------------------------------------------------------------------------
library(DESeq2)     # the analysis
library(ggplot2)    # plots
library(dplyr)      # handling tables
library(stringr)    # handling text
library(tidyr)      # reshaping tables
library(GEOquery)   # only needed if the data has to be downloaded


## ----load-counts------------------------------------------------------------------------------------
count_file <- "data/GSE191074_CountTable.txt.gz"

if (!file.exists(count_file)) {
  stop("Cannot find the count table.\n",
       "Open bulk-rnaseq-course.Rproj first, then reopen this document.\n",
       "Working directory is currently: ", getwd())
}

counts <- read.delim(count_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)

dim(counts)


## ----peek-------------------------------------------------------------------------------------------
counts[1:5, 1:4]


## ----sample-table-----------------------------------------------------------------------------------
design_df <- read.csv("data/sample_sheet.csv") |>
  mutate(
    condition = factor(condition,
                       levels = c("NK-only", "NK-PGE2", "NK-K562",
                                  "NK-Blasts", "NK-Blasts-K562")),
    replicate = factor(replicate)
  )

design_df


## ----check-order------------------------------------------------------------------------------------
identical(colnames(counts), as.character(design_df$sample))


## ----design-table-----------------------------------------------------------------------------------
table(design_df$condition, design_df$replicate)


## ----lib-sizes--------------------------------------------------------------------------------------
lib <- tibble(sample   = design_df$sample,
              millions = colSums(counts) / 1e6)

ggplot(lib, aes(x = sample, y = millions)) +
  geom_col(fill = "grey40") +
  labs(x = NULL, y = "Library size (millions)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----lcpm-------------------------------------------------------------------------------------------
# the one place edgeR is used: its cpm() is a convenient way to get
# normalized log counts for plotting. DESeq2 does everything else.
dge  <- edgeR::DGEList(counts)
dge  <- edgeR::calcNormFactors(dge)      # the composition correction
lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)


## ----log-transform----------------------------------------------------------------------------------
cpm_raw <- edgeR::cpm(dge)      # normalized, but not logged

bind_rows(
  as_tibble(cpm_raw) |>
    pivot_longer(everything(), names_to = "sample", values_to = "value") |>
    mutate(scale = "Raw CPM"),
  as_tibble(lcpm) |>
    pivot_longer(everything(), names_to = "sample", values_to = "value") |>
    mutate(scale = "log2 CPM")
) |>
  mutate(scale = factor(scale, levels = c("Raw CPM", "log2 CPM"))) |>
  ggplot(aes(x = sample, y = value)) +
  geom_boxplot(fill = "grey90", outlier.size = 0.3, outlier.alpha = 0.2) +
  facet_wrap(~ scale, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


## ----pca--------------------------------------------------------------------------------------------
# use the 500 most variable genes - the ones carrying the structure
top_var <- lcpm[order(apply(lcpm, 1, var), decreasing = TRUE)[1:500], ]

pca <- prcomp(t(top_var))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2))

pca_df <- as_tibble(pca$x[, 1:2]) |>
  mutate(condition = design_df$condition,
         replicate = design_df$replicate)

ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, shape = replicate)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1: ", pct[1], "% of variance"),
       y = paste0("PC2: ", pct[2], "% of variance")) +
  theme_minimal(base_size = 12)


## ----score-fn---------------------------------------------------------------------------------------
score_panel <- function(genes, mat) {
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[present, , drop = FALSE])
}


## ----sex-check--------------------------------------------------------------------------------------
y_genes <- c("RPS4Y1", "UTY", "KDM5D", "DDX3Y", "EIF1AY", "USP9Y")

sex <- tibble(
  sample  = colnames(lcpm),
  xist    = score_panel("XIST", lcpm),
  y_score = score_panel(y_genes, lcpm)
)

ggplot(sex, aes(x = xist, y = y_score)) +
  geom_point(size = 3, colour = "grey30") +
  expand_limits(x = c(-1, 8), y = c(0, 8)) +
  annotate("text", x = 6, y = 1, label = "female samples\nwould sit here",
           colour = "#C00000", size = 4) +
  labs(x = "XIST (log2 CPM)", y = "Y-gene score (log2 CPM)") +
  theme_minimal(base_size = 12)


## ----sex-raw----------------------------------------------------------------------------------------
counts[c("XIST", "RPS4Y1", "DDX3Y"), 1:6]


## ----panels-----------------------------------------------------------------------------------------
panels <- list(
  nk        = c("NCAM1", "KLRD1", "KLRF1", "NCR1", "NKG7", "GNLY", "PRF1"),
  t_cell    = c("TRAC"),
  monocyte  = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1"),
  b_cell    = c("MS4A1", "CD79A"),
  erythroid = c("HBG1", "HBG2", "HBZ", "GATA1", "GYPA", "KLF1"),
  platelet  = c("PPBP", "PF4")
)

scores <- as_tibble(lapply(panels, score_panel, mat = lcpm)) |>
  mutate(sample    = colnames(lcpm),
         condition = design_df$condition)


## ----nk-confirm-------------------------------------------------------------------------------------
round(range(scores$nk), 1)


## ----contamination----------------------------------------------------------------------------------
scores |>
  pivot_longer(cols = all_of(setdiff(names(panels), "nk")),
               names_to = "panel", values_to = "score") |>
  ggplot(aes(x = condition, y = score)) +
  geom_point(size = 2.5, colour = "grey30") +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = NULL, y = "Marker score (mean log2 CPM)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


## ----ery-by-cond------------------------------------------------------------------------------------
scores |>
  group_by(condition) |>
  summarise(erythroid = round(mean(erythroid), 2))


## ----hb-raw-----------------------------------------------------------------------------------------
counts[c("HBG2", "HBZ", "HBG1", "HBB"),
       c("NK-only-BR1.fastq", "NK-only-BR2.fastq",
         "NK-K562-BR1.fastq", "NK-K562-BR2.fastq")]


## ----dds--------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = design_df,
  design    = ~ replicate + condition
)


## ----filter-----------------------------------------------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)

dds <- dds[keep, ]


## ----run-deseq--------------------------------------------------------------------------------------
dds <- DESeq(dds)


## ----dispersion-------------------------------------------------------------------------------------
plotDispEsts(dds)


## ----results-names----------------------------------------------------------------------------------
resultsNames(dds)


## ----res-pge2---------------------------------------------------------------------------------------
res_pge2 <- results(dds, contrast = c("condition", "NK-PGE2", "NK-only"))

sum(res_pge2$padj < 0.05, na.rm = TRUE)


## ----top-pge2---------------------------------------------------------------------------------------
as.data.frame(res_pge2) |>
  tibble::rownames_to_column("gene") |>
  filter(!is.na(padj)) |>
  arrange(padj) |>
  select(gene, baseMean, log2FoldChange, padj) |>
  head(10)


## ----crem-------------------------------------------------------------------------------------------
match("CREM", rownames(res_pge2[order(res_pge2$padj), ]))


## ----pval-hist--------------------------------------------------------------------------------------
hist(res_pge2$pvalue, breaks = 50, col = "grey80", border = "white",
     main = "NK-PGE2 vs NK-only", xlab = "raw p-value")


## ----ma-plot----------------------------------------------------------------------------------------
# DESeq2:: spelled out: limma also has a plotMA, and if limma is loaded
# it masks this one and silently plots something else
DESeq2::plotMA(res_pge2, ylim = c(-4, 4))


## ----volcano-pge2-----------------------------------------------------------------------------------
volc <- as.data.frame(res_pge2) |>
  tibble::rownames_to_column("gene") |>
  filter(!is.na(padj)) |>
  mutate(sig = padj < 0.05)

# the ten most significant genes, to label
top10 <- volc |> arrange(padj) |> head(10)

ggplot(volc, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = sig), alpha = 0.6, size = 1) +
  ggrepel::geom_text_repel(data = top10, aes(label = gene),
                           size = 3.5, min.segment.length = 0,
                           max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c("FALSE" = "grey80", "TRUE" = "steelblue"),
                      guide = "none") +
  labs(x = "log2 fold change", y = "-log10 adjusted p") +
  theme_minimal(base_size = 12)


## ----res-k562---------------------------------------------------------------------------------------
res_k562 <- results(dds, contrast = c("condition", "NK-K562", "NK-only"))

sum(res_k562$padj < 0.05, na.rm = TRUE)


## ----top-k562---------------------------------------------------------------------------------------
as.data.frame(res_k562) |>
  tibble::rownames_to_column("gene") |>
  filter(!is.na(padj)) |>
  arrange(padj) |>
  select(gene, log2FoldChange, padj) |>
  head(10)


## ----passenger-ranks--------------------------------------------------------------------------------
k562_genes <- c("HBG1", "HBG2", "HBZ", "HBA1", "HBA2", "HBB",
                "GATA1", "GYPA", "GYPB", "KLF1")

as.data.frame(res_k562) |>
  tibble::rownames_to_column("gene") |>
  filter(!is.na(padj)) |>
  arrange(padj) |>
  mutate(rank = row_number()) |>
  filter(gene %in% k562_genes) |>
  select(gene, rank, log2FoldChange, padj)


## ----negative-control-------------------------------------------------------------------------------
as.data.frame(res_pge2) |>
  tibble::rownames_to_column("gene") |>
  filter(gene %in% k562_genes) |>
  select(gene, log2FoldChange, padj) |>
  arrange(padj)


## ----block-test-------------------------------------------------------------------------------------
# this refits the whole model, so it takes about 30 seconds
dds_nb <- DESeqDataSetFromMatrix(counts[keep, ], design_df, ~ condition)
dds_nb <- DESeq(dds_nb)
res_nb <- results(dds_nb, contrast = c("condition", "NK-PGE2", "NK-only"))

c(blocked   = sum(res_pge2$padj < 0.05, na.rm = TRUE),
  unblocked = sum(res_nb$padj   < 0.05, na.rm = TRUE))


## ----make-sheet-------------------------------------------------------------------------------------
#| eval: false
# library(stringr)
# 
# tibble(sample = colnames(counts)) |>
#   mutate(
#     # drop the "-BR1.fastq" ending, leaving the condition
#     condition = str_remove(sample, "-BR[0-9]+\\.fastq$"),
#     # pull out just the BR number
#     replicate = str_extract(sample, "BR[0-9]+")
#   ) |>
#   write.csv("data/sample_sheet.csv", row.names = FALSE)


## ----session-info-----------------------------------------------------------------------------------
sessionInfo()

