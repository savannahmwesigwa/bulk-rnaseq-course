# =============================================================
# BULK RNA-SEQ COURSE - SETUP SCRIPT
#
# Run this ONCE, a few days before the course.
#
# HOW TO RUN:
#   1. Open RStudio
#   2. File > Open File... and choose this file
#   3. Click "Source" (top right of this pane), or press Ctrl+Shift+S
#      (Cmd+Shift+S on Mac)
#   4. Wait. First install can take 10-20 minutes.
#   5. Read the PASS/FAIL report at the end.
#
# If anything says FAIL, copy the whole console output and email it to
# gushijima@mdanderson.org before the course.
# =============================================================

cat("\n========================================\n")
cat("  Bulk RNA-seq course - setup\n")
cat("========================================\n\n")


# ---- 1. Check R version -----------------------------------------------------

r_ver <- getRversion()
cat("R version:", as.character(r_ver), "\n")
cat("Platform :", R.version$platform, "\n\n")

if (r_ver < "4.4.0") {
  cat("*** YOUR R IS TOO OLD ***\n")
  cat("Please install R 4.4 or newer before continuing:\n")
  cat("  https://cloud.r-project.org\n")
  cat("Then restart RStudio and run this script again.\n\n")
  stop("R version too old. Nothing was installed.")
}


# ---- 2. Prefer pre-built binaries -------------------------------------------
# Compiling from source is where most installs fail. These settings tell R to
# download ready-made binaries instead.

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  install.packages.compile.from.source = "never",
  Ncpus = max(1, parallel::detectCores() - 1)
)


# ---- 3. Install BiocManager --------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}


# ---- 4. Install everything we need -------------------------------------------

cran_pkgs <- c(
  "ggplot2",      # plotting
  "dplyr",        # handling tables
  "stringr",      # handling text
  "tidyr",        # reshaping tables
  "ggrepel",      # non-overlapping labels on plots
  "pheatmap",     # heatmaps
  "gridExtra",    # arranging plots
  "RColorBrewer"  # color palettes
)

bioc_pkgs <- c(
  "DESeq2",       # differential expression - the main tool
  "edgeR",        # used for one normalization step
  "limma",        # comes with edgeR; not used directly
  "apeglm",       # shrinking fold changes
  "GEOquery",     # downloading data from GEO
  "AnnotationDbi",# gene ID conversion
  "org.Hs.eg.db"  # human gene annotation
)

missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]

if (length(missing_cran) > 0) {
  cat("\nInstalling from CRAN:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran)
}

if (length(missing_bioc) > 0) {
  cat("\nInstalling from Bioconductor:", paste(missing_bioc, collapse = ", "), "\n")
  cat("(This is the slow part. Go and make a coffee.)\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}


# ---- 5. Check Quarto ---------------------------------------------------------
# Quarto ships with recent RStudio. We only need it if you want to re-render
# the course document yourself - not required to follow along.

quarto_ok <- nzchar(Sys.which("quarto"))


# ---- 6. Report ---------------------------------------------------------------

cat("\n\n========================================\n")
cat("  RESULTS\n")
cat("========================================\n\n")

all_pkgs <- c(cran_pkgs, bioc_pkgs)
results <- data.frame(
  package = all_pkgs,
  status  = NA_character_,
  version = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(all_pkgs)) {
  p <- all_pkgs[i]
  ok <- suppressWarnings(suppressMessages(
    requireNamespace(p, quietly = TRUE)
  ))
  results$status[i]  <- if (ok) "PASS" else "FAIL"
  results$version[i] <- if (ok) as.character(packageVersion(p)) else "-"
}

for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-6s  %-16s %s\n",
              results$status[i], results$package[i], results$version[i]))
}

cat("\n  ", if (quarto_ok) "PASS  " else "note  ",
    "quarto          ",
    if (quarto_ok) system("quarto --version", intern = TRUE) else
      "not found (optional - not needed to follow along)", "\n", sep = "")


# ---- 7. Functional test ------------------------------------------------------
# Installed is not the same as working. Actually run something.

cat("\n----------------------------------------\n")
cat("  Functional test\n")
cat("----------------------------------------\n\n")

test_ok <- tryCatch({
  suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(ggplot2)
  })

  set.seed(1)
  m <- matrix(rnbinom(2000, mu = 100, size = 10), ncol = 4)
  rownames(m) <- paste0("gene", 1:500)
  colnames(m) <- paste0("s", 1:4)
  cd <- data.frame(grp = factor(c("a", "a", "b", "b")))

  d <- DESeqDataSetFromMatrix(m, cd, ~ grp)
  d <- DESeq(d, quiet = TRUE)
  r <- results(d)

  stopifnot(nrow(r) == 500)
  TRUE
}, error = function(e) {
  cat("FAILED with error:\n", conditionMessage(e), "\n")
  FALSE
})

cat(if (test_ok) "  PASS  DESeq2 ran a small test analysis successfully\n"
    else         "  FAIL  DESeq2 could not complete a test analysis\n")


# ---- 8. Verdict --------------------------------------------------------------

n_fail <- sum(results$status == "FAIL")

cat("\n========================================\n")
if (n_fail == 0 && test_ok) {
  cat("  ALL GOOD - you are ready for the course.\n")
} else {
  cat("  SOMETHING DID NOT WORK.\n\n")
  cat("  Please copy EVERYTHING in this console and\n")
  cat("  email it to gushijima@mdanderson.org before\n")
  cat("  the course. Do not wait until the day.\n")
}
cat("========================================\n\n")

cat("For the organizer:\n")
print(sessionInfo())
