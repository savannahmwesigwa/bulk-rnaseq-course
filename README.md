# Bulk RNA-seq: From Reads to Results

A one-hour introduction to bulk RNA-seq analysis, for biologists and
clinicians with little or no bioinformatics background.

Everything runs on a regular laptop. No prior R experience needed.

---

## Start here

**[SETUP.md](SETUP.md)** — install R and RStudio, and check it works.
Do this at least two days before the course.

---

## What this covers

- What a count matrix is, and why raw counts cannot be compared directly
- CPM, RPKM, TPM — what each corrects for, and which to use when
- Checking whether your samples are what the labels say
- Detecting contamination from co-cultured cells
- Running a differential expression analysis with DESeq2
- Reading the output, and knowing when not to trust it

Upstream processing (alignment, quantification) is covered as background
only — this lab runs [nf-core/rnaseq](https://nf-co.re/rnaseq) for that.

---

## The data

Real published data: **[GSE191074](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191074)**
— human NK cells from healthy donors, five conditions, three biological
replicates each.

| Condition | What was in the well |
|---|---|
| `NK-only` | NK cells alone |
| `NK-PGE2` | + prostaglandin E2 |
| `NK-K562` | + K562 leukemia cell line |
| `NK-Blasts` | + primary AML blasts |
| `NK-Blasts-K562` | + both |

The count table is included in `data/`, so you do not need to download
anything during the course.

---

## Contents

```
00_setup.R              Run this first. Installs everything and checks it works.
SETUP.md                Installation instructions.
hands_on_tutorial.qmd   The practical. Run top to bottom.
hands_on_tutorial.R     The same code as a plain R script, if you prefer.
data/                   Count table from GEO.
slides/                 Theory presentation, plus an appendix of reference material.
docs/                   The practical rendered as a web page.
```

---

## Following along without R

The rendered practical is at
**[savannahmwesigwa.github.io/bulk-rnaseq-course](https://savannahmwesigwa.github.io/bulk-rnaseq-course)**
— all the code, output and figures, readable in a browser.

Useful if your installation breaks on the day, or if you just want to read
through it afterward.

---

## Questions

gushijima@mdanderson.org

---

## Data credit

All data from GEO accession GSE191074, generated and deposited by the
original study authors and made public for reuse.

The contamination discussed in this course is a normal technical feature of
co-culture experiments, not a criticism of the original work.
