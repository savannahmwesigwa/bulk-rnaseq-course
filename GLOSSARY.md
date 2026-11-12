# Glossary

Terms from the bulk RNA-seq workshop, in plain language.

---

## The data

**Read** — one short sequence that came off the sequencer. A fragment of
a transcript, not a whole one.

**Library** — the prepared, sequenceable pool of fragments made from one
sample. "Preparing a library" is all the bench work between RNA and the
machine.

**Library size** — total reads for a sample. Also called sequencing
depth. Varies between samples for reasons that have nothing to do with
biology.

**Count** — how many reads landed in a given gene in a given sample. The
raw material for everything downstream.

**Count matrix** — genes down the side, samples across the top, a count
in every cell. What the pipeline hands you.

**Sample sheet** — a separate table saying what each sample is:
condition, donor, batch, sex, anything else that varied. The model can
only use structure it is told about.

**FASTQ** — the raw sequencer output file. Four lines per read: a name,
the sequence, a spacer, and a quality score per base.

**polyA selection** — keeping only RNA with a polyA tail, which is mostly
mRNA. Cheap and standard, but loses anything without a tail.

**Ribo-depletion** — removing ribosomal RNA specifically and keeping
everything else, including non-coding RNA. Better for degraded samples.

**Adapter** — a short known sequence attached to each fragment so the
sequencer can grab hold of it.

**Pipeline** — the automated steps from FASTQ to count matrix: quality
control, alignment, counting. This lab runs nf-core/rnaseq.

---

## Normalization

**Normalization** — putting samples on a comparable scale so differences
reflect biology rather than how much each was sequenced.

**Size factor** — the per-sample number DESeq2 works out and divides by.
It corrects both for depth and for one gene hogging the library.

**Composition effect** — when one gene takes a much bigger share of the
library in one sample, making every other gene look lower even though
nothing changed in the cell. Size factors exist to fix this.

**CPM** — counts per million. Divides by library size only. Fixes depth,
ignores gene length. Used for filtering and plots.

**RPKM / FPKM** — divides by depth *then* by gene length. FPKM is the
paired-end version. Flawed, because sample totals end up different, so
the same number means different things in different samples.

**TPM** — the same two divisions in the opposite order. Because the
scaling happens last, every sample sums to exactly one million, so a TPM
of 50 means the same thing everywhere. Preferred over RPKM.

**logCPM** — CPM on a log2 scale. Expression spans several orders of
magnitude, so logging is what makes it plottable.

**Spike-in** — a known quantity of synthetic RNA added to every sample,
used as an external reference when you expect a genuine global shift in
transcription that normal methods would mistake for an artifact.

---

## Design

**Biological replicate** — an independent instance of what you are
studying: a different donor, animal, or separately started culture. When
people say "n", this is what they mean.

**Technical replicate** — the same biological material measured twice.
Tells you how reproducible the machine is, which with modern sequencing
is "very", so they are usually unnecessary.

**Pseudoreplication** — treating measurements that are not independent as
if they were. Three wells from one donor counted as n=3 is the classic
case; thousands of cells from three mice counted as n=thousands is the
single-cell version.

**Batch** — anything shared by a group of samples that could affect the
measurement: processing day, reagent lot, operator, sequencing run,
extraction kit.

**Batch effect** — a systematic difference between samples caused by a
batch rather than by biology.

**Confounding** — when two variables carry the same information, so their
effects cannot be separated. All controls on Monday and all treated on
Tuesday means "treatment" and "day" are the same variable. No analysis
fixes this.

**Blocking** — spreading each condition across every batch, so the batch
effect hits all groups equally and can be estimated and removed.

**Paired design** — the same donor appears in every condition, so each
donor acts as their own control. Much more powerful than independent
groups.

**Covariate** — a variable other than your main question that you put in
the model so its effect is accounted for. Batch and donor are the usual
ones. A covariate only helps if it varies *within* your conditions — if
it tracks the condition exactly, it is confounding.

**Design matrix** — the table telling the model what each sample is. In R
you write it as a formula: `~ batch + condition`.

**Formula** — the `~ batch + condition` notation. Read right to left:
compare the conditions, after accounting for batch.

**Reference level** — the condition everything else is compared against.
Whichever level comes first in the factor, so the order you set matters.

**Coefficient** — a number the model estimated for every gene, usually
"how far is this condition from the reference". `resultsNames()` lists
them.

**Contrast** — the specific comparison you pull out of a fitted model,
such as treated versus control. Can be any pair, not just comparisons
against the reference.

**Interaction** — asking whether the effect of one factor depends on
another: does the drug work differently in knockouts than wild type?
Needs substantially more samples than a main effect.

---

## The statistics

**Poisson** — the distribution describing pure sampling randomness. Scoop
cookies from one uniform batch of dough and the chip counts vary by
luck alone. Its spread is fixed by its average: you do not get to choose
it.

**Overdispersion** — real replicates vary more than Poisson predicts,
because the underlying biology genuinely differs between them. Different
batches of dough, not just different scoops.

**Negative binomial** — Poisson sampling on top of a mean that itself
varies between replicates. One extra parameter (the dispersion) captures
that extra variation. This is the distribution DESeq2 and edgeR assume.

**Dispersion** — how much a gene's counts vary between replicates beyond
what sampling alone would cause. High dispersion means a noisy gene.

**Shrinkage** — pulling a shaky per-gene estimate toward the typical
behaviour of similar genes. Applied to dispersions, and separately to
fold changes.

**Empirical Bayes** — the formal name for shrinkage: using the whole
dataset to inform the estimate for each individual gene. It is what makes
n=3 workable.

**Fold change** — how much a gene changed, as a ratio.

**log2 fold change** — that ratio on a log2 scale, so +1 is two-fold up,
−1 is two-fold down, 0 is no change. Symmetric, which plain ratios are
not.

**baseMean** — a gene's average normalized count across all samples. One
number per gene, used as its expression level on plots.

**Wald test** — tests one specific comparison: is this effect different
from zero?

**Likelihood ratio test (LRT)** — compares a full model against a reduced
one and asks whether the dropped terms matter at all. Use it when the
question spans several coefficients, such as a time course.

**p-value** — the chance of seeing a result this extreme if nothing were
really happening.

**Multiple testing** — the problem that testing 20,000 genes at p < 0.05
gives about 1,000 false hits even when nothing is going on.

**FDR (false discovery rate)** — the expected proportion of your
significant list that is wrong. Controlling FDR at 5% means about 5% of
your hits are false.

**padj** — the adjusted p-value, corrected for multiple testing. The
column to use. Never filter on raw `pvalue`.

**Benjamini–Hochberg** — the standard method for FDR adjustment.

**Independent filtering** — DESeq2 automatically removing genes too
lowly expressed to reach significance, so the correction is less harsh
on the rest. Those genes come back with `padj = NA`. Not an error.

**lfcShrink** — shrinks fold change estimates toward zero when they are
unreliable, leaving well-supported ones alone. Does not change p-values,
only the effect size.

---

## Plots

**PCA** — reduces thousands of genes to a couple of axes capturing most
of the variation, so samples can be plotted as points. Used to spot
outliers, batch structure, and whether conditions separate at all.

**MA plot** — fold change against average expression. **M** = minus (a
log ratio is a subtraction of logs), **A** = average. Check the cloud
sits centred on zero.

**Volcano plot** — fold change on the x-axis, significance on the y-axis.
Shows size and confidence of change together.

**p-value histogram** — the distribution of raw p-values. The most
diagnostic plot in the analysis. Flat with a spike at zero is healthy; a
hill in the middle means something is missing from the model.

**Dispersion plot** — each gene's variability against its expression
level, showing the raw estimates, the fitted trend, and the shrunk final
values.

**Marker panel** — a set of genes a particular cell type expresses, used
as a score to check what is in a sample.

---

## Things you will meet next

**Gene set / pathway analysis** — asking which biological processes
changed, rather than which individual genes. The usual next step after a
gene list.

**Over-representation analysis** — takes your significant list and asks
whether any pathway appears more often than chance. Simple, but depends
on where you put the cutoff.

**GSEA / ranked methods** — uses the whole ranked gene list rather than a
cutoff. Generally preferred.

**Length bias** — longer genes collect more reads, so they reach
significance more easily, and pathways full of long genes can look
enriched for no biological reason. `goseq` corrects for it.

**Pseudobulk** — for single-cell data, summing counts within each sample
and cell type to get a bulk-like matrix, then running this exact
workflow. Your replicates are the samples, not the cells.

**Deconvolution** — estimating what mixture of cell types produced a bulk
sample.

**Vignette** — the worked tutorial that ships with a Bioconductor
package. Usually the best documentation available, and the place to look
for current behaviour rather than the original paper.
