# Questions and answers

Questions that come up about bulk RNA-seq, grouped by topic. Starred
ones (★) are the ones people ask most often.

If yours is not here, email gushijima@mdanderson.org.

---

## Design

**★ What if I can't do a paired design? My samples are from different patients and each only gives material for one condition.**

Then it is unpaired and that is fine — it just costs power, so you need
more samples to see the same effect. What you must not do is pretend
unpaired samples are paired. If some donors give both conditions and
others only one, say so: the model can use the paired ones as paired and
still include the rest.

**We already ran an experiment where all the controls were done first. Is it useless?**

For that comparison, largely yes — the treatment effect and the batch
effect cannot be separated. Three things can still be salvaged. Any
sub-comparison within a single batch is clean. The data gives you an
effect size estimate for powering a proper repeat. And if even two or
three samples break the pattern, those are worth a lot — bring them.

**What actually counts as a batch?**

Anything shared by a group of samples that could plausibly affect the
measurement: processing day, reagent lot, who did the prep, sequencing
run, flow cell, extraction kit, passage number. You do not need to model
all of them — you need to *record* all of them, then we can check which
ones matter.

**How do I know if 3 replicates is enough for my effect?**

Honestly, usually you do not until afterward. Rough guide: a clean cell
line with a strong stimulus, 3 is often plenty. Primary cells from
different donors with a subtle treatment, 3 will probably disappoint.
If you have pilot data we can estimate it properly. If the question is
an interaction, assume 3 is not enough.

**Should I randomize the order samples go on the sequencer?**

Yes, and it costs nothing. Same logic as the processing day: if all your
controls end up in one lane, lane effects become treatment effects.

**One of my samples looks like an outlier on PCA. Can I drop it?**

Only with a technical reason — low library size, failed QC, obvious
contamination, clear sample swap. Not because it is inconvenient.
Dropping samples that disagree with your hypothesis is how irreproducible
results happen. Whatever you decide, report it in the methods.

---

## Counts and normalization

**★ My collaborator sent me a TPM table and wants differential expression. What do I do?**

Ask for the counts. They exist — TPM is computed from them. If they
genuinely cannot be recovered, there are methods that work on
continuous data, but they are weaker and you should say so in the
methods.

**★ Can I compare my data to a published dataset?**

Carefully, and usually not by merging count matrices. Study-to-study
differences — library prep, pipeline, batch — are routinely larger than
the biological difference you care about. Two defensible options:
re-process both from raw data with the same pipeline and include study
as a batch term; or compare at the level of gene lists and pathways
rather than counts. Pasting two count matrices together and running
DESeq2 is not one of the options.

**What about spike-ins — do I need them?**

Only if you expect a genuine global shift in transcription, where the
usual normalization assumption breaks. Most experiments do not need
them. If you think yours might, that is a design-stage conversation.

**How many reads per sample?**

20–30 million for standard differential expression in human. More if you
care about lowly expressed genes, novel transcripts or splicing. Beyond
about 60 million you are buying very little for standard DE.

---

## Statistics and results

**★ What fold change cutoff should I use?**

The usual instinct is to filter the results table:

```r
res |> filter(padj < 0.05, abs(log2FoldChange) > 1)
```

The problem is that a fold change is an *estimate*, and some estimates
are far shakier than others. The filter treats them all as exact.

| | baseMean | estimated log2FC | 95% CI |
|---|---|---|---|
| GeneA | 20 | 1.3 | 0.1 to 2.5 |
| GeneB | 5,000 | 0.95 | 0.79 to 1.11 |

GeneA passes the filter, GeneB fails. But GeneA's true change could be
almost anything, while GeneB is measured precisely and sits right at the
boundary. The filter kept the unreliable one and dropped the reliable
one.

Build the threshold into the test instead:

```r
results(dds, contrast = ..., lfcThreshold = 1)
```

This changes the null hypothesis from "the change is zero" to "the change
is no bigger than two-fold", and returns p-values and FDR for *that*
question. Expect fewer genes than the naive filter gives — that is the
method being honest rather than strict.

If you only want to rank by fold change rather than threshold it, run
`lfcShrink()` first so the estimates are better behaved.

**★ I got 5,000 significant genes. Is that too many?**

Not necessarily. A strong perturbation genuinely changes a lot — our own
K562 comparison gave over 1,700. Before celebrating, check three things:
the p-value histogram looks healthy, the comparison is not confounded
with a batch, and something recognizable is near the top. If all three
hold, the number is the number.

**Nothing came out significant. What now?**

Look at the p-value histogram first. Flat with no spike means there is
genuinely nothing detectable — which is a result, not a failure. Then
check the PCA: if conditions do not separate at all, the effect is
smaller than the noise. Then check the design is right. What not to do
is relax the threshold until something appears.

**Should I use FDR 0.05 or 0.1?**

Either, decided before you look. 0.05 is conventional; 0.1 is defensible
for exploratory work where you will validate the hits anyway. Choosing
after seeing the results is not defensible.

**Do I still need qPCR validation?**

Reviewers often ask for it. Technically it is odd — RNA-seq is the more
quantitative measurement, so you are validating a better method with a
worse one. But it does independently confirm the RNA and the biology,
it is cheap, and it is expected. Do it for the handful of genes your
story rests on.

**How do I choose which genes to follow up?**

Not by p-value alone. Weigh effect size, whether the gene is expressed
at a level you could detect at the bench, whether it sits in a pathway
with other hits, and prior biology. A gene ranked 400th that fits a
coherent story often matters more than the top hit.

---

## Tools and next steps

**Should I use DESeq2 or edgeR?**

Either. They agree on the large majority of genes. Pick one and learn it
properly rather than switching between them. DESeq2's vignette is the
gentler read.

**What do I do with the gene list?**

Pathway or gene set analysis is the usual next step — it asks which
biological processes moved rather than which individual genes. There is
a slide on it in the appendix. It is a session of its own.

**★ Does any of this apply to single-cell data?**

More than you would think, and the main lesson is the same one as the
replicate slide.

The instinct in Seurat is to pool all cells from the treated samples,
all cells from the control, and run `FindMarkers`. With 5,000 cells a
side you get p-values around 1e-200. But cells from one mouse are not
independent of each other — 5,000 cells from 3 mice is n=3, not n=5,000.
It is the same error as counting three wells from one donor as three
replicates, at a scale where it is easier to miss.

The fix is **pseudobulk**: for each sample and cell type, sum the raw
counts across its cells. Three mice and four cell types gives a matrix
with twelve columns, which looks exactly like bulk RNA-seq. Run the
DESeq2 workflow from today on it. The mice are the replicates.

Where to draw the line:

| Question | Method |
|---|---|
| Which genes mark cluster 3 vs cluster 5? | cell-level tests are fine |
| Which genes change in NK cells between conditions? | pseudobulk |

This was genuinely contested until Squair et al. (2021, *Nature
Communications*) compared methods head to head and found that
single-cell-specific DE methods produce large numbers of false positives
for across-condition comparisons, while pseudobulk performs well.

**Can I use this for ATAC-seq?**

Largely yes. Peaks instead of genes, but the count matrix, the
normalization and the DESeq2/edgeR machinery are the same. Some details
differ — peak calling, and normalization is more delicate. The lab has a
Shiny app that covers both.

**★ Can you just run this for me?**

Yes, and talk to us before you design the experiment rather than after
you sequence. What to send: the raw counts, a sample sheet with every
piece of structure you know about (condition, donor, batch, sex,
processing date), and what you are actually trying to find out.

---

## From someone with a statistics or bioinformatics background

**Why block on batch in the design rather than batch-correcting the counts first with ComBat or `removeBatchEffect`?**

Because those tools are for *visualisation*, not for feeding into a test.
If you correct the data and then test it, the test has no idea that the
correction itself was estimated from the data, so the standard errors are
too small and the p-values are anticonservative. Putting batch in the
design matrix lets that uncertainty propagate properly.

The legitimate use of `removeBatchEffect` is producing a PCA or heatmap
that shows biology rather than batch. Never hand its output to DESeq2.

**When should I use `test = "LRT"` instead of the Wald test?**

Wald tests one coefficient or contrast: "is this particular effect
non-zero?" The LRT compares a full model against a reduced one and asks
whether the dropped terms collectively explain variance — the analogue of
an ANOVA F-test.

Use the LRT when the question spans several coefficients at once: does
condition matter at all across five levels, or does expression change
over a time course in any shape. Use Wald for specific pairwise
comparisons, which is most of what people want.

**Is the negative binomial actually right? What about zero inflation?**

For bulk RNA-seq, yes — bulk data shows little excess zero beyond what
the NB already predicts, and zero-inflated models mostly add parameters
without improving anything.

For single cell this was argued over for years. The current view
(Svensson 2020; Townes et al. 2019) is that UMI-based protocols are also
not zero-inflated relative to NB expectations, and the apparent
inflation in earlier data was largely an artefact of non-UMI protocols
and of normalisation choices.

**Does independent filtering bias the FDR?**

No, and the reason is specific: the filter statistic is `baseMean`,
which under the null hypothesis is independent of that gene's p-value.
You are not peeking at the test statistic, so filtering on it does not
distort the null distribution of what remains. Bourgon et al. (2010,
*PNAS*) is the argument. It raises power by shrinking the multiple
testing burden.

**Is BH valid when genes are correlated?**

Strictly, BH is proven under independence and under positive regression
dependence. Gene-gene correlations are mostly positive, since
co-regulated genes move together, so BH is generally accepted in
practice. Benjamini–Yekutieli is valid under arbitrary dependence but is
substantially more conservative and is rarely used for this.

**apeglm, ashr or normal for `lfcShrink`?**

`normal` is the original prior from the 2014 paper. It shrinks large,
well-measured fold changes too aggressively, which is why it is no
longer the default.

`apeglm` (Zhu et al. 2019) uses a heavy-tailed prior, so genuinely large
effects survive while noisy small ones are pulled in. It is the default
and it requires `coef`, not `contrast`.

`ashr` (Stephens 2017) uses adaptive shrinkage with a mixture prior and
accepts arbitrary contrasts. Reach for it when the comparison you want
cannot be written as a coefficient name.

**I have repeated measures — several samples per donor per condition. Shouldn't that be a mixed model?**

In principle yes. Three practical options, in increasing order of effort:

Sum the technical replicates within donor and condition
(`collapseReplicates`), which is conservative and simple. Use
`limma::duplicateCorrelation`, which estimates a single consensus
within-donor correlation and uses it — cheap and well established. Or
fit a proper GLMM with `glmmSeq` or NEBULA.

For a straightforward paired design, putting donor in the fixed-effects
design is usually enough and loses less power than people expect.

**How would I do a power calculation?**

`RNASeqPower` or `ssizeRNA`, or simulate from pilot data. All of them
need a dispersion estimate, which is the part you usually do not have —
so in practice it means pilot data, or borrowing dispersions from a
published dataset in a similar system. Anyone quoting a power
calculation without a dispersion estimate is quoting an assumption.

**What about SVA or RUV for batch effects I don't know about?**

Both estimate unwanted variation from the data and give you variables to
put in the design. SVA infers surrogate variables; RUV uses control
genes, either housekeeping genes or empirically derived negative
controls.

The caution is that they can absorb real signal. Pass your variable of
interest to `sva()` so it is protected, and always check what the
surrogate variables correlate with before trusting them — if a surrogate
tracks your condition, something has gone wrong.

**Does it matter whether counts came from Salmon or STAR + featureCounts?**

For gene-level DE they agree closely. Salmon distributes multi-mapping
reads probabilistically, which matters for gene families and paralogs
where featureCounts would discard them.

The nf-core file worth knowing about is
`salmon.merged.gene_counts_length_scaled.tsv`. The length scaling
corrects for effective transcript length differing between samples
because of differential isoform usage — the issue the tximport paper
(Soneson et al. 2015) was written about. That is why it is the
recommended input rather than the plain counts.

**How should I analyse a time course?**

An LRT against a reduced model without the time term answers "does this
gene change at all over time." If you care about the shape, put a spline
basis in the design (`splines::ns()`), or use a dedicated tool such as
maSigPro or ImpulseDE2. More time points with fewer replicates each is
often a better trade than the reverse.

**I ran five contrasts. Is my FDR still 5%?**

Per contrast, yes — each `results()` call adjusts independently. Across
all five jointly, no. Pooling and adjusting together would control the
family-wise error across contrasts, but it is unusual in practice and
quite conservative. Most papers report per-contrast FDR; the honest
thing is to say in the methods that is what you did.
