# Gut Microbiome Diversity in R

![R](https://img.shields.io/badge/R-4.5.3-276DC3?logo=r&logoColor=white)
![phyloseq](https://img.shields.io/badge/phyloseq-microbiome%20analysis-2E8B57)
![vegan](https://img.shields.io/badge/vegan-community%20ecology-green)
![ggplot2](https://img.shields.io/badge/ggplot2-visualization-orange)
![renv](https://img.shields.io/badge/renv-reproducible%20environment-blue)

Exploratory analysis of gut microbiome composition, alpha diversity and beta diversity in R, using the `enterotype` dataset included in the `phyloseq` package.

This is an **educational project**. Its purpose is to practice R tools for microbiome analysis and to build a portfolio.

---

## Project Goal

Explore how gut microbiome communities differ between the enterotype labels of the dataset, while checking how much of that difference is explained by the sequencing technology (`SeqTech`).

The project covers:

* Inspection of the dataset and its metadata
* Filtering of rare taxa
* Taxonomic composition at genus level
* Alpha diversity (Shannon and Simpson)
* Beta diversity (Bray-Curtis, PCoA, PERMANOVA, dispersion tests)
* Unsupervised clustering (PAM) compared against the existing labels

---

## Scientific Background

* **Relative abundance:** the abundance table stores proportions per sample, not raw counts.
* **Alpha diversity:** diversity within one sample. Shannon combines richness and evenness; Simpson gives more weight to abundant taxa.
* **Beta diversity:** differences in composition between samples, measured here with the Bray-Curtis distance and shown in 2D with a PCoA.
* **PERMANOVA (`adonis2`):** tests whether group centroids differ more than expected by chance.
* **Dispersion (`betadisper` + `permutest`):** tests whether groups differ in how spread out their samples are. A PERMANOVA can be significant because of unequal dispersion, so both tests are run.
* **PAM clustering:** unsupervised clustering on the Bray-Curtis matrix. The number of clusters is chosen with the average silhouette width.
* **Enterotype:** a sample label stored in the dataset metadata.

---

## Project Structure

```text
gut-microbiome-diversity/
│
├── data/                        # Intermediate .rds objects (contents not tracked by Git)
│
├── renv/                        # renv project environment
│
├── results/
│   └── figures/                 # Figures (.png) saved by the scripts
│
├── scripts/
│   ├── 00_setup.R               # Loads phyloseq, vegan and ggplot2
│   └── phase_1.R                # Dataset load and exploration. Taxa and sample filtering
│   └── phase_2.R                # Composition and alpha diversity
│   └── phase_3.R                # Beta diversity
│   └── phase_4.R                # Machine Learning: Clustering
│
├── .gitignore
├── .Rprofile
├── renv.lock                    # Locked package versions
└── README.md
```

The scripts use paths relative to the project root (for example `source("scripts/00_setup.R")` and `data/ps_cl.rds`), so they should be run from that folder.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/moralesgomez-dev/gut-microbiome-diversity
cd gut-microbiome-diversity
```

### 2. Restore the environment

`renv.lock` records R 4.5.3, Bioconductor 3.22 and the package versions used in the project.

From the project root, restore the package library recorded in `renv.lock`:

```r
renv::restore()
```

### 3. Data

No download is needed. The dataset is loaded from `phyloseq` with `data(enterotype)`.

### 4. Run the scripts

Run the phases in order, from the project root. Each phase reads the objects saved by the previous one.

```bash
Rscript scripts/phase_1.R
Rscript scripts/phase_2.R
Rscript scripts/phase_3.R
Rscript scripts/phase_4.R
```

Run each command only after the previous phase has completed successfully. `Rscript` must be available on your `PATH`.

---

## Pipeline Overview

### Phase 1: Data inspection and filtering

* Loads `enterotype` and inspects the abundance table, sample metadata and taxonomy table.
* Checks missing values in the metadata and how sequencing technologies and enterotypes are distributed.
* Builds a sensitivity table of detection threshold × minimum fraction of samples, and the abundance retained per sample.
* Filters taxa with `filter_taxa()`: relative abundance > 1e-4 in at least 5% of the samples (14 of 280). Taxa go from 553 to 131.
* No sample is removed. Abundances are renormalized to proportions.
* Saves `data/enterotype_filt.rds`.

### Phase 2: Composition and alpha diversity

* Works with genera only (the dataset has a single taxonomic rank, `Genus`).
* Removes the unclassified taxon `-1`, renormalizes and saves `data/ps_cl.rds`.
* Selects the 10 genera with the highest cumulative relative abundance.
* Plots mean composition by enterotype, split by sequencing technology.
* Computes Shannon, Simpson and Observed per sample, in two versions: classified genera only (main analysis) and including `-1` (sensitivity analysis).
* Saves five figures:
  * `phase_2_composicion_generos.png`
  * `phase_2_shannon_generos_clasificados.png`
  * `phase_2_simpson_generos_clasificados.png`
  * `phase_2_shannon_incluye_taxon_menos_1.png`
  * `phase_2_simpson_incluye_taxon_menos_1.png`

### Phase 3: Beta diversity

* Computes Bray-Curtis distances and a PCoA on `ps_cl`, coloured by `SeqTech` and by `Enterotype`.
* Removes samples with missing `Enterotype` and recalculates distances.
* Runs PERMANOVA: `adonis2(dist ~ Enterotype + SeqTech, permutations = 999, by = "margin")`, with `set.seed(123)`.
* Runs `betadisper` + `permutest` (999 permutations) for `Enterotype` and for `SeqTech`.
* Saves `phase_3_pcoa_seqtech.png`, `phase_3_pcoa_enterotype.png`, `phase_3_dispersion_enterotype.png` and `phase_3_dispersion_seqtech.png`.

### Phase 4: Clustering

* Computes the average silhouette width for k = 2 to 6 with `cluster::pam` on the Bray-Curtis matrix and selects k = 3.
* Compares the clusters with `Enterotype` and `SeqTech` using contingency tables and row proportions.
* Colours the PCoA by cluster and saves `phase_4_pcoa_clustering.png`.


---

## Dataset

| Parameter | Value |
|-----------|-------|
| Source | `enterotype` dataset from the `phyloseq` package |
| Samples | 280 |
| Taxa (before filtering) | 553 |
| Taxa (after filtering) | 131 |
| Taxonomic rank | Genus only |
| Abundance values | Relative abundance (proportions) |
| Sequencing technologies | Illumina, Sanger, Pyro454 |
| Smallest enterotype group | Enterotype 2 (35 samples) |

Metadata variables used or checked in the scripts: `Enterotype`, `SeqTech`, `Nationality`, `Gender`, `Age` and `ClinicalStatus`.

The `enterotype` dataset distributed with `phyloseq` is based on:

Arumugam, M. et al. (2011). Enterotypes of the human gut microbiome. *Nature, 473*(7346), 174–180. [https://doi.org/10.1038/nature09944](https://doi.org/10.1038/nature09944).

`SeqTech` is complete, but `Enterotype` is missing for 9 of the 280 samples. Phase 3 excludes those 9 samples before PERMANOVA, leaving 271 samples for that analysis.

---

## Results

The scripts show these results. They are the author's readings of the outputs.

**Filtering (Phase 1).** The median abundance retained per sample after taxa filtering is 0.9988 (minimum 0.82). The number of detected taxa depends on the technology (about 96 for Illumina, about 60 for Sanger and about 37 for Pyro454), not on sample quality.

**Composition (Phase 2).** The unclassified taxon `-1` is 0 in Pyro454 samples, but about half of the sample on average in Illumina and Sanger. This is why the main analysis removes it and renormalizes.

**Beta diversity (Phase 3).** The PCoA groups samples differently by `SeqTech` and by `Enterotype`. In the PERMANOVA both variables are associated with differences in Bray-Curtis composition, and sequencing technology has a slightly larger marginal effect in this model.

The PERMANOVA used 271 samples after excluding the 9 samples with missing `Enterotype`. With marginal effects (`by = "margin"`), both terms were significant: `Enterotype` (R² = 0.238, pseudo-F = 63.72, p = 0.001) and `SeqTech` (R² = 0.274, pseudo-F = 73.37, p = 0.001). The reported p-values are the minimum resolution available with 999 permutations.

Dispersion also differed among enterotypes (PERMDISP: F = 66.93, p = 0.001) and sequencing technologies (F = 8.37, p = 0.002). Therefore, interpret the PERMANOVA cautiously: unequal within-group dispersion may contribute to its significant results, and the tests do not establish that only group centroids differ.

**Clustering (Phase 4).** k = 3 has the best average silhouette, although it is close to k = 4. The clusters partly match the enterotypes, mainly cluster 3 and enterotype 2, but there is no clear one-to-one match. Cluster 2 is clearly dominated by Pyro454, while clusters 1 and 3 are spread across technologies. This suggests that Bray-Curtis clustering captures sequencing technology.

The average silhouette widths were 0.205 (`k = 2`), 0.230 (`k = 3`), 0.224 (`k = 4`), 0.161 (`k = 5`) and 0.148 (`k = 6`). Although `k = 3` was highest, its value is low and close to `k = 4`, indicating weak cluster separation.

For `k = 3`, the cluster-by-enterotype counts were:

| Cluster | Enterotype 1 | Enterotype 2 | Enterotype 3 | Missing |
|---------|-------------:|-------------:|-------------:|--------:|
| 1 | 96 | 0 | 15 | 7 |
| 2 | 72 | 9 | 38 | 1 |
| 3 | 0 | 26 | 15 | 1 |

Within-cluster proportions show that cluster 1 is mostly enterotype 1 (81%), cluster 2 is mixed but mostly enterotype 1 (60%), and cluster 3 is mostly enterotype 2 (62%) with a substantial enterotype 3 share (36%).

The cluster-by-technology counts were:

| Cluster | Illumina | Pyro454 | Sanger |
|---------|---------:|--------:|-------:|
| 1 | 63 | 25 | 30 |
| 2 | 0 | 118 | 2 |
| 3 | 22 | 11 | 9 |

Cluster 2 is almost entirely Pyro454 (98%), while clusters 1 and 3 contain a mix of technologies. These contingency tables are descriptive; they do not test statistical association.

![PCoA by sequencing technology](results/figures/phase_3_pcoa_seqtech.png)

![PCoA by enterotype](results/figures/phase_3_pcoa_enterotype.png)

![PCoA by PAM cluster](results/figures/phase_4_pcoa_clustering.png)

---

## Limitations

### Documented in the scripts

* **Demographic and clinical variables:** `Nationality`, `Gender` and `ClinicalStatus` are available only for the 39 Sanger samples, and `Age` has many missing values. They are not used in the main analysis.
* **Technology confounding:** the number of detected taxa varies with the sequencing technology, so richness comparisons are confounded by the technique.
* **Unbalanced design:** the three enterotypes appear in all three technologies, but not in equal proportions. In Sanger, enterotype 3 dominates and enterotype 1 is a minority. Comparisons by enterotype must account for `SeqTech`.
* **Filter threshold:** the amount of abundance retained measures how much abundance is lost. It does not prove that the removed taxa are noise.
* **Samples:** no sample is removed. The minimum retention (0.82) shows no clear separation from the rest.
* **Taxonomic resolution:** only the genus rank is available.
* **Unclassified taxon `-1`:** it is absent in Pyro454 and about half of the sample in Illumina and Sanger. Each technology seems to have been processed with a different classification method. It is treated as an aggregate category of uncertain identity.
* **Relative abundances:** Shannon and Simpson can be computed from proportions, but Chao1 and ACE are not interpreted as reliable richness estimates.
* **Diversity scope:** Observed, Shannon and Simpson describe only the genera kept after the Phase 1 filter, not the original richness. `estimate_richness` may warn that there are no singletons.
* **PERMANOVA:** it should be interpreted with caution if dispersion is unequal. A significant `permutest` shows unequal dispersion between groups. It does not say which groups differ or prove that centroids differ.
* **Clustering:** the clusters do not cleanly reproduce the enterotypes and seem influenced by sequencing technology. This alone does not prove that the associations are significant.

### Additional limitations

* **Descriptive alpha diversity:** Phase 2 only draws boxplots. No statistical test is applied to alpha diversity.
* **No pairwise tests:** PERMANOVA and the dispersion tests are global. There are no post-hoc comparisons between pairs of groups.
* **Possible batch effects:** the materials do not say whether `SeqTech` is also confounded with study or cohort. If it is, technology cannot be separated from those effects.
* **Missing `Enterotype` values:** the PERMANOVA and dispersion tests exclude the 9 of 280 samples (3.2%) with missing `Enterotype`.
* **Filter sensitivity:** the filter thresholds (1e-4 and 5%) were chosen by the author. Their effect on the beta diversity and clustering results was not evaluated. Only the inclusion of taxon `-1` has a sensitivity analysis, and only for alpha diversity.
* **Choice of k:** the silhouette values for k = 3 and k = 4 are close, so the choice of k = 3 is not strongly supported.
* **Single distance:** only Bray-Curtis is used. Other distances could give different groupings.
* **Reproducibility:** only the PERMANOVA step sets a random seed (`set.seed(123)`). The permutation tests are reproducible if Phase 3 is run from top to bottom.
* **Scope:** the analysis is exploratory and educational. There is no external validation, and the results do not support clinical or causal conclusions.

---

## Author

[moralesgomez-dev](https://github.com/moralesgomez-dev)

---

## License

[FALTA: licencia; no hay archivo LICENSE en el árbol de archivos]