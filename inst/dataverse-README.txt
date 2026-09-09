This README file was generated on 2026-03-19 (YYYY-MM-DD) by Daniel Hammarström.
Last updated: 2026-09-09.


-------------------
GENERAL INFORMATION
-------------------
// Title of Dataset: Simulated data sets for: seqwrap: an R package for flexible iterative fitting of high-dimensional data
// DOI: https://doi.org/10.18710/I7U71O
// Contact Information
     // Name: Daniel Hammarström
     // Institution: University of Inland Norway
     // Email: daniel.hammarstrom@inn.no
     // ORCID: https://orcid.org/0000-0001-8360-2100

// Authors: Chidimma Echebiri; Ellefsen, Stian; Ahmad, Rafi; Hammarström, Daniel
// Data Type: RDS files (R objects).
// Date of Collection: 2025-10 to 2026-09 [YYYY-MM]


// Description of dataset:
The datasets in this repository contain simulated data and results from analyses of them, together with the results of a DNA methylation case study, which can be used to replicate the results in the paper "seqwrap: an R package for flexible iterative fitting of high-dimensional data". The files are organised in two groups by their directory label:

1. Simulation results (`estimates/`, `estimates2/`, `evaluations/`, `evaluations2/`, `simdata/`, `simdata2/`). These belong in the folder `analysis/data/raw_data/` of the source repository available at https://www.github.com/trainome/seqwrap-paper.

2. Methylation case-study results (everything under `derived_data/`). These belong in the folder `analysis/data/derived_data/` of the same repository, with the `derived_data/` prefix removed.

The function `download_dataverse()` in the repository (`analysis/R/simulation-functions.R`) downloads all files and places them in the correct folders. See https://www.github.com/trainome/seqwrap-paper for details on how to extract and analyze the data.

// Version history:
Version 1 (2026-03): simulation results.
Version 2 (2026-09): added the methylation case-study results under `derived_data/`.



--------------------------
METHODOLOGICAL INFORMATION
--------------------------
The purpose of this repository is to enable the reproduction of results presented in the manuscript "seqwrap: an R package for flexible iterative fitting of high-dimensional data" (see https://www.github.com/trainome/seqwrap-paper). Data simulations, the analysis of generated data sets, and the methylation case study are computationally intensive; this repository contains the generated data sets and the model estimates presented in the manuscript.

// Simulations

Data sets were simulated based on estimates from an observed RNA-seq dataset (Pillon et al. 2022) comprising two independent groups and three time points, using a within-participant repeated-measures design. Simulations used a negative binomial distribution to generate count data for 10,000 gene targets across 20 data sets (10 data sets per dispersion scenario; see the manuscript for details). For a fixed proportion of targets, the parameter values were set to have non-null effects to evaluate sensitivity and false discovery rates across statistical models. Values for dispersion parameters were drawn from the observed mean-dispersion relationship in a preliminary model of observed data.

The generated data sets are stored in `simdata` and `simdata2`:
`/clean/`, contains the filtered counts after removal of low-expression targets; `/popeffect/`, contains the underlying population effects for each target; `/raw/`, contains raw, unfiltered target counts.

Filtered target counts were analyzed using five statistical models (see the manuscript for details). The estimates from each model and data set are stored in `estimates` and `estimates2`. The naming convention follows `[model]_estimates_[dataset].RDS`, where model refers to a specific model specification, and dataset refers to a data set within a dispersion scenario. Similarly, `evaluations` and `evaluations2` contain evaluations from each model (and target) fit.

// Methylation case study

The case study re-analyses the Illumina EPIC array data of Seaborne et al. (2018), GEO accession GSE114763: eight participants sampled at five time points (baseline, 30 min after a first session, after a 7-week loading period, after an unloading period, and after a reloading period). Raw IDAT files were downloaded from GEO, quantile normalized between arrays (minfi::preprocessQuantile) and probe-filtered, leaving 770 441 CpG sites. Two target models were fitted to every site with seqwrap: a normal model on M-values and a beta regression model on beta-values, both with a participant-level random intercept and REML estimation with Satterthwaite degrees of freedom (see the manuscript for details).

`derived_data/seaborne-gset-quantile.RDS` is the normalized, filtered array set (a minfi GenomicRatioSet) and `derived_data/seaborne-metadata.RDS` the sample metadata. Both were produced by `analysis/R/methylation-case-study-dataprep.R`.

`derived_data/permutation_v4/` contains 200 iterations of a within-participant permutation study (one file per iteration) used to estimate Type I error and power of four model arms fitted with glmmTMB (beta, beta_reml, m, m_reml). Produced by `analysis/R/methylation-error-permutation.R`.

`derived_data/full_v3/` contains the two REML arms (`beta_reml`, `m_reml`) fitted to all 770 441 sites on the quantile normalized arrays, their metadata files, and the wall-clock timing of the run. Produced by `analysis/R/methylation-case-full.R`.

`derived_data/full_v2/timing.RDS` contains the wall-clock timing of an earlier four-arm run (beta, beta_reml, m, m_reml) on functionally normalized arrays. Only the timing is included; the manuscript quotes it and the timing experiment validates its extrapolation against it. The arm files of that run are superseded by `full_v3/` and are not included.

`derived_data/timing_v1/` contains the timing experiment: the four arms timed on sub-samples of the same arrays while varying the number of targets, the model complexity, and the number of samples (one file per cell, `cell-*.RDS`), the experimental design (`design.RDS`), the sub-sampling draws (`draws.RDS`), the warm-up run (`warmup.RDS`), and the summaries used in the manuscript (`timing-cells.RDS`, `timing-scaling.RDS`). Produced by `analysis/R/methylation-timing-experiment.R`.

References:

Pillon NJ, Smith JAB, Alm PS, Chibalin AV et al. Distinctive exercise-induced inflammatory response and exerkine induction in skeletal muscle of people with type 2 diabetes. Sci Adv 2022 Sep 9;8(36):eabo3192. PMID: 36070371

Seaborne RA, Strauss J, Cocks M, Shepherd S et al. Human Skeletal Muscle Possesses an Epigenetic Memory of Hypertrophy. Sci Rep 2018 Jan 30;8(1):1898. PMID: 29382913

--------------------
DATA & FILE OVERVIEW
--------------------
// Folder structure and file list:
-- estimates/
|   |-- m1_estimates_1.RDS
|   |-- …
|   |-- m5_estimates_10.RDS
-- estimates2/
|   |-- m1_estimates_1.RDS
|   |-- …
|   |-- m5_estimates_10.RDS
-- evaluations/
|   |-- m1_evaluations_1.RDS
|   |-- …
|   |-- m5_evaluations_10.RDS
-- evaluations2/
|   |-- m1_evaluations_1.RDS
|   |-- …
|   |-- m5_evaluations_10.RDS
-- simdata/
|   |-- clean/
|   |   |-- clean_dataset_1.RDS
|   |   |-- …
|   |   |-- clean_dataset_10.RDS
|   |-- popeffect/
|   |   |-- population_effects_1.RDS
|   |   |-- …
|   |   |-- population_effects_10.RDS
|   |-- raw/
|   |   |-- dataset_1.RDS
|   |   |-- …
|   |   |-- dataset_10.RDS
-- simdata2/
|   |-- clean/
|   |   |-- clean_dataset_1.RDS
|   |   |-- …
|   |   |-- clean_dataset_10.RDS
|   |-- popeffect/
|   |   |-- population_effects_1.RDS
|   |   |-- …
|   |   |-- population_effects_10.RDS
|   |-- raw/
|   |   |-- dataset_1.RDS
|   |   |-- …
|   |   |-- dataset_10.RDS
-- derived_data/
|   |-- seaborne-gset-quantile.RDS
|   |-- seaborne-metadata.RDS
|   |-- permutation_v4/
|   |   |-- perm_0001.RDS
|   |   |-- …
|   |   |-- perm_0200.RDS
|   |-- full_v3/
|   |   |-- beta_reml.RDS
|   |   |-- beta_reml-meta.RDS
|   |   |-- m_reml.RDS
|   |   |-- m_reml-meta.RDS
|   |   |-- timing.RDS
|   |-- full_v2/
|   |   |-- timing.RDS
|   |-- timing_v1/
|   |   |-- cell-*.RDS (152 files)
|   |   |-- design.RDS
|   |   |-- draws.RDS
|   |   |-- timing-cells.RDS
|   |   |-- timing-scaling.RDS
|   |   |-- warmup.RDS

-----------------------------------------
DATA-SPECIFIC INFORMATION
-----------------------------------------
Scripts for reproducing results in the manuscript using these data are available at https://www.github.com/trainome/seqwrap-paper.

--------------------------
SHARING/ACCESS INFORMATION
--------------------------
// Licenses/Restrictions: See Terms tab.
// Links to publications that cite or use the data: See metadata field Related Publication.
// Links/relationships to related data sets: See metadata field Related Dataset.
// Data sources: See metadata field Data Source. The methylation case study uses raw data from GEO accession GSE114763 (Seaborne et al. 2018).
// Recommended citation: See citation generated by repository.
