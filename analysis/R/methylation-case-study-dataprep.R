# methylation-case-study-dataprep.R ####
#
#
#
# This script is shorthand for preparing the data from Seaborne et al. 2018.
# Raw data is saved together with processed data to intensity values.
# Resulting rgset (red-green set) is saved together with metadata.
#
# here::here("analysis/data/derived_data/seaborne-rgset.RDS")
# here::here("analysis/data/derived_data/seaborne-metadata.RDS")

#### Packages and settings ####################################################

library(GEOquery) # For accessing the raw data
library(minfi) # For pre-processing of raw data
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(maxprobes) # For filtering cross-reactive probes


# Load color scale
source(here::here("analysis/figures/figure-opts.R"))

# Prepare data download.
# Note for use in the present repo:
# raw_data and derived data are gitignored so safe for
# large files. If code-re-use, check the download folder before committing
# changes to git.

# The GEO ID
gse_id <- "GSE114763"

# Set download directories
raw_dir <- here::here("analysis/data/raw_data", gse_id)
idat_dir <- file.path(raw_dir, "idat")
der_dir <- here::here("analysis/data/derived_data")

dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(der_dir, recursive = TRUE, showWarnings = FALSE)


# A function for retrying connection to GEO for downloading large data.
# This was created after issues with downloading the large zipped archive
# from GEO.
retry_url <- function(expr, tries = 5, pause = 5) {
  call <- substitute(expr) # re-evaluated on each attempt
  env <- parent.frame()
  for (i in seq_len(tries)) {
    out <- try(eval(call, env), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    if (i < tries) Sys.sleep(pause * i)
  }
  stop(
    "failed after ",
    tries,
    " attempts: ",
    conditionMessage(attr(out, "condition"))
  )
}

# Series-level metadata
#
# The metadata contains file names and titles include
# time-points and participant information.
pheno_file <- file.path(raw_dir, "pheno_raw.rds")

if (!file.exists(pheno_file)) {
  soft_url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
    "?acc=",
    gse_id,
    "&targ=gsm&form=text&view=brief"
  )

  soft <- retry_url({
    res <- curl::curl_fetch_memory(soft_url)
    if (res$status_code != 200) stop("HTTP ", res$status_code)
    strsplit(rawToChar(res$content), "\r?\n")[[1]]
  })

  # one record per sample; "!Sample_<key> = <value>" lines within each
  kv <- tibble::tibble(
    record = cumsum(startsWith(soft, "^SAMPLE")),
    line = soft
  ) |>
    dplyr::filter(startsWith(line, "!")) |>
    tidyr::separate(
      line,
      into = c("key", "value"),
      sep = " = ",
      extra = "merge"
    ) |>
    dplyr::mutate(key = stringr::str_remove(key, "^!Sample_"))

  # characteristics are themselves "<key>: <value>" pairs
  characteristics <- kv |>
    dplyr::filter(key == "characteristics_ch1") |>
    tidyr::separate(
      value,
      into = c("key", "value"),
      sep = ": ",
      extra = "merge"
    ) |>
    dplyr::filter(key %in% c("gender", "age", "weight", "height"))

  # Cleaning the metadata for downstream processing.
  pheno_raw <- kv |>
    dplyr::filter(key %in% c("geo_accession", "title")) |>
    dplyr::bind_rows(characteristics) |>
    tidyr::pivot_wider(
      id_cols = record,
      names_from = key,
      values_from = value
    ) |>
    dplyr::select(-record) |>

    mutate(title_sep = title) |>
    separate(title_sep, into = c("seq_sample_id", "timeid"), sep = ":") |>
    mutate(timeid = gsub(" Muscle_", "", timeid)) |>
    mutate(timeid = gsub("Participant_", "", timeid)) |>
    mutate(timeid = gsub("7wk_", "7wk", timeid)) |>

    separate(timeid, into = c("time", "participant", "replicate"), sep = "_") |>
    # Making time-point more obvious
    mutate(
      time = case_when(
        time == "Baseline" ~ "1_baseline",
        time == "Acute" ~ "2_acute",
        time == "7wkLoading" ~ "3_loading",
        time == "7wkUnloading" ~ "4_unloading",
        time == "7wkRe-oading" ~ "5_reloading"
      )
    ) |>
    # Fixing meta-data details
    mutate(
      weight = gsub(" kg", "", weight),
      height = gsub("cm", "", height),
      weight = as.numeric(weight),
      height = as.numeric(height),
      age = as.numeric(age)
    ) |>
    # Fix a misprinted participant id
    mutate(participant = gsub("P006", "006", participant))

  saveRDS(pheno_raw, pheno_file)
}

pheno_raw <- readRDS(pheno_file)

# Double check 40 expected samples
stopifnot(nrow(pheno_raw) == 40)

# Download and extract raw IDAT files.
# IDATs are distributed gzipped; minfi reads decompressed files
tar_file <- file.path(raw_dir, paste0(gse_id, "_RAW.tar"))

if (!file.exists(tar_file)) {
  # ~722 MB; increase the timeout for slow connections
  options(timeout = 60 * 60)
  retry_url(GEOquery::getGEOSuppFiles(
    gse_id,
    baseDir = dirname(raw_dir),
    makeDirectory = TRUE
  ))
}

if (length(list.files(idat_dir, pattern = "\\.idat$")) == 0) {
  utils::untar(tar_file, exdir = idat_dir)

  # Unzip idat files
  gz <- list.files(idat_dir, pattern = "\\.idat\\.gz$", full.names = TRUE)
  invisible(lapply(gz, R.utils::gunzip, overwrite = TRUE, remove = TRUE))
}


idat_stems <- list.files(idat_dir, pattern = "_Grn\\.idat$") |>
  stringr::str_remove("_Grn\\.idat$")

array_info <- tibble::tibble(stem = idat_stems) |>
  tidyr::separate(
    stem,
    into = c("geo_accession", "slide", "array"),
    sep = "_",
    remove = FALSE
  ) |>
  dplyr::mutate(Basename = file.path(idat_dir, stem))

metadata <- pheno_raw |>
  # Removing outlier sample, see "Quality Control of DNA Methylation Probes"
  # in @seaborne2018b
  filter(!grepl("SkM_Epi_Mem_1:", title), !grepl("SkM_Epi_Mem_39:", title)) |>

  tibble::as_tibble() |>
  dplyr::select(geo_accession, title, participant, gender, age, time) |>
  dplyr::mutate(
    label = stringr::str_remove(title, "^SkM_Epi_Mem_\\d+:\\s*")
  ) |>
  dplyr::left_join(array_info, by = "geo_accession")

stopifnot(!anyNA(metadata$time), !anyNA(metadata$participant))


rgset_file <- file.path(der_dir, "seaborne-rgset.RDS")
meta_file <- file.path(der_dir, "seaborne-metadata.RDS")
gset_file <- file.path(der_dir, "seaborne-gset-normalized.RDS")

## Read the arrays #########################################################
#
# Guarded: read.metharray() parses all 80 IDAT files and the resulting rgset
# is ~386 MB on disk, so re-reading and rewriting it on every source() is
# minutes of work and a large write for an identical object. The downloads
# above are guarded the same way.

if (!file.exists(rgset_file)) {
  rgset <- minfi::read.metharray(
    metadata$Basename,
    extended = TRUE,
    force = TRUE
  )

  # Check that orders correspond
  stopifnot(identical(colnames(rgset), basename(metadata$Basename)))

  colnames(rgset) <- metadata$geo_accession

  saveRDS(rgset, rgset_file)
}


## Quality control, normalization and probe filtering #######################
#
# This was previously carried out in analysis/R/methylation-case-full.R. It
# lives here instead because the normalized gset is an input to BOTH the full
# model run and the permutation study (analysis/R/methylation-error-
# permutation.R), and the permutation study has no reason to depend on a
# script whose main job is fitting models over every site.

if (!file.exists(gset_file)) {
  if (!exists("rgset")) rgset <- readRDS(rgset_file)

  # Determine detection p-values ###
  # See ?minfi::detectionP for details.
  det_p <- minfi::detectionP(rgset)

  # Quality control (detection limits)
  qc <- tibble::tibble(
    geo_accession = colnames(det_p),
    mean_det_p = colMeans(det_p),
    frac_failed = colMeans(det_p > 0.01)
  )

  # Below is the suggested workflow for confirming sex. We expect only males
  # and the function indicate that only one sex is present which should
  # suffice as confirmation.
  # predicted   <- minfi::getSex(minfi::mapToGenome(msetraw))
  # qc$pred_sex <- predicted$predictedSex

  # All samples pass quality control
  keep_samples <- qc$geo_accession[qc$mean_det_p < 0.01]

  rgset <- rgset[, keep_samples]
  det_p <- det_p[, keep_samples]

  # Drop any sample the QC removed from the metadata as well. In this data set
  # every sample passes, so this is a no-op -- but when it is not, the gset and
  # the metadata must not disagree about which samples exist. The previous
  # version of this code (in methylation-case-full.R) computed the filtered
  # metadata and then discarded it, leaving the saved metadata describing
  # samples the gset no longer contained.
  metadata <- metadata |> dplyr::filter(geo_accession %in% keep_samples)

  # Functional normalization
  gset <- minfi::preprocessFunnorm(rgset, ratioConvert = FALSE)

  # Initial number of probes
  nprobes_init <- dim(gset)[1]

  ## (i) detection
  det_p <- det_p[match(minfi::featureNames(gset), rownames(det_p)), ]
  keep <- rowSums(det_p < 0.01) == ncol(gset)
  gset <- gset[keep, ]

  # Number of probes after filtering
  nprobes_detectionfilter <- dim(gset)[1]

  ## (ii) SNP-affected probes
  gset <- minfi::dropLociWithSnps(gset, snps = c("SBE", "CpG"), maf = 0)

  nprobes_snp <- dim(gset)[1]

  ## (iii) cross-reactive probes
  gset <- maxprobes::dropXreactiveLoci(gset)

  nprobes_cross <- dim(gset)[1]

  ## (iv) sex chromosomes
  anno <- minfi::getAnnotation(gset)
  gset <- gset[!anno$chr %in% c("chrX", "chrY"), ]

  nprobes_final <- dim(gset)[1]

  # Print number of probes
  cat(paste0(
    "The total number of probes in the data set after filtering: ",
    nprobes_final,
    "
"
  ))

  ## Clean up memory
  rm(rgset)
  gc()

  saveRDS(gset, gset_file)
}


## Saving data for downstream use ##########################################
#
# Written last so that it reflects any QC filtering applied above.

saveRDS(metadata, meta_file)
