# Can the published analysis be reproduced from the deposited arrays?
#
# @seaborne2018 report SWAN-normalized beta values analysed with an ANOVA, and
# their supplementary tables (analysis/data/resource_data/seaborne2018-suppl.xlsx,
# sheets 2B-2D) list the positions called differentially methylated for the
# three contrasts against baseline. This script rebuilds that procedure on the
# arrays deposited with the paper and compares the result, position by
# position, with the published lists.
#
# The reproduction is worth doing for a specific reason. The case study argues
# that the published lists come from a test that ignores the repeated-measures
# structure -- a two-group comparison in which the five arrays contributed by a
# participant are treated as five independent samples. That argument currently
# rests on a degrees-of-freedom check: the published p-values are reproduced by
# an F on 1 and (n_baseline + n_condition - 2) df. This script asks the
# stronger question. Not "is the reported p consistent with that test", but
# "does that test, run on these arrays, return these numbers".
#
# What is being varied, and why each is here:
#
#   swan_all    SWAN, every probe on the array. The normalization the paper
#               states, with no filtering of our own -- the closest thing to
#               the published pipeline, and the set against which probe
#               COVERAGE should be read. A supplementary position missing here
#               is a genuine disagreement about the array; a position missing
#               only from the filtered sets was dropped by us, on purpose.
#
#   swan_filt   SWAN, filtered exactly as the case study filters (detection p,
#               SNP-overlapping, cross-reactive, sex chromosomes). This is the
#               comparison that says whether the filtering, rather than the
#               normalization or the test, moves the answer.
#
#   funnorm     Functional normalization, the same filtering. The case study
#               is computed on this. Its role here is to separate a
#               normalization effect from a test effect: if the published
#               numbers reproduce under SWAN and not under funnorm, the
#               difference is the normalization; if they reproduce under both,
#               the normalization is not what the published lists turn on.
#
# The test itself is the two-group comparison, unpaired: baseline against the
# later time point, participants treated as independent. With two groups the
# one-way ANOVA F on (1, n1 + n2 - 2) and the equal-variance t test are the
# same test, and the supplement reports the F, so both are computed and the F
# is what is compared against the published column.
#
# NOTE ON WHAT THIS SCRIPT DOES NOT SETTLE. Agreement here is evidence that
# the published pipeline has been identified correctly. Disagreement is not
# by itself evidence of an error in the original: the paper's pipeline ran in
# different software, on a manifest and with defaults that are not fully
# recoverable from the text, and beta values differ in the third decimal
# between implementations. The quantity to read is therefore the SHAPE of the
# agreement -- whether the same positions are called, in the same direction,
# at the same order of magnitude -- and not whether p-values agree to the
# digit.

## Packages and settings ###################################################
library(minfi)
library(dplyr)
library(ggplot2)

# preprocessSWAN() samples probes internally. The same seed as
# analysis/R/methylation-normalization-check.R, so the SWAN set built there
# and the one built here are the same set.
set.seed(2)

out_dir <- here::here("analysis/data/derived_data/published_anova")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

der_dir <- here::here("analysis/data/derived_data")
norm_dir <- file.path(der_dir, "normalization_check")

sup_file <- here::here("analysis/data/resource_data/seaborne2018-suppl.xlsx")

stopifnot(file.exists(sup_file))


## The published lists #####################################################
#
# Read exactly as in the case study (case-study-methylation.qmd,
# `seaborne-supplement`): sheet 2A is a summary, 2B-2D are the per-contrast
# lists, the header is on the second row, and only the first 13 columns carry
# data.
sup_sheets <- c(
  time3_loading = "Supp 2B; Loading",
  time4_unloading = "Supp 2C; Unloading",
  time5_reloading = "Supp 2D; Reloading"
)

sup_levels <- c(
  time3_loading = "3_loading",
  time4_unloading = "4_unloading",
  time5_reloading = "5_reloading"
)

read_sup <- function(sheet) {
  x <- suppressMessages(
    readxl::read_excel(sup_file, sheet = sheet, skip = 1, col_names = FALSE)
  )[, 1:13]
  names(x) <- c(
    "column", "id", "gene", "p_time", "p_pair", "meanratio",
    "fc", "fc_desc", "estimate", "f_time", "ss_time", "ss_error", "df_time"
  )
  dplyr::mutate(x, dplyr::across(
    c(column, p_time, p_pair, meanratio, fc, estimate, f_time, ss_time,
      ss_error, df_time),
    as.numeric
  ))
}

sup <- bind_rows(lapply(sup_sheets, read_sup), .id = "term") |>
  # The two group means the row was computed from, recovered from the two
  # quantities that determine them:
  #
  #   Estimate  = b_base - b_cond
  #   MeanRatio = b_base / b_cond
  #
  # so b_cond = Estimate / (MeanRatio - 1) and b_base = MeanRatio * b_cond.
  # These are the published beta values themselves, not a derived statistic,
  # and they are what makes a value-level comparison possible at all.
  mutate(
    b_cond = estimate / (meanratio - 1),
    b_base = b_cond * meanratio
  )

# The direction convention. `Estimate` is the baseline mean MINUS the
# condition mean on the beta scale, so a POSITIVE estimate is a position that
# is hypomethylated at the later time point. Everything below compares our
# `diff` (condition minus baseline) against `-estimate`.
stopifnot(
  all(grepl("^Baseline up", sup$fc_desc[sup$estimate > 0])),
  all(grepl("^Baseline down", sup$fc_desc[sup$estimate < 0])),
  identical(sup$p_time, sup$p_pair)
)

# The published summary sheet, read as the counts to be reproduced rather than
# recomputed from the per-contrast sheets.
sup_summary <- suppressMessages(
  readxl::read_excel(sup_file, sheet = "Supp 2A; Summary", col_names = FALSE)
)
names(sup_summary) <- c("contrast", "hyper", "hypo")
sup_summary <- sup_summary[-1, ] |>
  mutate(
    term = names(sup_sheets),
    hyper = as.numeric(hyper),
    hypo = as.numeric(hypo)
  )

metadata <- readRDS(file.path(der_dir, "seaborne-metadata.RDS"))


## The published test ######################################################
#
# Numerator df 1 (`df_time` is 1 in every row) and denominator df from the
# arrays actually deposited for the contrast. This is the check the case study
# makes; it is repeated here because everything below is built on it.
sup_ddf <- vapply(
  sup_levels,
  function(lv) {
    sum(metadata$time == "1_baseline") + sum(metadata$time == lv) - 2
  },
  numeric(1)
)

# The sum of squares for a two-group comparison is determined by the group
# difference alone: SS_time = (n1 n2 / (n1 + n2)) * diff^2. Checking the
# published SS against the published Estimate is a second, independent
# confirmation that the rows describe an unpaired two-group ANOVA on n1 + n2
# arrays -- and, because it also fixes the group sizes, that the arrays it was
# computed on are the ones deposited.
sup_n <- vapply(
  sup_levels,
  function(lv) {
    c(sum(metadata$time == "1_baseline"), sum(metadata$time == lv))
  },
  numeric(2)
)

published_test <- sup |>
  mutate(
    ddf = sup_ddf[term],
    ss_pred = (sup_n[1, term] * sup_n[2, term] /
      (sup_n[1, term] + sup_n[2, term])) * estimate^2
  ) |>
  summarise(
    .by = term,
    n = n(),
    df_time = unique(df_time),
    ddf = first(ddf),
    max_dev = max(abs(pf(f_time, 1, ddf, lower.tail = FALSE) - p_pair)),
    # Relative, because SS spans orders of magnitude across positions.
    max_ss_dev = max(abs(ss_pred - ss_time) / ss_time),
    max_p = max(p_pair)
  )

print(as.data.frame(published_test))


## Beta values #############################################################
#
# swan_all is built here and cached; the other two are read from what the
# case-study pipeline and the normalization check have already written.

swan_all_file <- file.path(out_dir, "beta-swan-all.RDS")

if (!file.exists(swan_all_file)) {
  message("building unfiltered SWAN beta values (minutes)")

  rgset <- readRDS(file.path(der_dir, "seaborne-rgset.RDS"))
  gset_fun <- readRDS(file.path(der_dir, "seaborne-gset-normalized.RDS"))
  rgset <- rgset[, colnames(gset_fun)]
  rm(gset_fun)
  gc()

  # No detection-p, SNP, cross-reactive or sex-chromosome filtering: the point
  # of this set is to be everything the array measures.
  gset_swan_all <- minfi::mapToGenome(minfi::preprocessSWAN(rgset))
  beta_swan_all <- minfi::getBeta(gset_swan_all, offset = 100)

  saveRDS(beta_swan_all, swan_all_file)

  rm(rgset, gset_swan_all)
  gc()
} else {
  beta_swan_all <- readRDS(swan_all_file)
}

gset_swan <- readRDS(file.path(norm_dir, "gset-swan.RDS"))
beta_swan_filt <- minfi::getBeta(gset_swan, offset = 100)
rm(gset_swan)
gc()

gset_fun <- readRDS(file.path(der_dir, "seaborne-gset-normalized.RDS"))
beta_funnorm <- minfi::getBeta(gset_fun, offset = 100)
rm(gset_fun)
gc()

beta_sets <- list(
  swan_all = beta_swan_all,
  swan_filt = beta_swan_filt,
  funnorm = beta_funnorm
)

# All three must carry the same arrays in the same order, or the group indices
# built once below would address different samples in different sets.
beta_sets <- lapply(beta_sets, function(m) m[, metadata$geo_accession])

message(
  "positions: ",
  paste(sprintf("%s %d", names(beta_sets), vapply(beta_sets, nrow, integer(1))),
        collapse = ", ")
)


## The test ################################################################
#
# The published test, vectorised over positions. Identical in form to
# two_group() in analysis/R/methylation-normalization-check.R, with the F
# statistic added so it can be compared against the supplement's own column.
two_group <- function(mat, i1, i2) {
  n1 <- length(i1)
  n2 <- length(i2)

  m1 <- rowMeans(mat[, i1, drop = FALSE])
  m2 <- rowMeans(mat[, i2, drop = FALSE])

  ss <- rowSums((mat[, i1, drop = FALSE] - m1)^2) +
    rowSums((mat[, i2, drop = FALSE] - m2)^2)

  ddf <- n1 + n2 - 2
  tstat <- (m1 - m2) / sqrt((ss / ddf) * (1 / n1 + 1 / n2))

  data.frame(
    id = rownames(mat),
    # The group means are kept, not just their difference: the supplement
    # reports both a difference and a ratio, which together determine the two
    # means it was computed from (see `sup` above), and comparing those means
    # against these is what separates a disagreement about the TEST from a
    # disagreement about the VALUES.
    base_mean = m1,
    cond_mean = m2,
    # condition minus baseline, i.e. -1 times the supplement's `Estimate`
    diff = m2 - m1,
    f = tstat^2,
    ddf = ddf,
    p = 2 * pt(-abs(tstat), ddf),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

tm <- metadata$time[match(colnames(beta_sets[[1]]), metadata$geo_accession)]
i_base <- which(tm == "1_baseline")

fits <- bind_rows(lapply(names(beta_sets), function(nm) {
  bind_rows(lapply(names(sup_levels), function(t) {
    r <- two_group(beta_sets[[nm]], i_base, which(tm == sup_levels[[t]]))
    r$term <- t
    r$set <- nm
    r
  }))
}))

saveRDS(fits, file.path(out_dir, "twogroup-fits.RDS"))

rm(beta_sets, beta_swan_all, beta_swan_filt, beta_funnorm)
gc()


## Comparison 1: coverage ##################################################
#
# How many of the published positions are present at all. Read on swan_all
# this is about the array; read on the filtered sets it is about our filtering,
# and a position we dropped cannot agree or disagree -- it is simply not
# comparable, which is why every agreement statistic below is computed on the
# matched positions only.
sup_ids <- distinct(sup, term, id)

# A position listed twice within one contrast would inflate every count below
# and would silently turn the joins into many-to-many.
stopifnot(nrow(sup_ids) == nrow(sup))

coverage <- fits |>
  distinct(set, term, id) |>
  inner_join(sup_ids, by = c("term", "id")) |>
  count(set, term, name = "n_matched") |>
  right_join(
    expand.grid(
      set = unique(fits$set),
      term = names(sup_sheets),
      stringsAsFactors = FALSE
    ),
    by = c("set", "term")
  ) |>
  left_join(count(sup, term, name = "n_sup"), by = "term") |>
  mutate(
    n_matched = tidyr::replace_na(n_matched, 0L),
    n_missing = n_sup - n_matched,
    pct_matched = 100 * n_matched / n_sup
  ) |>
  arrange(set, term)


## Comparison 2: agreement on the published positions ######################
#
# Joined on position and contrast. `p_sup`/`f_sup`/`est_sup` are the published
# columns, `p`/`f`/`diff` ours.
#
# The p-values are compared on the log10 scale, since they span orders of
# magnitude and the question is whether the same test was run, not whether the
# fourth digit matches. `recovered` is the practical reading: of the positions
# the paper lists, how many does this reproduction also call at the threshold
# the paper applied.
joined <- fits |>
  inner_join(
    select(sup, term, id, p_sup = p_pair, f_sup = f_time, est_sup = estimate),
    by = c("term", "id")
  ) |>
  # The supplement's estimate is baseline minus condition; ours is the
  # reverse. Flipped here so the two are on one convention.
  mutate(est_sup = -est_sup)

sup_alpha <- max(sup$p_pair)

agreement <- joined |>
  summarise(
    .by = c(set, term),
    n = n(),
    # p-values
    cor_logp = cor(log10(p), log10(p_sup)),
    med_abs_logp = median(abs(log10(p) - log10(p_sup))),
    recovered = mean(p < sup_alpha),
    # F statistics, compared directly against the published column
    cor_f = cor(f, f_sup),
    # effect sizes on the beta scale
    cor_est = cor(diff, est_sup),
    med_abs_est = median(abs(diff - est_sup)),
    sign_agree = mean(sign(diff) == sign(est_sup))
  ) |>
  arrange(set, term)


## Comparison 3: the counts the paper reports ##############################
#
# The reverse direction of comparison 2. Running the published test on these
# arrays, at the threshold the published lists were cut at, how many positions
# come out -- and are they split between hyper- and hypomethylated the way the
# summary sheet reports?
#
# `hyper` is condition above baseline (our diff > 0), which is the sheet's
# "Compared to Base Hypermethylated".
called <- fits |>
  summarise(
    .by = c(set, term),
    n_tested = n(),
    n_called = sum(p < sup_alpha, na.rm = TRUE),
    hyper = sum(p < sup_alpha & diff > 0, na.rm = TRUE),
    hypo = sum(p < sup_alpha & diff < 0, na.rm = TRUE)
  ) |>
  left_join(
    select(sup_summary, term, sup_hyper = hyper, sup_hypo = hypo),
    by = "term"
  ) |>
  mutate(
    sup_total = sup_hyper + sup_hypo,
    ratio = n_called / sup_total
  ) |>
  arrange(set, term)


## Comparison 4: is the disagreement about the test, or about the values? ###
#
# Comparisons 1-3 can only say whether the same positions come out. If they do
# not, that could still be a difference of test, of normalization, of
# filtering, or of which arrays went into which group. This section removes all
# of those at once by comparing the BETA VALUES rather than the test applied to
# them: `b_base` recovered from the supplement against `base_mean` computed
# here, position by position.
#
# The comparison needs a benchmark, because a correlation is only interpretable
# against the value it would take if everything were right. `within_ours` is
# that benchmark: the mean beta of a position at baseline against its mean beta
# at the later time point, in our data. Those are different arrays, different
# participants at two of the five time points, and a real biological contrast
# -- and they still agree at r > 0.99, because what a position's mean beta is
# mostly reflects which position it is. That is the scale on which the
# published-versus-ours correlation should be read.

value_check <- fits |>
  inner_join(select(sup, term, id, b_base, b_cond), by = c("term", "id")) |>
  filter(is.finite(b_base), is.finite(b_cond)) |>
  summarise(
    .by = c(set, term),
    n = n(),
    # the benchmark, computed on exactly these positions
    within_ours = cor(base_mean, cond_mean),
    # the same quantity, published against reproduced
    cor_base = cor(b_base, base_mean),
    cor_cond = cor(b_cond, cond_mean),
    med_abs_base = median(abs(b_base - base_mean))
  ) |>
  arrange(set, term)


# The structural check. Mean methylation is not arbitrary across the array: CpG
# islands are hypomethylated and open-sea positions are methylated, in every
# tissue, in every data set. A set of beta values that is correctly matched to
# its probe identifiers reproduces that pattern; one that is not is flat,
# because it is then a random position's value under each label.
#
# Selection is the confound to control here, and it is a real one: the
# published rows are the positions that reached p < 0.05, and selecting for a
# detectable difference favours intermediate beta values, which compresses the
# pattern on its own. `ours_selected` applies the SAME kind of selection to our
# data -- the positions this reproduction calls, at the same threshold, in the
# same contrast -- so the two columns are subject to the same compression and
# can be read against each other.
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

anno_all <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
island <- setNames(anno_all$Relation_to_Island, rownames(anno_all))
gene_anno <- setNames(toupper(anno_all$UCSC_RefGene_Name), rownames(anno_all))

island_ref <- fits |>
  filter(set == "swan_all", term == "time3_loading") |>
  transmute(
    id,
    island = island[id],
    base_mean,
    selected = p < sup_alpha
  ) |>
  filter(!is.na(island))

island_check <- sup |>
  filter(term == "time3_loading", is.finite(b_base)) |>
  transmute(id, island = island[id], theirs = b_base) |>
  filter(!is.na(island)) |>
  summarise(.by = island, n_theirs = n(), theirs = mean(theirs)) |>
  full_join(
    island_ref |>
      summarise(.by = island, n_ours_all = n(), ours_all = mean(base_mean)),
    by = "island"
  ) |>
  full_join(
    island_ref |>
      filter(selected) |>
      summarise(
        .by = island,
        n_ours_sel = n(),
        ours_selected = mean(base_mean)
      ),
    by = "island"
  ) |>
  arrange(island)


# Whether the supplement is self-consistent, which decides what a mismatch can
# and cannot be attributed to.
#
# `pair_stable`  is each probe identifier paired with one and only one of the
#                file's own row identifiers.
# `cross_sheet`  the baseline group is the SAME arrays in all three contrasts,
#                so a position listed in two sheets must carry the same
#                recovered baseline mean in both. This is computed from the
#                file alone and involves none of our data.
# `gene_agree`   does the gene the sheet names for a position agree with the
#                gene the manifest annotates that position to.
#
# These matter because they separate two very different explanations. A
# spreadsheet-level accident -- a column sorted out of step with its
# neighbours -- would break them. An upstream mismatch, applied once, before
# the tables were written, would not: it would be carried identically into
# every sheet and would leave the file perfectly self-consistent.
sup_pairs <- distinct(sup, column, id)

cross_sheet <- sup |>
  filter(is.finite(b_base)) |>
  summarise(.by = id, n_sheets = n(), spread = max(b_base) - min(b_base)) |>
  filter(n_sheets > 1)

gene_rows <- sup |>
  filter(!is.na(gene), gene != "", !is.na(gene_anno[id])) |>
  mutate(ours = gene_anno[id])

gene_agree <- mean(mapply(
  function(a, b) {
    length(intersect(
      unique(strsplit(toupper(a), ";")[[1]]),
      unique(strsplit(b, ";")[[1]])
    )) > 0
  },
  gene_rows$gene,
  gene_rows$ours
))

internal <- data.frame(
  n_rows = nrow(sup),
  n_ids = n_distinct(sup$id),
  n_pairs = nrow(sup_pairs),
  pair_stable = nrow(sup_pairs) == n_distinct(sup$id),
  n_multi_sheet = nrow(cross_sheet),
  max_baseline_spread = if (nrow(cross_sheet)) max(cross_sheet$spread) else NA,
  gene_rows = nrow(gene_rows),
  gene_agree = gene_agree
)


## Save ####################################################################

saveRDS(
  list(
    published_test = published_test,
    coverage = coverage,
    agreement = agreement,
    called = called,
    value_check = value_check,
    island_check = island_check,
    internal = internal,
    sup_alpha = sup_alpha,
    sup_summary = sup_summary
  ),
  file.path(out_dir, "published-anova-check.RDS")
)


## Figure ##################################################################
#
# Published p-value against reproduced p-value on the positions the paper
# lists, one panel per contrast, one colour per normalization. The diagonal is
# exact agreement; the horizontal line is the threshold the published lists
# were cut at, so points below it are positions this reproduction also calls.
plot_dat <- filter(joined, set %in% c("swan_all", "funnorm"))
plot_dat <- slice_sample(plot_dat, n = min(30000L, nrow(plot_dat)))

fig <- ggplot(plot_dat, aes(-log10(p_sup), -log10(p), colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(
    yintercept = -log10(sup_alpha),
    colour = "grey60",
    linewidth = 0.3,
    linetype = 2
  ) +
  geom_point(alpha = 0.15, size = 0.5) +
  facet_wrap(~term) +
  labs(
    x = "Published -log10 p",
    y = "Reproduced -log10 p",
    colour = NULL,
    caption = paste0(
      "positions listed in the supplement; grey diagonal is exact agreement, ",
      "dashed line the published threshold"
    )
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "published-anova-check.png"),
  fig,
  width = 9,
  height = 3.5,
  dpi = 300
)

# The structural figure. Mean baseline methylation by CpG-island relation, the
# published values against ours -- both on the positions the supplement lists,
# with our own equally-selected positions as the third column to show that the
# pattern survives selection.
fig_island <- island_check |>
  select(island, theirs, ours_all, ours_selected) |>
  tidyr::pivot_longer(-island, names_to = "source", values_to = "beta") |>
  ggplot(aes(island, beta, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0.5, colour = "grey60", linewidth = 0.3) +
  labs(
    x = NULL,
    y = "Mean beta at baseline",
    fill = NULL,
    caption = paste(
      "islands are hypomethylated and open sea methylated in any array data;",
      "a flat profile is what values not matched to their probe give"
    )
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "island-structure.png"),
  fig_island,
  width = 8,
  height = 3.5,
  dpi = 300
)


## Report ##################################################################

cat("\nPublished threshold: p <", format(sup_alpha), "\n")

cat("\nCoverage of the published positions:\n")
print(as.data.frame(coverage))

cat("\nAgreement on the published positions:\n")
print(as.data.frame(agreement))

cat("\nPositions called by the reproduction vs the published summary:\n")
print(as.data.frame(called))

cat("\nBeta values, published vs reproduced (within_ours is the benchmark):\n")
print(as.data.frame(value_check))

cat("\nMean baseline beta by CpG-island relation:\n")
print(as.data.frame(island_check))

cat("\nInternal consistency of the supplement:\n")
print(as.data.frame(internal))
