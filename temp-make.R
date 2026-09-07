# temp-make.R -- run the methylation pipeline end to end on another machine.
#
# TEMPORARY. This is a runner for one re-run, not a build system: it hard-codes
# the three scripts that have to be re-executed after the normalization change
# and does nothing else. Delete it once the results are in.
#
# WHY A RE-RUN. The pipeline moved from functional to quantile normalization.
# Functional normalization leaves a shift behind that depends on the
# methylation level of the site, and the mean of that shift is read as signal
# by the adaptive shrinkage step, which drove pi0 under the reporting threshold
# and declared every position differentially methylated. See the note at
# `gset_file` in analysis/R/methylation-case-study-dataprep.R. Everything
# downstream of the arrays therefore has to be recomputed.
#
# Each step runs in its OWN R process rather than being source()d into this
# one. The full run holds a results object of the order of a gigabyte, and a
# fresh process per step means memory is returned between steps and a failure
# in one does not poison the next.
#
# EVERY STEP IS RESUMABLE. All four scripts skip work whose output file already
# exists, so re-running this script after an interruption continues where it
# stopped rather than starting over. It is safe to kill and restart.
#
# Expect roughly 45 hours in total on 16 cores, dominated by two steps:
#
#   dataprep     ~15 min   (plus ~722 MB of downloads if the IDATs are absent)
#   permutation  ~25 h     (200 iterations at ~7.6 min each)
#   full run     ~20 h     (two REML arms: ~10.6 h beta, ~9.0 h M)
#
# The machine should be otherwise idle and should not sleep.


## Preflight ###############################################################

stopifnot(file.exists("DESCRIPTION"), file.exists("renv.lock"))

# ON THE NEW MACHINE, RUN renv::restore() FIRST. The lockfile now carries the
# whole pipeline, but it did not until recently: renv is in `explicit` snapshot
# mode (renv/settings.json), which records only what DESCRIPTION declares, and
# DESCRIPTION declared neither minfi, nor maxprobes, nor the EPIC annotation
# packages, nor ashr. They were present in the local library, so everything ran
# and `renv::status()` reported a consistent project, while a restore on a
# fresh machine would have installed none of them. They are declared now.
#
# The check below is the backstop for that class of mistake, and it is here
# rather than at the end because the failure it catches would otherwise arrive
# 45 hours in, at the last step. Note maxprobes is GitHub-only
# (markgene/maxprobes), so it cannot come from a CRAN mirror.
required <- c(
  "ashr",
  "etrunct",
  "glmmTMB",
  "here",
  "maxprobes",
  "minfi",
  "pbapply",
  "seqwrap",
  "GEOquery",
  "IlluminaHumanMethylationEPICmanifest"
)

missing <- required[!vapply(
  required,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]

if (length(missing)) {
  stop(
    "missing packages: ",
    paste(missing, collapse = ", "),
    ".\n  Run renv::restore() in this project first.",
    call. = FALSE
  )
}

# The core count is not configured anywhere: both long scripts call
# parallel::detectCores() and pass it to every seqwrap() call, so the run uses
# whatever this machine has. It is printed because it is the single largest
# determinant of how long the next two days take -- the published timings
# (~20 h for the full run, ~25 h for the permutation) were measured at 16
# cores, and they scale roughly inversely with this number.
cores <- parallel::detectCores()
if (is.na(cores) || cores < 1L) cores <- 1L

cat("cores detected:", cores, "\n")
cat("timings quoted in this script were measured at 16 cores.\n")


## The steps, in order #####################################################
#
# The permutation runs before the full run because it is what licenses the
# tests the full run reports. They are independent, so the order can be
# swapped if the full-run results are wanted sooner.
steps <- c(
  "analysis/R/methylation-case-study-dataprep.R",
  "analysis/R/methylation-error-permutation.R",
  "analysis/R/methylation-case-full.R"
)

stopifnot(all(file.exists(steps)))

rscript <- file.path(R.home("bin"), "Rscript")

log_dir <- "temp-make-logs"
dir.create(log_dir, showWarnings = FALSE)

for (s in steps) {
  log_file <- file.path(log_dir, paste0(basename(s), ".log"))

  cat("\n==== ", s, "  started ", format(Sys.time()), "\n", sep = "")
  cat("     log: ", log_file, "\n", sep = "")
  utils::flush.console()

  t0 <- Sys.time()

  # stdout and stderr both go to the log; the console keeps only the step
  # banners, so an unattended run leaves something readable behind.
  status <- system2(rscript, s, stdout = log_file, stderr = log_file)

  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (status != 0) {
    stop(
      s,
      " failed with status ",
      status,
      " after ",
      sprintf("%.1f", mins),
      " min. See ",
      log_file,
      ". Fix the cause and re-run this script -- completed work is skipped.",
      call. = FALSE
    )
  }

  cat(sprintf("==== %s  done in %.1f min\n", s, mins))
  utils::flush.console()
}

cat("\nAll steps completed.\n")


## NOT run here, and why ###################################################
#
# analysis/R/methylation-timing-experiment.R
#   Measures compute cost, not data values, so the normalization change does
#   not affect its conclusions. It reads the FUNCTIONALLY normalized gset and
#   full_v2/timing.RDS, which a fresh machine will not have; both are kept in
#   derived_data on the machine that produced timing_v1/.
#
# analysis/R/methylation-shrinkage.R
#   Removed from the pipeline; the revised paper does not report the
#   shrinkage step. The script and its results (shrinkage_v1/, shrinkage_v2/)
#   are in analysis/archive/.
