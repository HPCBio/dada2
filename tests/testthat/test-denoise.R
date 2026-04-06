# Tests for dada() and learnErrors()
#
# INPUT DATA:
#   - test_fnF / test_fnR defined in helper-data.R (raw FASTQs)
#   - test_err (tperr1) defined in helper-data.R
#     Unit tests use tperr1 to avoid the cost of learnErrors().
#     Replace test_err in helper-data.R with learnErrors() output for real data.

# ---------------------------------------------------------------------------
# dada() — sample inference using a pre-built error matrix
# ---------------------------------------------------------------------------

test_that("dada returns a valid dada object for a single sample", {
  # INPUT DATA: derepFastq on test_fnF[1]; error rates from test_err (tperr1)
  derep <- derepFastq(test_fnF[1], verbose = FALSE)
  dd    <- dada(derep, err = test_err, multithread = FALSE, verbose = FALSE)

  expect_s4_class(dd, "dada")
  expect_true(!is.null(dd$denoised),    label = "has $denoised")
  expect_true(!is.null(dd$clustering),  label = "has $clustering")

  # $denoised: named integer vector — names are inferred ASV sequences
  expect_true(length(dd$denoised) > 0, label = "at least one ASV inferred")
  expect_true(all(dd$denoised >= 1),   label = "all ASV abundances >= 1")
  expect_true(all(nchar(names(dd$denoised)) > 0),
              label = "ASV sequences are non-empty strings")
  expect_true(all(grepl("^[ACGTacgt]+$", names(dd$denoised))),
              label = "ASV sequences contain only A/C/G/T")
})

test_that("dada processes a list of derep objects and returns a list of dada objects", {
  # INPUT DATA: derepFastq on test_fnF; error rates from test_err (tperr1)
  dereps <- derepFastq(test_fnF, verbose = FALSE)
  dds    <- dada(dereps, err = test_err, multithread = FALSE, verbose = FALSE)

  expect_true(is.list(dds))
  expect_equal(length(dds), length(test_fnF))
  for (i in seq_along(dds)) {
    expect_s4_class(dds[[i]], "dada")
    expect_true(length(dds[[i]]$denoised) > 0,
                label = sprintf("element %d has inferred ASVs", i))
  }
})

test_that("dada infers fewer unique sequences than input unique reads", {
  # Denoising should collapse errors — ASV count must be <= unique read count.
  # INPUT DATA: derepFastq on test_fnF[1]; error rates from test_err (tperr1)
  derep <- derepFastq(test_fnF[1], verbose = FALSE)
  dd    <- dada(derep, err = test_err, multithread = FALSE, verbose = FALSE)

  expect_lte(length(dd$denoised), length(derep$uniques),
             label = "ASVs <= unique input reads")
})

# ---------------------------------------------------------------------------
# learnErrors() — error rate estimation
# ---------------------------------------------------------------------------

test_that("learnErrors returns a valid error rate matrix", {
  # INPUT DATA: raw forward FASTQs (test_fnF) — see helper-data.R
  # Note: learnErrors reads directly from FASTQ files; replace test_fnF
  # in helper-data.R to learn error rates from your own sequencing data.
  errs <- learnErrors(test_fnF, multithread = FALSE, verbose = FALSE)

  err_mat <- errs$err_out

  # Should be a 16-row numeric matrix (one row per nucleotide transition)
  expect_true(is.matrix(err_mat),    label = "error rates are a matrix")
  expect_equal(nrow(err_mat), 16,   label = "16 rows (one per nt transition)")
  expect_true(ncol(err_mat)  > 1,    label = "multiple quality score columns")

  # All values are probabilities
  expect_true(all(err_mat >= 0), label = "all error rates >= 0")
  expect_true(all(err_mat <= 1), label = "all error rates <= 1")
})

test_that("learnErrors output is usable as dada() err argument", {
  # Round-trip: learn errors from the bundled FASTQs then run dada().
  # INPUT DATA: test_fnF for both learnErrors and derepFastq — see helper-data.R
  err   <- learnErrors(test_fnF, multithread = FALSE, verbose = FALSE)
  derep <- derepFastq(test_fnF[1], verbose = FALSE)
  dd    <- dada(derep, err = err, multithread = FALSE, verbose = FALSE)

  expect_s4_class(dd, "dada")
  expect_true(length(dd$denoised) > 0, label = "ASVs inferred with learned errors")
})
