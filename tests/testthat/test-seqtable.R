# Tests for makeSequenceTable() and removeBimeraDenovo()
#
# INPUT DATA: test_fnF / test_fnR and test_err defined in helper-data.R
# Unit tests use tperr1 (test_err) to avoid the cost of learnErrors().
# To test with real data, replace test_fnF/R and test_err in helper-data.R.

# Shared setup: build a minimal sequence table from both samples
local({
  derepsF <- derepFastq(test_fnF, verbose = FALSE)
  derepsR <- derepFastq(test_fnR, verbose = FALSE)
  dadsF   <- dada(derepsF, err = test_err, multithread = FALSE, verbose = FALSE)
  dadsR   <- dada(derepsR, err = test_err, multithread = FALSE, verbose = FALSE)

  mergers <<- lapply(seq_along(test_fnF), function(i) {
    mergePairs(
      dadaF  = dadsF[[i]], derepF = derepsF[[i]],
      dadaR  = dadsR[[i]], derepR = derepsR[[i]],
      verbose = FALSE
    )
  })
  names(mergers) <<- paste0("sample", seq_along(test_fnF))
})

# ---------------------------------------------------------------------------
# makeSequenceTable()
# ---------------------------------------------------------------------------

test_that("makeSequenceTable returns a matrix with correct dimensions", {
  # INPUT DATA: mergers list built from test_fnF/R + test_err — see helper-data.R
  seqtab <- makeSequenceTable(mergers)

  expect_true(is.matrix(seqtab),   label = "sequence table is a matrix")
  expect_true(is.integer(seqtab) || is.numeric(seqtab),
              label = "sequence table contains counts")
  expect_equal(nrow(seqtab), length(mergers),
               label = "one row per sample")
  expect_true(ncol(seqtab) > 0, label = "at least one ASV column")
  expect_equal(rownames(seqtab), names(mergers),
               label = "row names match sample names")
})

test_that("makeSequenceTable column names are valid DNA sequences", {
  # Column names of the sequence table are ASV sequences.
  # INPUT DATA: mergers list — see helper-data.R
  seqtab <- makeSequenceTable(mergers)

  expect_true(all(grepl("^[ACGTacgt]+$", colnames(seqtab))),
              label = "all ASV column names are valid DNA sequences")
})

test_that("makeSequenceTable abundances are non-negative integers", {
  # INPUT DATA: mergers list — see helper-data.R
  seqtab <- makeSequenceTable(mergers)

  expect_true(all(seqtab >= 0), label = "no negative counts")
  expect_true(any(seqtab > 0),  label = "at least one non-zero count")
})

test_that("makeSequenceTable total reads match sum of accepted merger abundances", {
  # The total read count in the table should equal the sum of accepted merger
  # abundances across all samples.
  # INPUT DATA: mergers list — see helper-data.R
  seqtab <- makeSequenceTable(mergers)

  expected_total <- sum(vapply(mergers, function(m) {
    sum(m$abundance[m$accept])
  }, numeric(1)))

  expect_equal(sum(seqtab), expected_total,
               label = "sequence table totals match accepted merger abundances")
})

# ---------------------------------------------------------------------------
# removeBimeraDenovo()
# ---------------------------------------------------------------------------

test_that("removeBimeraDenovo returns a matrix with same row structure", {
  # INPUT DATA: sequence table from mergers — see helper-data.R
  seqtab       <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = FALSE)

  expect_true(is.matrix(seqtab.nochim), label = "result is a matrix")
  expect_equal(nrow(seqtab.nochim), nrow(seqtab),
               label = "same number of samples after chimera removal")
  expect_equal(rownames(seqtab.nochim), rownames(seqtab),
               label = "row names preserved")
})

test_that("removeBimeraDenovo does not increase column count", {
  # INPUT DATA: sequence table from mergers — see helper-data.R
  seqtab        <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = FALSE)

  expect_lte(ncol(seqtab.nochim), ncol(seqtab),
             label = "chimera removal cannot add ASVs")
})

test_that("removeBimeraDenovo retains non-chimeric read counts", {
  # All remaining abundance should be <= the original total.
  # INPUT DATA: sequence table from mergers — see helper-data.R
  seqtab        <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = FALSE)

  expect_lte(sum(seqtab.nochim), sum(seqtab),
             label = "total reads after chimera removal <= total reads before")
  expect_true(sum(seqtab.nochim) > 0,
              label = "at least some reads survive chimera removal")
})
