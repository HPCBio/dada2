# Tests for mergePairs()
#
# INPUT DATA: test_fnF / test_fnR and test_err defined in helper-data.R
# Unit tests use tperr1 (test_err) to avoid the cost of learnErrors().
# To test with real data, replace test_fnF/R and test_err in helper-data.R.

# Shared setup: dereplicate and denoise both samples (forward + reverse)
local({
  derepsF <<- derepFastq(test_fnF, verbose = FALSE)
  derepsR <<- derepFastq(test_fnR, verbose = FALSE)
  dadsF   <<- dada(derepsF, err = test_err, multithread = FALSE, verbose = FALSE)
  dadsR   <<- dada(derepsR, err = test_err, multithread = FALSE, verbose = FALSE)
})

test_that("mergePairs returns a data.frame with expected columns for one sample", {
  # INPUT DATA: derepsF[[1]], derepsR[[1]], dadsF[[1]], dadsR[[1]] — see helper-data.R
  merger <- mergePairs(
    dadaF  = dadsF[[1]], derepF = derepsF[[1]],
    dadaR  = dadsR[[1]], derepR = derepsR[[1]],
    verbose = FALSE
  )

  expect_true(is.data.frame(merger))
  expect_true(all(c("sequence", "abundance", "accept") %in% names(merger)),
              label = "required columns present")
  expect_true(nrow(merger) > 0, label = "at least one merger row")
})

test_that("mergePairs produces accepted mergers with valid sequences", {
  # INPUT DATA: derepsF[[1]], derepsR[[1]], dadsF[[1]], dadsR[[1]] — see helper-data.R
  merger <- mergePairs(
    dadaF  = dadsF[[1]], derepF = derepsF[[1]],
    dadaR  = dadsR[[1]], derepR = derepsR[[1]],
    verbose = FALSE
  )

  accepted <- merger[merger$accept, ]
  expect_true(nrow(accepted) > 0, label = "at least one accepted merger")
  expect_true(all(accepted$abundance >= 1), label = "accepted abundances >= 1")
  expect_true(all(grepl("^[ACGTacgt]+$", accepted$sequence)),
              label = "merged sequences contain only A/C/G/T")
})

test_that("mergePairs handles multiple samples when called in a loop", {
  # INPUT DATA: derepsF, derepsR, dadsF, dadsR (both samples) — see helper-data.R
  mergers <- lapply(seq_along(test_fnF), function(i) {
    mergePairs(
      dadaF  = dadsF[[i]], derepF = derepsF[[i]],
      dadaR  = dadsR[[i]], derepR = derepsR[[i]],
      verbose = FALSE
    )
  })

  expect_equal(length(mergers), length(test_fnF))
  for (i in seq_along(mergers)) {
    expect_true(is.data.frame(mergers[[i]]),
                label = sprintf("sample %d: merger is a data.frame", i))
    accepted <- mergers[[i]][mergers[[i]]$accept, ]
    expect_true(nrow(accepted) > 0,
                label = sprintf("sample %d: has accepted mergers", i))
  }
})

test_that("mergePairs minOverlap parameter filters short overlaps", {
  # With a very high minOverlap, fewer (or zero) mergers should be accepted.
  # INPUT DATA: derepsF[[1]], derepsR[[1]], dadsF[[1]], dadsR[[1]] — see helper-data.R
  merger_strict <- mergePairs(
    dadaF  = dadsF[[1]], derepF = derepsF[[1]],
    dadaR  = dadsR[[1]], derepR = derepsR[[1]],
    minOverlap = 100,
    verbose = FALSE
  )
  merger_lenient <- mergePairs(
    dadaF  = dadsF[[1]], derepF = derepsF[[1]],
    dadaR  = dadsR[[1]], derepR = derepsR[[1]],
    minOverlap = 10,
    verbose = FALSE
  )

  n_accepted_strict  <- sum(merger_strict$accept)
  n_accepted_lenient <- sum(merger_lenient$accept)
  expect_lte(n_accepted_strict, n_accepted_lenient,
             label = "stricter minOverlap yields <= accepted mergers")
})
