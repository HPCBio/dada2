# Tests for derepFastq()
#
# INPUT DATA: test_fnF defined in helper-data.R
# To use external FASTQ files set DADA2_TEST_FWD_1 / DADA2_TEST_FWD_2.

test_that("derepFastq returns a valid derep object for a single file", {
  # INPUT DATA: first forward FASTQ (test_fnF[1]) — see helper-data.R
  derep <- derepFastq(test_fnF[1], verbose = FALSE)

  expect_s4_class(derep, "derep")
  expect_true(!is.null(derep$uniques),  label = "has $uniques")
  expect_true(!is.null(derep$quals),    label = "has $quals")
  expect_true(!is.null(derep$map),      label = "has $map")

  # $uniques: named integer vector — names are sequences, values are abundances
  expect_true(is.integer(derep$uniques) || is.numeric(derep$uniques))
  expect_true(length(derep$uniques) > 0, label = "at least one unique sequence")
  expect_true(all(derep$uniques >= 1),   label = "all abundances >= 1")
  expect_true(all(nchar(names(derep$uniques)) > 0), label = "sequence names non-empty")

  # $map: integer vector mapping each read to a unique sequence index
  expect_true(is.integer(derep$map) || is.numeric(derep$map))
  expect_equal(sum(derep$uniques), length(derep$map),
               label = "total abundance equals number of reads")
})

test_that("derepFastq processes multiple files and returns a list", {
  # INPUT DATA: both forward FASTQs (test_fnF) — see helper-data.R
  dereps <- derepFastq(test_fnF, verbose = FALSE)

  expect_true(is.list(dereps))
  expect_equal(length(dereps), length(test_fnF))
  for (i in seq_along(dereps)) {
    expect_s4_class(dereps[[i]], "derep",
                    label = sprintf("element %d is a derep object", i))
    expect_true(length(dereps[[i]]$uniques) > 0,
                label = sprintf("element %d has unique sequences", i))
  }
})

test_that("derepFastq unique sequences contain only valid DNA bases", {
  # INPUT DATA: first forward FASTQ (test_fnF[1]) — see helper-data.R
  derep <- derepFastq(test_fnF[1], verbose = FALSE)
  seqs  <- names(derep$uniques)

  expect_true(all(grepl("^[ACGTacgt]+$", seqs)),
              label = "all unique sequences contain only A/C/G/T")
})

test_that("derepFastq quality matrix dimensions match unique sequences", {
  # INPUT DATA: first forward FASTQ (test_fnF[1]) — see helper-data.R
  derep <- derepFastq(test_fnF[1], verbose = FALSE)

  expect_equal(nrow(derep$quals), length(derep$uniques),
               label = "one quality row per unique sequence")
  expect_equal(ncol(derep$quals), nchar(names(derep$uniques)[1]),
               label = "quality columns match read length")
})
