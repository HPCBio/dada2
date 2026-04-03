# Tests for filterAndTrim()
#
# INPUT DATA: test_fnF / test_fnR defined in helper-data.R
# To use external FASTQ files set DADA2_TEST_FWD_1, DADA2_TEST_REV_1, etc.

test_that("filterAndTrim handles single-end reads", {
  # INPUT DATA: first forward FASTQ (test_fnF[1]) — see helper-data.R
  filt <- tempfile(fileext = ".fastq.gz")
  on.exit(unlink(filt))

  out <- filterAndTrim(
    fwd  = test_fnF[1],
    filt = filt,
    truncLen = 230,
    maxN  = 0,
    maxEE = 2,
    compress = TRUE,
    verbose  = FALSE
  )

  expect_true(is.matrix(out))
  expect_equal(colnames(out), c("reads.in", "reads.out"))
  expect_true(out[, "reads.in"]  > 0, label = "input reads > 0")
  expect_true(out[, "reads.out"] > 0, label = "some reads pass filter")
  expect_lte(out[, "reads.out"], out[, "reads.in"])
  expect_true(file.exists(filt))
})

test_that("filterAndTrim handles paired-end reads", {
  # INPUT DATA: both samples forward + reverse (test_fnF, test_fnR) — see helper-data.R
  filtF <- c(tempfile(fileext = ".fastq.gz"), tempfile(fileext = ".fastq.gz"))
  filtR <- c(tempfile(fileext = ".fastq.gz"), tempfile(fileext = ".fastq.gz"))
  on.exit(unlink(c(filtF, filtR)))

  out <- filterAndTrim(
    fwd      = test_fnF,
    filt     = filtF,
    rev      = test_fnR,
    filt.rev = filtR,
    truncLen = c(230, 200),
    maxN  = 0,
    maxEE = c(2, 2),
    compress = TRUE,
    verbose  = FALSE
  )

  expect_true(is.matrix(out))
  expect_equal(nrow(out), length(test_fnF))
  expect_equal(colnames(out), c("reads.in", "reads.out"))
  expect_true(all(out[, "reads.in"]  > 0), label = "all samples have input reads")
  expect_true(all(out[, "reads.out"] > 0), label = "all samples pass filter")
  expect_true(all(out[, "reads.out"] <= out[, "reads.in"]))
  expect_true(all(file.exists(c(filtF, filtR))))
})

test_that("filterAndTrim maxN=0 removes reads with Ns", {
  # INPUT DATA: first forward FASTQ (test_fnF[1]) — see helper-data.R
  filt <- tempfile(fileext = ".fastq.gz")
  on.exit(unlink(filt))

  out <- filterAndTrim(
    fwd  = test_fnF[1],
    filt = filt,
    maxN = 0,
    verbose = FALSE
  )

  # All passing reads should have zero Ns — verify via derepFastq unique seqs
  if (out[, "reads.out"] > 0) {
    derep <- derepFastq(filt)
    seqs  <- names(derep$uniques)
    expect_false(any(grepl("N", seqs, ignore.case = TRUE)),
                 label = "no N bases in filtered reads")
  }
})

test_that("filterAndTrim truncLen clips reads to expected length", {
  # INPUT DATA: first forward FASTQ (test_fnF[1]) — see helper-data.R
  trunc_len <- 200L
  filt <- tempfile(fileext = ".fastq.gz")
  on.exit(unlink(filt))

  out <- filterAndTrim(
    fwd      = test_fnF[1],
    filt     = filt,
    truncLen = trunc_len,
    verbose  = FALSE
  )

  if (out[, "reads.out"] > 0) {
    derep <- derepFastq(filt)
    read_lens <- nchar(names(derep$uniques))
    expect_true(all(read_lens == trunc_len),
                label = "all reads truncated to truncLen")
  }
})
