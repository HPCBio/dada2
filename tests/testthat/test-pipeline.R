# End-to-end paired-end pipeline integration test
#
# This test runs the complete DADA2 pipeline from raw FASTQs to a chimera-free
# sequence table. It exercises every major step in sequence and verifies
# that outputs chain correctly between steps.
#
# INPUT DATA: test_fnF / test_fnR defined in helper-data.R
#   Step 1 (filterAndTrim)  — reads test_fnF, test_fnR
#   Step 2 (derepFastq)     — reads filtered FASTQ files written in Step 1
#   Step 3 (learnErrors)    — reads filtered FASTQ files written in Step 1
#   Step 4 (dada)           — uses derep objects and learned error rates
#   Step 5 (mergePairs)     — uses dada and derep objects from Steps 2 & 4
#   Step 6 (makeSequenceTable) — uses merged pairs from Step 5
#   Step 7 (removeBimeraDenovo) — uses sequence table from Step 6
#
# To run the pipeline with your own data, update test_fnF / test_fnR in
# helper-data.R (or set DADA2_TEST_FWD_1, DADA2_TEST_REV_1, etc.).
#
# This test is marked slow and will be skipped when the environment variable
# DADA2_SKIP_SLOW=1 is set (e.g., for quick smoke-test CI runs).

test_that("full paired-end pipeline produces a non-chimeric sequence table", {
  skip_if(
    Sys.getenv("DADA2_SKIP_SLOW") == "1",
    "Skipping slow integration test (DADA2_SKIP_SLOW=1)"
  )

  n_samples <- length(test_fnF)
  sample_names <- paste0("sample", seq_len(n_samples))

  # ── Step 1: Filter and trim ───────────────────────────────────────────────
  # INPUT DATA: test_fnF (forward), test_fnR (reverse) — see helper-data.R
  filtF <- replicate(n_samples, tempfile(fileext = ".fastq.gz"))
  filtR <- replicate(n_samples, tempfile(fileext = ".fastq.gz"))
  on.exit(unlink(c(filtF, filtR)), add = TRUE)

  filter_out <- filterAndTrim(
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

  expect_true(all(filter_out[, "reads.out"] > 0),
              label = "Step 1: all samples have reads after filtering")

  # Keep only samples with surviving reads
  keep       <- filter_out[, "reads.out"] > 0
  filtF      <- filtF[keep]
  filtR      <- filtR[keep]
  sample_names <- sample_names[keep]

  # ── Step 2: Dereplicate ───────────────────────────────────────────────────
  # INPUT DATA: filtered FASTQs from Step 1
  derepsF <- derepFastq(filtF, verbose = FALSE)
  derepsR <- derepFastq(filtR, verbose = FALSE)

  if (length(filtF) == 1) {
    derepsF <- list(derepsF)
    derepsR <- list(derepsR)
  }

  expect_equal(length(derepsF), length(filtF),
               label = "Step 2: one forward derep object per sample")
  expect_equal(length(derepsR), length(filtR),
               label = "Step 2: one reverse derep object per sample")

  # ── Step 3: Learn error rates ─────────────────────────────────────────────
  # INPUT DATA: filtered FASTQs from Step 1 (learnErrors reads them internally)
  errF <- learnErrors(filtF, multithread = FALSE, verbose = FALSE)
  errR <- learnErrors(filtR, multithread = FALSE, verbose = FALSE)

  expect_equal(nrow(errF$err_out), 16L, label = "Step 3: forward error matrix has 16 rows")
  expect_equal(nrow(errR$err_out), 16L, label = "Step 3: reverse error matrix has 16 rows")

  # ── Step 4: Sample inference ──────────────────────────────────────────────
  # INPUT DATA: derep objects from Step 2; error matrices from Step 3
  dadsF <- dada(derepsF, err = errF, multithread = FALSE, verbose = FALSE)
  dadsR <- dada(derepsR, err = errR, multithread = FALSE, verbose = FALSE)

  if (!is.list(dadsF)) dadsF <- list(dadsF)
  if (!is.list(dadsR)) dadsR <- list(dadsR)

  expect_equal(length(dadsF), length(filtF),
               label = "Step 4: one forward dada object per sample")
  expect_true(all(vapply(dadsF, function(d) length(d$denoised) > 0, logical(1))),
              label = "Step 4: all samples have inferred forward ASVs")

  # ── Step 5: Merge paired reads ────────────────────────────────────────────
  # INPUT DATA: dada + derep objects from Steps 2 & 4
  mergers <- lapply(seq_along(filtF), function(i) {
    mergePairs(
      dadaF  = dadsF[[i]], derepF = derepsF[[i]],
      dadaR  = dadsR[[i]], derepR = derepsR[[i]],
      verbose = FALSE
    )
  })
  names(mergers) <- sample_names

  expect_equal(length(mergers), length(filtF),
               label = "Step 5: one merger data.frame per sample")
  expect_true(any(vapply(mergers, function(m) any(m$accept), logical(1))),
              label = "Step 5: at least one sample has accepted mergers")

  # ── Step 6: Build sequence table ─────────────────────────────────────────
  # INPUT DATA: mergers from Step 5
  seqtab <- makeSequenceTable(mergers)

  expect_true(is.matrix(seqtab),
              label = "Step 6: sequence table is a matrix")
  expect_equal(nrow(seqtab), length(filtF),
               label = "Step 6: one row per sample in sequence table")
  expect_true(ncol(seqtab) > 0,
              label = "Step 6: at least one ASV in sequence table")
  expect_true(sum(seqtab) > 0,
              label = "Step 6: sequence table has non-zero counts")

  # ── Step 7: Remove chimeras ───────────────────────────────────────────────
  # INPUT DATA: sequence table from Step 6
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose = FALSE)

  expect_true(is.matrix(seqtab.nochim),
              label = "Step 7: chimera-free table is a matrix")
  expect_equal(nrow(seqtab.nochim), nrow(seqtab),
               label = "Step 7: row count unchanged after chimera removal")
  expect_lte(ncol(seqtab.nochim), ncol(seqtab),
             label = "Step 7: chimera removal cannot add ASVs")
  expect_true(sum(seqtab.nochim) > 0,
              label = "Step 7: reads survive chimera removal")

  # ── Read tracking summary ─────────────────────────────────────────────────
  reads_in      <- filter_out[keep, "reads.in"]
  reads_filtered <- filter_out[keep, "reads.out"]
  reads_final   <- rowSums(seqtab.nochim)

  expect_true(all(reads_final <= reads_filtered),
              label = "final read counts <= filtered read counts")
  expect_true(all(reads_filtered <= reads_in),
              label = "filtered read counts <= raw read counts")
})
