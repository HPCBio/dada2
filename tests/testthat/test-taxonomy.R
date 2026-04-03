# Tests for assignTaxonomy() and assignSpecies()
#
# INPUT DATA: test_train_fa, test_species_fa defined in helper-data.R
# Sequences used as query are taken from the example_seqs.fa bundled file, which
# was designed to pair with the bundled reference databases.
#
# To use external reference databases, set:
#   DADA2_TEST_TRAIN_FA   — path to your assignTaxonomy reference FASTA
#   DADA2_TEST_SPECIES_FA — path to your assignSpecies reference FASTA
# And supply matching query sequences (e.g., ASVs from your own pipeline run).

# Load bundled example query sequences
.example_seqs_path <- system.file("extdata", "example_seqs.fa", package = "dada2")
# INPUT DATA: example_seqs.fa — replace with your own ASV sequences if using
#             external reference databases
test_seqs <- getSequences(.example_seqs_path)

# ---------------------------------------------------------------------------
# assignTaxonomy()
# ---------------------------------------------------------------------------

test_that("assignTaxonomy returns a character matrix with correct structure", {
  # INPUT DATA: test_seqs (query) + test_train_fa (reference) — see helper-data.R
  taxa <- assignTaxonomy(
    seqs     = test_seqs,
    refFasta = test_train_fa,
    minBoot  = 50,
    verbose  = FALSE
  )

  expect_true(is.matrix(taxa),         label = "taxonomy result is a matrix")
  expect_true(is.character(taxa),      label = "taxonomy result contains characters")
  expect_equal(nrow(taxa), length(test_seqs),
               label = "one row per query sequence")
  expect_true(ncol(taxa) >= 2,         label = "at least two taxonomic rank columns")
})

test_that("assignTaxonomy row names are the query sequences", {
  # INPUT DATA: test_seqs + test_train_fa — see helper-data.R
  taxa <- assignTaxonomy(
    seqs     = test_seqs,
    refFasta = test_train_fa,
    minBoot  = 50,
    verbose  = FALSE
  )

  expect_equal(rownames(taxa), test_seqs,
               label = "row names are the input sequences")
})

test_that("assignTaxonomy assigns at least one sequence to a known taxon", {
  # The bundled example_seqs.fa contains sequences that match the bundled
  # training set, so at least the highest rank should be assigned.
  # INPUT DATA: test_seqs + test_train_fa — see helper-data.R
  taxa <- assignTaxonomy(
    seqs     = test_seqs,
    refFasta = test_train_fa,
    minBoot  = 0,           # low threshold to ensure something is assigned
    verbose  = FALSE
  )

  expect_true(any(!is.na(taxa[, 1])),
              label = "at least one assignment at the highest rank")
})

test_that("assignTaxonomy outputBootstraps returns list with taxa and boot matrices", {
  # INPUT DATA: test_seqs + test_train_fa — see helper-data.R
  result <- assignTaxonomy(
    seqs             = test_seqs,
    refFasta         = test_train_fa,
    minBoot          = 50,
    outputBootstraps = TRUE,
    verbose          = FALSE
  )

  expect_true(is.list(result),             label = "result is a list when outputBootstraps=TRUE")
  expect_true("tax"  %in% names(result),   label = "result has $tax")
  expect_true("boot" %in% names(result),   label = "result has $boot")
  expect_true(is.matrix(result$tax),       label = "$tax is a matrix")
  expect_true(is.matrix(result$boot),      label = "$boot is a matrix")
  expect_equal(dim(result$tax), dim(result$boot),
               label = "taxa and bootstrap matrices have same dimensions")
  expect_true(all(result$boot[!is.na(result$boot)] >= 0 &
                  result$boot[!is.na(result$boot)] <= 100),
              label = "bootstrap values are percentages (0–100)")
})

# ---------------------------------------------------------------------------
# assignSpecies()
# ---------------------------------------------------------------------------

test_that("assignSpecies returns a two-column character matrix", {
  # INPUT DATA: test_seqs (query) + test_species_fa (reference) — see helper-data.R
  species <- assignSpecies(
    seqs     = test_seqs,
    refFasta = test_species_fa,
    verbose  = FALSE
  )

  expect_true(is.matrix(species),      label = "result is a matrix")
  expect_true(is.character(species),   label = "result contains characters")
  expect_equal(nrow(species), length(test_seqs),
               label = "one row per query sequence")
  expect_equal(ncol(species), 2L,      label = "exactly two columns")
  expect_equal(colnames(species), c("Genus", "Species"),
               label = "columns named Genus and Species")
})

test_that("assignSpecies assigns at least one exact match", {
  # The bundled species reference was designed to pair with example_seqs.fa.
  # INPUT DATA: test_seqs + test_species_fa — see helper-data.R
  species <- assignSpecies(
    seqs     = test_seqs,
    refFasta = test_species_fa,
    verbose  = FALSE
  )

  expect_true(any(!is.na(species[, "Genus"])),
              label = "at least one genus assignment")
})
