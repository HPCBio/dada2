# =============================================================================
# Test data configuration — edit this file to use external datasets
# =============================================================================
#
# By default all tests use the small example files bundled in inst/extdata.
# To run the suite against your own data, set these environment variables
# before launching R (or before calling devtools::test() / R CMD check):
#
#   DADA2_TEST_FWD_1   path to first sample forward FASTQ (gzipped OK)
#   DADA2_TEST_REV_1   path to first sample reverse FASTQ (paired-end only)
#   DADA2_TEST_FWD_2   path to second sample forward FASTQ
#   DADA2_TEST_REV_2   path to second sample reverse FASTQ
#   DADA2_TEST_TRAIN_FA    path to assignTaxonomy reference FASTA (gzipped OK)
#   DADA2_TEST_SPECIES_FA  path to assignSpecies reference FASTA (gzipped OK)
#
# Any variable that is unset falls back to the corresponding bundled file.
# =============================================================================

.extdata <- system.file("extdata", package = "dada2")

# ---- Paired-end Illumina FASTQ inputs ---------------------------------------
# Bundled files: ~16S V4 amplicon reads, two mock community samples.
# Substitution: provide FASTQ files from your own sequencing run.

.get_path <- function(env_var, bundled_file) {
  val <- Sys.getenv(env_var, unset = "")
  if (nzchar(val)) val else file.path(.extdata, bundled_file)
}

test_fnF <- c(
  .get_path("DADA2_TEST_FWD_1", "sam1F.fastq.gz"),  # Sample 1 forward reads
  .get_path("DADA2_TEST_FWD_2", "sam2F.fastq.gz")   # Sample 2 forward reads
)

test_fnR <- c(
  .get_path("DADA2_TEST_REV_1", "sam1R.fastq.gz"),  # Sample 1 reverse reads
  .get_path("DADA2_TEST_REV_2", "sam2R.fastq.gz")   # Sample 2 reverse reads
)

# ---- Taxonomy reference databases -------------------------------------------
# Bundled files: small example references (genus-level NBC + exact-match species).
# Substitution: use full SILVA/RDP/GTDB releases for realistic taxonomy tests.

test_train_fa   <- .get_path("DADA2_TEST_TRAIN_FA",   "example_train_set.fa.gz")
test_species_fa <- .get_path("DADA2_TEST_SPECIES_FA", "example_species_assignment.fa.gz")

# ---- Pre-built error rate matrix --------------------------------------------
# tperr1 is bundled with dada2: a 16x41 error matrix (transitions × Q-scores).
# Using it in unit tests avoids the cost of learnErrors() while still exercising
# the dada() inference engine.
# Substitution: replace with errF / errR produced by learnErrors() on your data.
data("tperr1", envir = environment())
test_err <- tperr1
