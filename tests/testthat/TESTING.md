# Running the dada2 test suite

## Quick start

From an R session with the package installed or loaded via `devtools`:

```r
devtools::test()
```

Or from the shell:

```sh
Rscript -e "devtools::test('.')"
```

To run a single test file:

```r
devtools::test(filter = "filter")   # runs test-filter.R
devtools::test(filter = "pipeline") # runs test-pipeline.R
```

To run via `R CMD check` (as CI does):

```sh
R CMD check --no-manual --no-vignettes .
```

---

## Test data

By default all tests use the small example files bundled in `inst/extdata`.
To run the suite against your own FASTQ files or reference databases, set
any of the following environment variables before launching R:

| Variable | Description | Bundled fallback |
|---|---|---|
| `DADA2_TEST_FWD_1` | Sample 1 forward FASTQ (gzip OK) | `sam1F.fastq.gz` |
| `DADA2_TEST_REV_1` | Sample 1 reverse FASTQ | `sam1R.fastq.gz` |
| `DADA2_TEST_FWD_2` | Sample 2 forward FASTQ | `sam2F.fastq.gz` |
| `DADA2_TEST_REV_2` | Sample 2 reverse FASTQ | `sam2R.fastq.gz` |
| `DADA2_TEST_TRAIN_FA` | `assignTaxonomy` reference FASTA | `example_train_set.fa.gz` |
| `DADA2_TEST_SPECIES_FA` | `assignSpecies` reference FASTA | `example_species_assignment.fa.gz` |

Example (bash):

```sh
export DADA2_TEST_FWD_1=/data/sample1_R1.fastq.gz
export DADA2_TEST_REV_1=/data/sample1_R2.fastq.gz
export DADA2_TEST_FWD_2=/data/sample2_R1.fastq.gz
export DADA2_TEST_REV_2=/data/sample2_R2.fastq.gz
export DADA2_TEST_TRAIN_FA=/refs/silva_nr_v138.fa.gz
export DADA2_TEST_SPECIES_FA=/refs/silva_species_assignment_v138.fa.gz
Rscript -e "devtools::test('.')"
```

All path resolution happens in `helper-data.R`, which is the single file to
edit if you prefer hard-coded paths over environment variables.

---

## Skipping the slow integration test

`test-pipeline.R` runs the full pipeline end-to-end including `learnErrors()`
and is the slowest test. Skip it by setting:

```sh
export DADA2_SKIP_SLOW=1
```

The remaining unit tests use the bundled `tperr1` error matrix instead of
calling `learnErrors()` and complete in seconds.

---

## Test file overview

| File | Pipeline step(s) covered |
|---|---|
| `helper-data.R` | Shared data paths (edit here to use external data) |
| `test-filter.R` | `filterAndTrim()` |
| `test-derep.R` | `derepFastq()` |
| `test-denoise.R` | `learnErrors()`, `dada()` |
| `test-paired.R` | `mergePairs()` |
| `test-seqtable.R` | `makeSequenceTable()`, `removeBimeraDenovo()` |
| `test-taxonomy.R` | `assignTaxonomy()`, `assignSpecies()` |
| `test-pipeline.R` | Full end-to-end integration (slow) |
