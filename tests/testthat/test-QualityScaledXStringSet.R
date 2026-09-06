# library(testthat); library(alabaster.string); source("test-QualityScaledXStringSet.R")

library(Biostrings)
stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
scores <- NumericList(lapply(width(stuff), FUN=runif, min=0, max=0.01))

test_that("saving and loading works for quality scaled DNAStringSets", {
    qstuff <- QualityScaledDNAStringSet(stuff, PhredQuality(scores))

    tmp <- tempfile()
    saveObject(qstuff, tmp)
    expect_identical(readObject(tmp), qstuff)
})

test_that("saving and loading acknowledges the quality encoding type", {
    # Works with the Illumina encoding:
    qstuff <- QualityScaledDNAStringSet(stuff, IlluminaQuality(scores))

    tmp <- tempfile()
    saveObject(qstuff, tmp)
    expect_identical(readObject(tmp), qstuff)

    # Works with the Solexa encoding:
    qstuff <- QualityScaledDNAStringSet(stuff, SolexaQuality(scores))

    tmp <- tempfile()
    saveObject(qstuff, tmp)
    expect_identical(readObject(tmp), qstuff)
})
