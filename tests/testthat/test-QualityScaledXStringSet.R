# This checks whether we can save and load DNAStringSets correctly.
# library(testthat); library(alabaster.string); source("test-QualityScaledXStringSet.R")

library(Biostrings)
stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
scores <- NumericList(lapply(width(stuff), FUN=runif, min=0, max=0.01))

test_that("saving and loading works for quality scaled DNAStringSets", {
    tmp <- tempfile()
    dir.create(tmp)

    qstuff <- QualityScaledDNAStringSet(stuff, PhredQuality(scores))
    info <- stageObject(qstuff, tmp, path="dna_thing")
    .writeMetadata(info, tmp)
    expect_match(info$sequence_string_set$sequence_file$resource$path, ".fastq.gz$")

    seq.meta <- alabaster.base::acquireMetadata(tmp, info$sequence_string_set$sequence_file$resource$path)
    expect_identical(seq.meta$sequence_file$quality_encoding, "phred")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, qstuff)
})

test_that("saving and loading acknowledges the quality encoding type", {
    tmp <- tempfile()
    dir.create(tmp)

    # Works with the Illumina encoding:
    qstuff <- QualityScaledDNAStringSet(stuff, IlluminaQuality(scores))
    info <- stageObject(qstuff, tmp, path="dna_thing")
    .writeMetadata(info, tmp)
    seq.meta <- alabaster.base::acquireMetadata(tmp, info$sequence_string_set$sequence_file$resource$path)
    expect_identical(seq.meta$sequence_file$quality_encoding, "illumina")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, qstuff)

    # Works with the Solexa encoding:
    qstuff <- QualityScaledDNAStringSet(stuff, SolexaQuality(scores))
    info <- stageObject(qstuff, tmp, path="dna_thing_solexa")
    .writeMetadata(info, tmp)
    seq.meta <- alabaster.base::acquireMetadata(tmp, info$sequence_string_set$sequence_file$resource$path)
    expect_identical(seq.meta$sequence_file$quality_encoding, "solexa")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, qstuff)
})
