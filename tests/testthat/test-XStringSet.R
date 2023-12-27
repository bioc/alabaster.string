# This checks whether we can save and load DNAStringSets correctly.
# library(testthat); library(alabaster.string); source("test-XStringSet.R")

library(Biostrings)
stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
                                                   
test_that("saving and loading works without names", {
    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(stuff, tmp, path="sequence_thing")
    .writeMetadata(info, tmp)

    seq.path <- info$sequence_string_set$sequence_file$resource$path
    expect_match(seq.path, ".fa.gz$")

    seq.meta <- acquireMetadata(tmp, seq.path)
    expect_identical(seq.meta$fasta_file$type, "DNA")
    expect_identical(seq.meta$fasta_file$compression, "gzip")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, stuff)

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with names", {
    tmp <- tempfile()
    dir.create(tmp)
    names(stuff) <- paste0("SEQ_", seq_along(stuff))

    info <- stageObject(stuff, tmp, path="sequence_thing")
    .writeMetadata(info, tmp)
    expect_match(info$sequence_string_set$sequence_file$resource$path, ".fa.gz$")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, stuff)

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with mcols", {
    tmp <- tempfile()
    dir.create(tmp)
    mcols(stuff)$foo <- sample(letters, length(stuff))
    mcols(stuff)$bar <- runif(length(stuff))

    info <- stageObject(stuff, tmp, path="sequence_thing")
    .writeMetadata(info, tmp)
    expect_match(info$sequence_string_set$sequence_data$resource$path, ".csv.gz$")

    roundtrip <- loadXStringSet(info, tmp)
    expect_equal(roundtrip, stuff)

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with metadata", {
    tmp <- tempfile()
    dir.create(tmp)
    metadata(stuff)$name <- "Aaron is the greatest"

    info <- stageObject(stuff, tmp, path="sequence_thing")
    .writeMetadata(info, tmp)
    expect_match(info$sequence_string_set$other_data$resource$path, ".$")

    roundtrip <- loadXStringSet(info, tmp)
    expect_equal(roundtrip, stuff)

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with RNA", {
    rstuff <- RNAStringSet(c("AAA", "CC", "G", "UUUU"))
                                                   
    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(rstuff, tmp, path="sequence_thing")
    .writeMetadata(info, tmp)
    seq.path <- info$sequence_string_set$sequence_file$resource$path
    expect_match(seq.path, ".fa.gz$")

    seq.meta <- acquireMetadata(tmp, seq.path)
    expect_identical(seq.meta$fasta_file$type, "RNA")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, rstuff)

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with proteins", {
    astuff <- AAStringSet(c("AAA", "PPPP", "FFFF", "IIII"))
                                                   
    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(astuff, tmp, path="sequence_thing")
    .writeMetadata(info, tmp)
    seq.path <- info$sequence_string_set$sequence_file$resource$path
    expect_match(seq.path, ".fa.gz$")

    seq.meta <- acquireMetadata(tmp, seq.path)
    expect_identical(seq.meta$fasta_file$type, "AA")
    expect_identical(seq.meta$fasta_file$compression, "gzip")

    roundtrip <- loadXStringSet(info, tmp)
    expect_identical(roundtrip, astuff)

    # Trying in the new world.
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with custom things", {
    stuff <- BStringSet(c("XXX", "uuu", "???", "acgt"))
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})
