# library(testthat); library(alabaster.string); source("test-XStringSet.R")

library(Biostrings)
stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
                                                   
test_that("saving and loading works without names", {
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with names", {
    names(stuff) <- paste0("SEQ_", seq_along(stuff))

    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with mcols", {
    mcols(stuff)$foo <- sample(letters, length(stuff))
    mcols(stuff)$bar <- runif(length(stuff))

    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with metadata", {
    metadata(stuff)$name <- "Aaron is the greatest"

    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})

test_that("saving and loading works with RNA", {
    rstuff <- RNAStringSet(c("AAA", "CC", "G", "UUUU"))
                                                   
    tmp <- tempfile()
    saveObject(rstuff, tmp)
    expect_identical(readObject(tmp), rstuff)
})

test_that("saving and loading works with proteins", {
    astuff <- AAStringSet(c("AAA", "PPPP", "FFFF", "IIII"))

    tmp <- tempfile()
    saveObject(astuff, tmp)
    expect_identical(readObject(tmp), astuff)
})

test_that("saving and loading works with custom things", {
    stuff <- BStringSet(c("XXX", "uuu", "???", "acgt"))
    tmp <- tempfile()
    saveObject(stuff, tmp)
    expect_identical(readObject(tmp), stuff)
})
