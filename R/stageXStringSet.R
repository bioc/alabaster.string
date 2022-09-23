#' Stage a XStringSet
#'
#' Stage a XStringSet by saving it to the appropriate file format.
#'
#' @param x A \linkS4class{XStringSet} or any of its subclasses, in particular a \linkS4class{QualityScaledXStringSet}.
#' @inheritParams alabaster.base::stageObject
#'
#' @return A list containing metadata for \code{x}.
#' A subdirectory is created at \code{path} inside \code{dir} and the contents of \code{x} are saved to various files within that subdirectory.
#' If \code{x} is a QualityScaledXStringSet, a FASTQ file is created instead of a FASTA file.
#'
#' @author Aaron Lun
#' @examples
#' library(Biostrings)
#' stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
#'
#' tmp <- tempfile()
#' dir.create(tmp)
#' stageObject(stuff, tmp, path="dna_thing")
#' list.files(tmp, recursive=TRUE)
#'
#' @export
#' @rdname stageXStringSet
#' @importFrom Biostrings writeXStringSet writeQualityScaledXStringSet quality
#' @importClassesFrom Biostrings XStringSet QualityScaledXStringSet 
#' @import alabaster.base methods
setMethod("stageObject", "XStringSet", function(x, dir, path, child=FALSE, ...) {
    full <- file.path(dir, path)
    dir.create(full, showWarnings=FALSE)

    # Writing our thing to FASTA or FASTQ.
    copy <- x
    if (is.null(names(copy))) {
        names(copy) <- seq_along(copy)
    }
   
    type <- NULL
    if (is(x, "DNAStringSet")) {
        type <- "DNA"
    } else if (is(x, "RNAStringSet")) {
        type <- "RNA"
    } else if (is(x, "AAStringSet")) {
        type <- "AA"
    } else {
        stop("unknown sequence type for '", class(x)[1], "'")
    }

    meta <- list(
        `$schema`="sequence_file/v1.json",
        is_child=TRUE,
        sequence_file=list(
            type=type,
            compression="gzip"
        )
    )

    if (is(x, "QualityScaledXStringSet")) {
        target <- "sequence.fastq.gz"
        meta$path <- file.path(path, target)
        meta$sequence_file$format <- "FASTQ"

        Q <- quality(copy)
        meta$sequence_file$quality_encoding <- switch(as.character(class(Q)),
            PhredQuality="phred",
            SolexaQuality="solexa",
            IlluminaQuality="illumina",
            stop("unrecognized quality string encoding")
        )

        writeQualityScaledXStringSet(copy, filepath=file.path(full, target), compress=TRUE)
    } else {
        target <- "sequence.fa.gz"
        meta$path <- file.path(path, target)
        meta$sequence_file$format <- "FASTA"
        writeXStringSet(copy, filepath=file.path(full, target), compress=TRUE)
    }

    out <- .writeMetadata(meta, dir=dir)

    # Processing mcols and metadata.
    mcol.data <- .processMcols(x, dir=dir, path=path, mcols.name="sequence_data")
    meta.data <- .processMetadata(x, dir=dir, path=path, meta.name="other_data")

    list(
        `$schema`="sequence_string_set/v1.json",
        path=file.path(path, 'set.json'),
        is_child=child,
        sequence_string_set=list(
            sequence_file=list(resource=out),
            sequence_data=mcol.data,
            other_data=meta.data,
            names=!is.null(names(x))
        )
    )
})
