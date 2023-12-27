#' Save a XStringSet
#'
#' Save a \linkS4class{XStringSet} to its on-disk representation.
#'
#' @param x A \linkS4class{XStringSet} or any of its subclasses, in particular a \linkS4class{QualityScaledXStringSet}.
#' @inheritParams alabaster.base::saveObject
#'
#' @return The contents of \code{x} are saved into a \code{path}, and \code{NULL} is invisibly returned.
#'
#' @author Aaron Lun
#' @examples
#' library(Biostrings)
#' stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
#'
#' tmp <- tempfile()
#' saveObject(stuff, tmp)
#' list.files(tmp, recursive=TRUE)
#'
#' @export
#' @rdname saveXStringSet
#' @aliases stageObject,XStringSet-method
#' @importFrom Biostrings writeXStringSet writeQualityScaledXStringSet quality
#' @importClassesFrom Biostrings XStringSet QualityScaledXStringSet 
#' @importFrom S4Vectors mcols
#' @import alabaster.base methods
#' @importFrom utils write.table
setMethod("saveObject", "XStringSet", function(x, path, ...) {
    dir.create(path, showWarnings=FALSE)

    details <- list(version="1.0", length=length(x))

    if (is(x, "DNAStringSet")) {
        details$sequence_type <- "DNA"
    } else if (is(x, "RNAStringSet")) {
        details$sequence_type <- "RNA"
    } else if (is(x, "AAStringSet")) {
        details$sequence_type <- "AA"
    } else {
        details$sequence_type <- "custom"
    }

    has.qual <- FALSE
    if (is(x, "QualityScaledXStringSet")) {
        has.qual <- TRUE
        Q <- quality(x)
        if (is(Q, "PhredQuality")) {
            details$quality_type <- "phred"
            details$quality_offset <- 33
        } else if (is(Q, "SolexaQuality")) {
            details$quality_type <- "solexa"
        } else if (is(Q, "IlluminaQuality")) {
            details$quality_type <- "phred"
            details$quality_offset <- 64
        } else {
            stop("unrecognized quality string encoding")
        }
    }

    saveObjectFile(
        path, 
        "sequence_string_set", 
        list(sequence_string_set=details)
    )

    copy <- x
    names(copy) <- seq_along(copy) - 1L
    if (has.qual) {
        writeQualityScaledXStringSet(copy, filepath=file.path(path, "sequences.fastq.gz"), compress=TRUE)
    } else {
        writeXStringSet(copy, filepath=file.path(path, "sequences.fasta.gz"), compress=TRUE)
    }

    if (!is.null(names(x))) {
        con <- gzfile(file.path(path, "names.txt.gz"), open="wb")
        on.exit(close(con))
        write.table(file=con, x=data.frame(names(x)), col.names=FALSE, row.names=FALSE, quote=TRUE)
    }

    saveMetadata(x, mcols.path=file.path(path, "sequence_annotations"), metadata.path=file.path(path, "other_annotations"), ...)
    invisible(NULL)
})

##################################
######### OLD STUFF HERE #########
##################################

#' @export
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

    meta <- list(is_child=TRUE)
    inner <- list(compression="gzip", type=type)

    if (is(x, "QualityScaledXStringSet")) {
        target <- "sequence.fastq.gz"
        meta$path <- file.path(path, target)
        meta[["$schema"]] <- "fastq_file/v1.json"

        Q <- quality(copy)
        encoding <- switch(as.character(class(Q)),
            PhredQuality="phred",
            SolexaQuality="solexa",
            IlluminaQuality="illumina",
            stop("unrecognized quality string encoding")
        )
        inner$quality_encoding <- encoding
        meta$fastq_file <- inner

        writeQualityScaledXStringSet(copy, filepath=file.path(full, target), compress=TRUE)
    } else {
        target <- "sequence.fa.gz"
        meta$path <- file.path(path, target)
        meta[["$schema"]] <- "fasta_file/v1.json"
        meta$fasta_file <- inner
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
