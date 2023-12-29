#' Read an XStringSet from disk.
#'
#' Read a \linkS4class{XStringSet} object from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{saveObject}} method for \linkS4class{XStringSet} objects.
#' @param metadata Named list of metadata for this object, see \code{\link{readObjectFile}} for details.
#' @param ... Further arguments passed to internal \code{\link{altReadObject}} calls.
#'
#' @return An \linkS4class{XStringSet} subclass containing DNA, RNA, protein or custom sequences.
#' This may also be a \linkS4class{QualityScaledDNAStringSet} with quality scores.
#'
#' @seealso
#' \code{"\link{saveObject,XStringSet-method}"}, to save an \linkS4class{XStringSet} to disk.
#'
#' @examples
#' library(Biostrings)
#' stuff <- DNAStringSet(c("AAA", "CC", "G", "TTTT"))
#'
#' tmp <- tempfile()
#' saveObject(stuff, tmp)
#' readObject(tmp)
#' 
#' @export
#' @aliases loadXStringSet
#' @importFrom Biostrings readDNAStringSet readRNAStringSet readAAStringSet readBStringSet readQualityScaledDNAStringSet
#' @importFrom alabaster.base acquireFile acquireMetadata .restoreMetadata
#' @importFrom utils read.table
readXStringSet <- function(path, metadata, ...) {
    info <- metadata$sequence_string_set

    if (!is.null(info$quality_type) && info$quality_type != "none") {
        if (info$sequence_type != "DNA") {
            stop("unknown quality scaled sequence type '", info$sequence_type, "'")
        }
        if (info$quality_type == "phred") {
            if (info$quality_offset == 33) {
                qscore <- "phred"
            } else {
                qscore <- "illumina"
            }
        } else {
            qscore <- "solexa"
        }
        x <- readQualityScaledDNAStringSet(file.path(path, "sequences.fastq.gz"), quality.scoring=qscore)

    } else {
        seq.file <- file.path(path, "sequences.fasta.gz")
        if (info$sequence_type == "DNA") {
            x <- readDNAStringSet(seq.file)
        } else if (info$sequence_type == "RNA") {
            x <- readRNAStringSet(seq.file)
        } else if (info$sequence_type == "AA") {
            x <- readAAStringSet(seq.file)
        } else {
            x <- readBStringSet(seq.file)
        }
    }

    name.path <- file.path(path, "names.txt.gz")
    if (file.exists(name.path)) {
        names(x) <- read.table(name.path)[,1]
    } else {
        names(x) <- NULL
    }

    readMetadata(x, mcols.path=file.path(path, "sequence_annotations"), metadata.path=file.path(path, "other_annotations"), ...)
}

##################################
######### OLD STUFF HERE #########
##################################

#' @export
loadXStringSet <- function(seq.info, project) {
    seq.meta <- acquireMetadata(project, path=seq.info$sequence_string_set$sequence_file$resource$path)
    seq.file <- acquireFile(project, path=seq.info$sequence_string_set$sequence_file$resource$path)

    if ("fastq_file" %in% names(seq.meta)) {
        type <- seq.meta$fastq_file$type
        if (is.null(type) || type == "DNA") {
            x <- readQualityScaledDNAStringSet(seq.file, quality.scoring=tolower(seq.meta$fastq_file$quality_encoding))
        } else {
            stop("unknown quality scaled sequence type '", type, "'")
        }
    } else {
        type <- seq.meta$fasta_file$type
        if (is.null(type) || type == "DNA") {
            x <- readDNAStringSet(seq.file)
        } else if (type == "RNA") {
            x <- readRNAStringSet(seq.file)
        } else if (type == "AA") {
            x <- readAAStringSet(seq.file)
        } else {
            stop("unknown sequence type '", type, "'")
        }
    }

    if (!isTRUE(seq.info$sequence_string_set$names)) {
        names(x) <- NULL
    }

    .restoreMetadata(x, mcol.data=seq.info$sequence_string_set$sequence_data, meta.data=seq.info$sequence_string_set$other_data, project=project)
}
