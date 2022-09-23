# Save biological sequence data to file

The **alabaster.string** package implements methods for saving and loading `XStringSet` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing biological sequences along with any sequence-specific metadata.
To get started, install the package and its dependencies from GitHub:

```r
devtools::install_github("ArtifactDB/alabaster.schemas")
devtools::install_github("ArtifactDB/alabaster.base")
devtools::install_github("ArtifactDB/alabaster.string")
```

In the example below, we save a `QualityScaledDNAStringSet` object to file:

```r
library(Biostrings)
example(QualityScaledDNAStringSet, echo=FALSE) # can't be bothered to copy it here.
qdna1
##   A QualityScaledDNAStringSet instance containing:
## 
## DNAStringSet object of length 2:
##     width seq
## [1]     4 TTGA
## [2]     4 CTCN
## 
## PhredQuality object of length 2:
##     width seq
## [1]     4 *+,-
## [2]     4 6789

library(alabaster.string)
tmp <- tempfile()
dir.create(tmp)
meta <- stageObject(qdna1, tmp, "dna")
meta[["$schema"]]
## [1] "sequence_string_set/v1.json"

roundtrip <- loadObject(meta, tmp)
class(roundtrip)
## [1] "QualityScaledDNAStringSet"
## attr(,"package")
## [1] "Biostrings"
```
