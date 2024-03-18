# Save biological sequence data to file

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/alabaster.string.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/alabaster.string.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/alabaster.string/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/alabaster.string.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/alabaster.string.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/alabaster.string/)|

The **alabaster.string** package implements methods for saving and loading `XStringSet` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing biological sequences along with any sequence-specific metadata.
To get started, install the package and its dependencies from [Bioconductor](https://bioconductor.org/packages/alabaster.string):

```r
# install.packages("BiocManager")
BiocManager::install("alabaster.string")
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
saveObject(qdna1, tmp)

roundtrip <- readObject(tmp)
class(roundtrip)
## [1] "QualityScaledDNAStringSet"
## attr(,"package")
## [1] "Biostrings"
```
