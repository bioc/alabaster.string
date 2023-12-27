.onLoad <- function(libname, pkgname) {
    registerReadObjectFunction("sequence_string_set", readXStringSet)
}

.onUnload <- function(libname, pkgname) {
    registerReadObjectFunction("sequence_string_set", NULL)
}
