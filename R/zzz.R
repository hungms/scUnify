.onLoad <- function(libname, pkgname) {
  dir <- system.file("essentials", package = pkgname)
  scripts <- list.files(dir, full.names = TRUE)
  for(script in scripts) {
    source(script)
  }
}
