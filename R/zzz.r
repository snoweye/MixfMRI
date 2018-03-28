### Lastest load into a package.

.onLoad <- function(libname, pkgname){
  library.dynam("MixfMRI", pkgname, libname)
  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  library.dynam.unload("MixfMRI", libpath)
  invisible()
} # End of .onUnload().

