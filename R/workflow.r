### For work flows only.

get.workflow <- function(dir.name = "workflow/",
    parallel = c("spmd"), pkg = "MixfMRI"){
  dir.name <- paste(dir.name, parallel[1], "/", sep = "")
  file.path <- tools::file_path_as_absolute(
                 system.file(dir.name, package = "MixfMRI"))
  file.path
} # End of get.workflow().

cp.workflow <- function(pkg = "MixfMRI", to = "."){
  path.current <- paste(to, "/tmp_workflow", sep = "")
  if(!file.exists(path.current)){
    dir.create(path.current, mode = "0755")
  }

  path.workflow <- get.workflow(pkg = pkg)

  case <- c("2d_I", "2d_V", "2d_noX", "3d_I", "3d_V", "3d_noX",
            "simu_2d")
  for(i.case in case){
    path.dir <- paste(path.workflow, "/", i.case, sep = "")
    file.copy(path.dir, path.current, overwrite = TRUE, recursive = TRUE)
  }
} # End of cp.workflow().
