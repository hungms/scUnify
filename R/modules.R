#' load_module
#'
#' Initializes modules
#' 
#' @param module module name
#' @export
load_module <- function(module, pkgname = "scUnify"){
    dir <- system.file(package = pkgname)
    scripts <- list.files(paste0(dir, "/modules/", module), full.names = TRUE)
    print(scripts)
    module_env <- new.env(parent = emptyenv())
    for(script in scripts){
      source(script, local = module_env)}}