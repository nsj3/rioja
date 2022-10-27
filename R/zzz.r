.onAttach <- function(lib, pkg)  {
   line1 <- paste("This is rioja ", utils::packageDescription("rioja", fields="Version"))
   packageStartupMessage(paste(line1), appendLF = TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("rioja", libpath)
}
