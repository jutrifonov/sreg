startupMessage <- function() {
  version <- packageVersion("sreg")
  msg <- c(paste0(
    "  ____  ____  _____ ____      Stratified Randomized\n",
    " / ___||  _ \\| ____/ ___|     Experiments\n",
    " \\___ \\| |_) |  _|| |  _  \n",
    "  ___) |  _ <| |__| |_| |  \n",
    " |____/|_| \\_\\_____\\____| version ", version, "\n",
    "                           \n"
  ))
  return(msg)
}




.onAttach <- function(lib, pkg) {
  msg <- startupMessage()
  if (!interactive()) {
    msg[1] <- paste("Package 'sreg' version", packageVersion("sreg"))
  }
  packageStartupMessage(msg)
  invisible()
}
