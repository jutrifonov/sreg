startupMessage <- function() {
  version <- packageVersion("sreg")
  msg <- c(paste0(
   col_blue("  ____  ____  _____ ____      Stratified Randomized\n"),
   col_blue(" / ___||  _ \\| ____/ ___|     Experiments\n"),
   col_blue(" \\___ \\| |_) |  _|| |  _  \n"),
   col_blue("  ___) |  _ <| |__| |_| |  \n"),
   col_blue(" |____/|_| \\_\\_____\\____| "), col_cyan("version "), col_cyan(version), "\n",
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



