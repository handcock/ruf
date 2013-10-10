.onAttach <- function(lib, pkg){
  info <- packageDescription("ruf")

  packageStartupMessage(
    paste('\nruf: version ', info$Version, ', created on ', info$Date, '\n',
          "Copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",
          'For citation type citation("ruf").\n', sep="")
                        )
  
}
