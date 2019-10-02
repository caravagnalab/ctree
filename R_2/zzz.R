# save loading of quartz for windows

.onLoad <- function(libname, pkgname) {
  
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c('tidyverse', 'tidygraph', 'ggraph', 'crayon', 'igraph', 'ggrepel')
  
  suppressMessages(sapply(requirements, require, character.only = TRUE))
  
  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-
  
  ctree_welcome_message =  getOption('ctree_welcome_message', default = TRUE)
  
  if(ctree_welcome_message)
  {
    pio::pioHdr('ctree - A package for clone trees')
    pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    pio::pioStr("GitHub : ", "caravagn/ctree", suffix = '\n')
    
    options(CNAqc_welcome_message = FALSE)
  }
  
  invisible()
}

