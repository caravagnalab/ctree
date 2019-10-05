.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c(
    'pio',
    'easypar',
    'tidyverse',
    'tidygraph',
    'ggraph',
    'crayon',
    'igraph',
    'ggrepel',
    'RColorBrewer',
    'clisymbols',
    'entropy',
    'matrixcalc',
    'reshape2'
  )
  
  suppressMessages(sapply(requirements, require, character.only = TRUE))
  
  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-
  
  ctree_welcome_message =  getOption('ctree_welcome_message', default = TRUE)
  
  if (ctree_welcome_message)
  {
    pio::pioHdr('ctree - Clone Trees in cancer')
    pio::pioStr("Author : ",
                "Giulio Caravagna <gcaravagn@gmail.com>",
                suffix = '\n')
    pio::pioStr("GitHub : ", "caravagn/ctree", suffix = '\n')
    
    options(ctree_welcome_message = FALSE)
  }
  
  invisible()
}
