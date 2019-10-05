.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c(
    'tidyverse',
    'tidygraph',
    'igraph',
    'ggraph',
    'pio',
    'easypar',
    'crayon',
    'clisymbols',
    'RColorBrewer',
    'entropy',
    'ggrepel'
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
