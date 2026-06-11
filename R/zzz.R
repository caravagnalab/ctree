#' @importFrom tibble as_tibble
#' @importFrom utils globalVariables
NULL

# Suppress R CMD check NOTEs for tidyverse-style non-standard
# evaluation column names used throughout dplyr/ggplot2 pipelines.
utils::globalVariables(c(
  "A", "CCF", "attachment", "cluster", "driver", "edges", "from", "to",
  "is.clonal", "is.driver", "label", "nDrivers", "nMuts", "name",
  "nodes", "pvalue", "region", "significant", "tot", "value",
  "variable", "variantID", "y"
))

.onLoad <- function(libname, pkgname)
{
  options(pio.string_fg_colour = crayon::bgYellow$black)

  invisible()
}

.onAttach <- function(libname, pkgname)
{
  packageStartupMessage(
    "Loading ctree, 'Clone trees in cancer'. Support : https://caravagnalab.github.io/ctree/"
  )
}
