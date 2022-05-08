
# ctree <a href="https://caravagnalab.github.io/ctree/"><img src="man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/ctree/workflows/R-CMD-check/badge.svg)](https://github.com/caravagnalab/ctree/actions)
[![pkgdown](https://github.com/caravagnalab/ctree/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/ctree/actions/workflows/pkgdown.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

The `ctree` package implements clones trees for cancer evolutionary
studies. These models are built from Cancer Cell Franctions (CCFs)
clusters computed via tumour subclonal deconvolution, using either one
or more tumour biopsies at once. They can be used to model evolutionary
trajectories from bulk sequencing data, especially if whole-genome
sequencing is available. The package implements S3 objects for the
mutation trees, as well as a Monte Carlo sampler to generate them, as
well as functions to plot and analyze the trees. The sibling of a clone
tree is a mutation tree, which is built from binary mutation profiles;
refer to the [mtree package](https://caravagn.github.io/mtree) for
mutation trees.

#### Citation

[![](https://img.shields.io/badge/doi-10.1038/s41592--018--0108--x-red.svg)](https://doi.org/10.1038/s41592-018-0108-x)

Please cite this if you use `ctree`:

-   G. Caravagna, Y. Giarratano, D. Ramazzoti, I. Tomlinson, T.A.
    Graham, G. Sanguinetti, A. Sottoriva. *Detecting repeated cancer
    evolution from multi-region tumor sequencing data.* Nature Methods
    15, 707–714 (2018).

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/ctree/-yellow.svg)](https://caravagnalab.github.io/ctree)

------------------------------------------------------------------------

### Installation

You can install the released version of `ctree` with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/ctree")
```

------------------------------------------------------------------------

#### Copyright and contacts

Giulio Caravagna. Cancer Data Science (CDS) Laboratory.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
