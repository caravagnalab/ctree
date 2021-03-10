
# ctree <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/ctree/workflows/R-CMD-check/badge.svg)](https://github.com/caravagnalab/ctree/actions)

[![pkgdown](https://github.com/caravagnalab/ctree/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/ctree/actions/workflows/pkgdown.yaml)

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)

<!-- badges: end -->

Thee `ctree` package provides a simple implementation of *clones trees*
in cancer.

clone trees that can be build of Cancer Cell Franctions (CCFs) clusters
computed by tumour subclonal deconvolution. This type of models can be
used to study the evolutionary trajectories of a tumour from bulk
sequencing data.

Clone trees can be built from Cancer Cell Franctions (CCFs) data
reporting fraction of cancer cells harbouring a somatic mutation
(usually, a single-nucleotide variant) across a number of tumour
regions. This type of models can be used to study the evolutionary
trajectories of a tumour from bulk sequencing data, especially in setups
where one can generate whole-genome sequencing data.

-   `ctree` provides an S3-based implementation of mutation trees, as
    well as a Monte Carlo sampler that has been discussed in [Caravagna
    et al; PMID:
    30171232](https://www.ncbi.nlm.nih.gov/pubmed/30171232). The package
    provides also functions to plot and analyze the trees.

-   The sibling of a clone tree is a *mutation tree*, which is obtained
    from the analysis of binary absence/ presence profiles of somatic
    variants from bulk targeted panels; mutation trees are implemented
    in the [mtree package](https://caravagn.github.io/mtree).

`ctree` is part of the `evoverse`, a [set of R
packages](https://caravagn.github.io/evoverse) to implement Cancer
Evolution analyses.

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagn.github.io/ctree/-yellow.svg)](https://caravagn.github.io/ctree)

------------------------------------------------------------------------

### Installation

You can install the released version of `ctree` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagn/ctree")
```

------------------------------------------------------------------------

#### Copyright and contacts

Giulio Caravagna, PhD. *Institute of Cancer Research, London, UK*.

[![](https://img.shields.io/badge/Email-gcaravagn@gmail.com-seagreen.svg)](mailto:gcaravagn@gmail.com)
[![](https://img.shields.io/badge/Github-caravagn-seagreen.svg)](https://github.com/caravagn)
[![](https://img.shields.io/badge/Twitter-@gcaravagna-steelblue.svg)](https://twitter.com/gcaravagna)
[![](https://img.shields.io/badge/Personal%20webpage-https://bit.ly/2kc9E6Y-red.svg)](https://sites.google.com/site/giuliocaravagna/)
