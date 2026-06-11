# Example ctree dataset.

Example ctree dataset.

## Usage

``` r
data(ctree_input)
```

## Format

A list of parameters to build a clone tree,

## Examples

``` r
data(ctree_input)
ctree_input
#> $CCF_clusters
#> # A tibble: 7 × 7
#>   cluster nMuts is.driver is.clonal    R1    R2    R3
#>   <chr>   <int> <lgl>     <lgl>     <dbl> <dbl> <dbl>
#> 1 1          72 TRUE      FALSE      0     0.92  0   
#> 2 2          69 TRUE      TRUE       0.99  0.98  0.99
#> 3 3          48 FALSE     FALSE      0     0     0.49
#> 4 4          29 FALSE     FALSE      0.01  0.01  0.93
#> 5 5          24 TRUE      FALSE      0.78  0     0   
#> 6 6          23 TRUE      FALSE      0.98  0.03  0.98
#> 7 7          15 FALSE     FALSE      0     0.41  0   
#> 
#> $samples
#> [1] "R1" "R2" "R3"
#> 
#> $drivers
#> # A tibble: 7 × 8
#>   patientID variantID is.driver is.clonal cluster    R1    R2    R3
#>   <chr>     <chr>     <lgl>     <lgl>     <chr>   <dbl> <dbl> <dbl>
#> 1 CRUK0002  RB1       TRUE      FALSE     1        0     0.92  0   
#> 2 CRUK0002  IKZF1     TRUE      FALSE     1        0     0.92  0   
#> 3 CRUK0002  KRAS      TRUE      FALSE     1        0     0.93  0   
#> 4 CRUK0002  MET       TRUE      TRUE      2        0.99  0.98  0.99
#> 5 CRUK0002  TERT      TRUE      TRUE      2        0.99  0.98  0.99
#> 6 CRUK0002  NF1       TRUE      FALSE     5        0.78  0     0   
#> 7 CRUK0002  EP300     TRUE      FALSE     6        0.96  0.03  0.98
#> 
#> $sspace.cutoff
#> [1] 10000
#> 
#> $n.sampling
#> [1] 5000
#> 
#> $store.max
#> [1] 100
#> 
#> $patient
#> [1] "CUK12345"
#> 
```
