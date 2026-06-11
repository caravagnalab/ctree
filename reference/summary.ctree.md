# Summary of a `"ctree"` object.

Reports some summary statistics for a `"rev_phylo"` object, which is a
bit more than just using `print.`

## Usage

``` r
# S3 method for class 'ctree'
summary(object, ...)
```

## Arguments

- object:

  A `ctree` tree.

- ...:

  Extra parameters

## Value

None.

## Examples

``` r
data(ctree_input)

x = ctrees(
ctree_input$CCF_clusters, 
ctree_input$drivers,
ctree_input$samples,
ctree_input$patient,
ctree_input$sspace.cutoff,
ctree_input$n.sampling,
ctree_input$store.max
)
#>  [ ctree ~ clone trees generator for CUK12345 ] 
#> 
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
#> ✔ Trees per region 1, 3, 1
#> ℹ Total 3 tree structures - search is exahustive
#> 
#> ── Ranking trees 
#> ✔ 3  trees with non-zero score, storing 3

summary(x[[1]])
#>  [ ctree - ctree rank 1/3 for CUK12345 ] 
#> 
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
#> Tree shape (drivers annotated)  
#> 
#>   \-GL
#>    \-2 :: MET, TERT
#>     |-1 :: RB1, IKZF1, KRAS
#>     | \-7
#>     \-6 :: EP300
#>      |-4
#>      | \-3
#>      \-5 :: NF1
#> 
#> Information transfer  
#> 
#>    MET ---> RB1 
#>    MET ---> IKZF1 
#>    MET ---> KRAS 
#>    TERT ---> RB1 
#>    TERT ---> IKZF1 
#>    TERT ---> KRAS 
#>    GL ---> MET 
#>    GL ---> TERT 
#>    EP300 ---> NF1 
#>    MET ---> EP300 
#>    TERT ---> EP300 
#> 
#> Tree score 0.6 
#> 
#> CCF clusters:  
#> 
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
#> Drivers:  
#> 
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
#> Pigeonhole principle: 12 0 
#> 
#>     R1   R2   R3
#> 1 TRUE TRUE TRUE
#> 2 TRUE TRUE TRUE
#> 4 TRUE TRUE TRUE
#> 6 TRUE TRUE TRUE
#> 
#> Goodness-of-fit: 1  
#> 
#> 
```
