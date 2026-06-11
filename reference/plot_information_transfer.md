# Plot the information transfer of a tree (drivers ordering).

The information transfer of the tree is the set of orderings associated
to the internal annotated driver events.

This function plots these orderings following a topological sort of the
node of the corresponding clone tree.

## Usage

``` r
plot_information_transfer(
  x,
  node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
  tree_layout = "tree",
  ...
)
```

## Arguments

- x:

  A `ctree` tree.

- node_palette:

  A function that can return, for an input number, a number of colours.

- tree_layout:

  Layout of this model, as of `ggraph`.

- ...:

  Other parameters, not used in this case.

## Value

A `ggplot` object for the plot.

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

plot_information_transfer(x[[1]])


# Change layout -- use igraph's "kk" layout
plot_information_transfer(x[[1]], tree_layout = 'kk')
```
