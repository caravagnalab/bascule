# Plot model selection scores

Function to plot the model selection scores for NMF on each variant type
and clustering. Reported scores are BIC and negative log-likelihood.

## Usage

``` r
plot_scores(x, types = get_types(x), remove_outliers = FALSE)
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to visualize.

- remove_outliers:

  Logical. If \`TRUE\`, outliers in each score will be removed.

## Value

ggplot2 object.
