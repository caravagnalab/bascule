# Plot clustering centroids

Plot clustering centroids

## Usage

``` r
plot_centroids(
  x,
  types = get_types(x),
  clusters = get_cluster_labels(x),
  cls = NULL,
  sort_by = NULL,
  exposure_thr = 0,
  quantile_thr = 0,
  signatures_list = get_signames(x),
  ...
)
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to visualize.

- clusters:

  List of clusters to visualize.

- cls:

  Custom color palette for signatures.

- sort_by:

  Signature to sort patients' exposures by.

- exposure_thr:

  Only signatures with exposures greater than \`exposure_thr\` in all
  samples will be highlighted.

- quantile_thr:

  add

- signatures_list:

  add

- ...:

  Additional arguments

## Value

ggplot object.
