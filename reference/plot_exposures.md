# Function to visualize the estimated exposures

Function to visualize the estimated exposures

## Usage

``` r
plot_exposures(
  x,
  types = get_types(x),
  samples = get_samples(x),
  clusters = get_cluster_labels(x),
  sample_name = FALSE,
  color_palette = NULL,
  add_centroid = FALSE,
  sort_by = NULL,
  exposure_thr = 0,
  quantile_thr = 0,
  signatures_list = get_signames(x)
)
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to visualize.

- samples:

  List of samples to visualize.

- clusters:

  List of clusters to visualize.

- sample_name:

  Logical. If \`TRUE\`, sample names will be reported on the x axis.

- color_palette:

  Custom color palette for signatures.

- add_centroid:

  Logical. If \`TRUE\`, also clustering's centroids will be plotted.

- sort_by:

  Signature to sort patients' exposures by.

- exposure_thr:

  Only signatures with exposures greater than \`exposure_thr\` in all
  samples will be highlighted.

- quantile_thr:

  add

- signatures_list:

  add

## Value

ggplot object.
