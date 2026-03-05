# Get exposures

Get exposures

## Usage

``` r
get_exposure(
  x,
  types = get_types(x),
  samples = get_samples(x),
  clusters = get_cluster_labels(x),
  add_groups = FALSE,
  matrix = FALSE
)
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to retrieve signames for.

- samples:

  List of samples to report exposures for.

- clusters:

  List of cluster labels to report exposures for.

- add_groups:

  Logical. If \`TRUE\` it will add a column with the sample's group
  label.

- matrix:

  Logical. If \`TRUE\`, it will return the signatures in wide format.

## Value

Exposures matrix in long or wide format names with names equal to
\`types\`.
