# Function to visualize the estimated and reference signatures.

Function to visualize the estimated and reference signatures.

## Usage

``` r
plot_signatures(
  x,
  types = get_types(x),
  context = T,
  cls = NULL,
  signames = get_signames(x)
)
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to visualize.

- context:

  Logical. If \`TRUE\`, context names will be reported on the x axis.

- cls:

  Custom color palette for signatures.

- signames:

  List of signatures to visualize.

## Value

ggplot object.
