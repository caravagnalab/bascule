# Function to visualize the similarity of denovo and reference signatures.

Function to visualize the similarity of denovo and reference signatures.

## Usage

``` r
plot_similarity_reference(
  x,
  reference = NULL,
  type = "SBS",
  similarity_cutoff = 0.8,
  context = T,
  add_pheatmap = T
)
```

## Arguments

- x:

  Bascule object.

- reference:

  External reference catalogue to compare the denovo with.

- similarity_cutoff:

  add

- context:

  Logical. If set to `TRUE`, the context labels are reported on the x
  axis.

- add_pheatmap:

  Logical. If set to `TRUE`, the heatmap with the similarity values
  among signatures will be reported.

- by_subs:

  Logical. If set to `TRUE`, the similarity is computed separately for
  each substitution.

## Value

add
