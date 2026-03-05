# Get de novo (reference) signatures

Get de novo (reference) signatures

## Usage

``` r
get_denovo_signatures(x, types = get_types(x), matrix = FALSE)
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to retrieve signames for.

- matrix:

  Logical. If \`TRUE\`, it will return the signatures in wide format.

## Value

De novo signature matrix in long or wide format names with names equal
to \`types\`.
