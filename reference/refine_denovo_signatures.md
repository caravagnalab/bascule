# De novo signatures refinement.

Function to refine the inferred de novo signatures. The function
performs a linear combination between de novo and reference signatures.
If a de novo signature can be explained as a linear combination of one
or more reference signatures, it will be removed and its exposures will
be distributed among the similar signatures.

## Usage

``` r
refine_denovo_signatures(x, types = get_types(x))
```

## Arguments

- x:

  bascule object.

- types:

  List of variant types to perform de novo refinement on.

## Value

bascule object.
