# Function to map de novo signatures to input catalogues

Function to map de novo signatures to input catalogues

## Usage

``` r
convert_dn_names(x, x.simul = NULL, reference_cat = NULL, cutoff = 0.8)
```

## Arguments

- x:

  Object of class "bascule_obj"

- x.simul:

  Another object of class "bascule_obj". If present, \`x\` de novo
  signatures will be mapped to \`x_simul\` ones

- reference_cat:

  List of reference catalogues. If not \`NULL\`, \`x\` de novo
  signatures will be mapped to these catalogues.

- cutoff:

  Threshold for the cosine similarity.

## Value

Modified version of \`x\` with the mapped de novo names renamed.
