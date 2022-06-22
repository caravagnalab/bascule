
# basilica <a href="https://caravagnalab.github.io/basilica/"><img src="man/figures/logo.png" align="right" height="139" /></a>

Basilica is a hierarchical Bayesian model to fit single-nucleotide
substitution signatures from multiple groups of patients, leveraging a
pre-existing catalogue of known signatures such as COSMIC. Basilica
searches for known signatures from the input catalogue as well as for
new signatures that outside the catalogue, accounting for hidden
structure in the input data (e.g., distinct tumour types). The model
uses non-negative matrix factorisation and variational inference
implemented in the
[pybasilica](https://github.com/caravagnalab/pybasilica) Python package.

#### Citation

[![](https://img.shields.io/badge/doi-10.1101/2021.02.02.429335-red.svg)](https://doi.org/....)

If you use `basilica`, please cite:

-   ……

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/basilica/-steelblue.svg)](https://caravagnalab.github.io/basilica)

### Installation

You can install the released version of `basilica` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/basilica")
```

------------------------------------------------------------------------

#### Copyright and contacts

Azad Sadr, Giulio Caravagna. Cancer Data Science (CDS) Laboratory,
University of Trieste, Italy.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
