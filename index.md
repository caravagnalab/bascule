# BASCULE

BASCULE is a Bayesian model to fit multiple signature types from
multiple patients, leveraging a pre-existing catalogue of known
signatures such as COSMIC. BASCULE searches for known signatures from
the input catalogue as well as for new signatures that outside the
catalogue, accounting for hidden structure in the input data (e.g.,
distinct tumour types). Moreover, bascule performs tensor clustering to
retrieve latent groups in the input cohort from the exposures of
multiple signature types jointly. The model uses non-negative matrix
factorisation and variational inference implemented in the
[pybascule](https://github.com/caravagnalab/pybascule) Python package.

#### Citation

[![doi.org/10.1186/s13059-025-03835-9](http://img.shields.io/badge/doi-10.1101/2024.09.16.613266-red.svg)](https:%3A/caravagnalab.github.io/bascule/doi.org/10.1186/s13059-025-03835-9)

If you use `bascule`, please cite:

- Buscaroli, E., Sadr, A., Bergamin, R. et al. BASCULE: bayesian
  inference and clustering of mutational signatures leveraging
  biological priors. [*Genome Biol*](NA) (2025).

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/bascule/-steelblue.svg)](https://caravagnalab.github.io/bascule)

### Installation

You can install the released version of `bascule` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/bascule")
```

If you want to fit your data with `bascule`, you also need to install
the required Python dependencies. Please follow the instructions in the
[Get
Started](https://caravagnalab.github.io/bascule/articles/bascule.html)
article.

------------------------------------------------------------------------

#### Copyright and contacts

Elena Buscaroli, Azad Sadr, Giulio Caravagna. Cancer Data Science (CDS)
Laboratory, University of Trieste, Italy.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
