# bclogit <img src="bclogit/man/figures/logo.png" align="right" height="139" alt="bclogit hex logo" />

**Conditional Logistic Regression with Optional Concordant Pairs**

[![pkgdown](https://github.com/tennenbaumj/bclogit_package_and_paper_repo/actions/workflows/pkgdown.yaml/badge.svg)](https://tennenbaumj.github.io/bclogit_package_and_paper_repo/)

Full documentation and function reference: **https://tennenbaumj.github.io/bclogit_package_and_paper_repo/**

## Overview

`bclogit` is an R package for fitting conditional logistic regression models that can optionally incorporate information from concordant pairs (and a reservoir of additional controls) to improve estimation of discordant pairs. This approach effectively uses a Bayesian prior derived from the concordant data to inform the discordant pairs analysis, suitable for matched case-control studies where data sparsity might be an issue.

## Installation

You can install the development version of `bclogit` from source:

```r
# install.packages("devtools")
devtools::install()
```

or use the command line

```
R CMD INSTALL bclogit
```

## Usage

Here is a basic example of how to use `bclogit`:

```r
library(bclogit)

# Simulate some matched pair data
set.seed(123)
n_pairs <- 50
x1 <- rnorm(n_pairs * 2)
x2 <- rbinom(n_pairs * 2, 1, 0.5)
strata <- rep(1:n_pairs, each = 2)
treatment <- rep(c(0, 1), n_pairs) # Case-Control
# ... (simulate response y) ...

# Fit the model
# Using Formula Interface
fit <- bclogit(y ~ x1 + x2, data = mydata, treatment = treatment, strata = strata)

# Summary of results
summary(fit)

# Coefficients
coef(fit)

# Confidence Intervals
confint(fit)
```


## Features

- **Standard Formula Interface**: Works like `glm` or `lm`.
- **Bayesian Prior Integration**: Uses `rstan` to fit models with informative priors.
- **Multiple Prior Types**: Supports "naive", "G prior", "PMP", and "hybrid".
- **Concordant Pairs Method**: Supports "GLM", "GEE", and "GLMM" for the concordant pairs step.
- We also implemented vanilla `clogit` which runs with a 5–10x speedup over `survival::clogit`. If you only want inference on the `j`th variable, you can set `do_inference_on_var = j` which gives the best speedup at large n: 9–11x at 10k–50k pairs.

## Citation

See [`CITATION.cff`](CITATION.cff), or in R: `citation("bclogit")`.

## License

This package is licensed under GPL3.
