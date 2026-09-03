# bclogit <img src="man/figures/logo.png" align="right" height="139" alt="bclogit hex logo" />

`bclogit` fits conditional logistic regression models for matched
case-control data. It can optionally incorporate information from
concordant pairs (and a reservoir of additional controls) as an
informative Bayesian prior, improving estimation for the discordant
pairs when data are sparse. It also ships a fast vanilla `clogit()`
drop-in replacement for `survival::clogit()` with a 5-11x speedup.

## Installation

Install the latest version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github(
  "tennenbaumj/bclogit_package_and_paper_repo",
  subdir = "bclogit"
)
```

The package compiles Stan models on installation, so a C++ toolchain
(Rtools on Windows, Xcode command line tools on macOS) is required.

## Usage

```r
library(bclogit)

fit <- bclogit(y ~ x1 + x2, data = mydata, treatment = treatment, strata = strata)
summary(fit)
coef(fit)
confint(fit)
```

See the [function reference](reference/index.html) for the full API.

## Citation

Kapelner, A. and Tennenbaum, J. (2026). "Improved Conditional Logistic
Regression using Information in Concordant Pairs with Software."
<https://doi.org/10.48550/arXiv.2602.08212>. In R: `citation("bclogit")`.
