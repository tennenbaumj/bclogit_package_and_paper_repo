## Resubmission

This is an update of bclogit (1.1 on CRAN -> 1.1.1). It addresses CRAN's
notice that 'glmmTMB' is likely to be archived.

* 'glmmTMB' has been moved from Imports to Suggests.
* 'glmmTMB' is used only for the optional `concordant_method = "GLMM"`, and
  only when it is available (checked with `requireNamespace()`). All other
  functionality (the default `"GLM"` and the `"GEE"` concordant methods, all
  prior types, and `clogit()`) works without it.
* If a user requests `concordant_method = "GLMM"` and 'glmmTMB' is not
  installed, bclogit stops with an informative error that says where to get
  it:

      install.packages("glmmTMB",
        repos = c("https://glmmtmb.r-universe.dev", "https://cloud.r-project.org"))

  This command works whether or not 'glmmTMB' is on CRAN. The same
  instructions are in `?bclogit` (argument `concordant_method`) and the
  README.
* Tests that use 'glmmTMB' are skipped when it is not installed.
* The minimum R version is now 4.1.0, the lowest that the current CRAN
  versions of bclogit's dependencies support.

## Test environments

* local Linux (WSL2), R 4.7.0

## R CMD check results

0 errors | 0 warnings | 0 notes
