# bclogit 1.1.1

* Moved `glmmTMB` from Imports to Suggests. `concordant_method = "GLMM"` now requires `glmmTMB` to be installed and gives an informative error otherwise; if `glmmTMB` is not available on CRAN it can be installed from <https://glmmtmb.r-universe.dev>.
* Now requires R >= 4.1.0, the minimum that bclogit's dependencies (via rstan) support on current CRAN.
* During GEE and GLMM prior construction, returns an array of robustness fallbacks used.

# bclogit 1.1

* Fixed bug in naive prior due to array naming.
* Implemented more robustness checks for GLMM.
* Documented the robustness checks for both GEE and GLMM.
* Improved other functions' documentation.

# bclogit 1.0

* Initial release.
