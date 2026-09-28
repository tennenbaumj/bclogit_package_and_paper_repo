#' @importFrom methods new
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel CxxFlags
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("Welcome to bclogit v",
            utils::packageVersion("bclogit"),
            sep = ""
        )
    )
}

# glmmTMB is in Suggests (it may be archived on CRAN); wrapped so tests can mock it
glmmTMB_available <- function() {
    requireNamespace("glmmTMB", quietly = TRUE)
}
