#' Initialize a new bclogit model
#'
#' This function fits a Bayesian conditional logistic regression model, incorporating
#' information from concordant pairs to improve estimation.
#'
#' @param formula For the formula method, a symbolic description of the model to be fitted. 
#' @param y For the default method, a binary (0,1) vector containing the response of each subject.
#' @param data A data.frame, data.table, or model.matrix containing the variables (optional for formula method).
#' @param treatment Optional vector specifying the treatment variable (required for default method, or can be specified in formula method).
#' @param strata Vector specifying the strata (matched pairs).
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action A function which indicates what should happen when the data contain NAs. 
#' @param X A data.frame, data.table, or model.matrix containing the variables (optional for formula method, required for default method).
#' @param concordant_method The method to use for fitting the concordant pairs and reservoir. Options are "GLM", "GEE", and "GLMM".
#' @param prior_type The type of prior to use for the discordant pairs. Options are "Naive", "G prior", "PMP", and "Hybrid".
#' @param chains Number of chains for Stan sampling. Default is 4.
#' @param treatment_name Optional string name for the treatment variable.
#' @param return_raw_stan_output Logical; if \code{TRUE}, the raw Stan posterior samples
#'   (iterations x chains x parameters) are stored in the returned object. Default \code{FALSE}.
#' @param prior_variance_treatment Prior variance for the treatment coefficient in the
#'   covariance matrix \code{Sigma_con}. Default is 100.
#' @param stan_refresh How often Stan reports sampling progress (in iterations).
#'   Default is 0 (silent). Set to a positive integer (e.g., 1 or 100) to see progress.
#' @param call Optional call object to store in the result.
#' @param ... Additional arguments passed to \code{rstan::sampling} (e.g., \code{iter}, \code{warmup}, \code{thin}, \code{seed}, \code{control}).
#' @return A list of class `bclogit` containing:
#'   \item{coefficients}{Estimated coefficients (posterior means).}
#'   \item{var}{Variance-covariance matrix of coefficients.}
#'   \item{model}{The fitted Stan model object.}
#'   \item{posterior_samples}{Raw posterior samples as a 3D array (iterations x chains x parameters) from \code{rstan::extract(model, permuted = FALSE)}. Only populated when \code{return_raw_stan_output = TRUE}; \code{NULL} otherwise.}
#'   \item{concordant_model}{The fitted model object for the concordant pairs/reservoir (GLM/GEE/GLMM).}
#'   \item{matched_data}{The processed matched pairs data from the premodeling step.}
#'   \item{prior_info}{A list with elements \code{mu} (prior mean vector), \code{Sigma}
#'     (prior covariance matrix), and \code{fallbacks} (a character vector naming every
#'     robustness fallback that was triggered during prior construction, or
#'     \code{character(0)} if none were needed). Possible fallback names are:
#'     \code{"insufficient_concordant_pairs"} (too few concordant pairs; diffuse prior used),
#'     \code{"na_in_prior_mean"} (NAs in prior mean replaced with 0),
#'     \code{"na_in_prior_covariance"} (NAs in prior covariance replaced with \code{diag(100)}),
#'     \code{"glm_fallback"} (non-positive-definite GEE/GLMM covariance; re-fitted via GLM),
#'     \code{"diffuse_prior"} (non-positive-definite covariance replaced with \code{diag(100)}).
#'   }
#'   \item{call}{The function call.}
#'   \item{terms}{The model terms.}
#'   \item{num_discordant}{Number of discordant pairs used.}
#'   \item{num_concordant}{Number of concordant pairs/reservoir entries used.}
#' @section Numerical Stability in GEE Pre-modeling:
#' When \code{concordant_method = "GEE"}, the function applies a layered set of
#' safeguards to ensure the prior covariance \code{Sigma_con} passed to Stan is
#' valid (positive definite, symmetric, and free of \code{NA} values).
#'
#' \describe{
#'   \item{Naive (model-based) covariance}{
#'     GEE's default covariance is the sandwich (robust) estimator, which can
#'     produce degenerate (e.g. negative) variances when the number of independent
#'     clusters is small. The code instead uses \code{geese$vbeta.naiv}, the
#'     model-based covariance that assumes the working correlation structure is
#'     exactly correct. It is less efficient asymptotically but always
#'     well-structured for small samples.
#'   }
#'   \item{Cholesky positive-definiteness check}{
#'     After extracting the covariance, a Cholesky decomposition is attempted
#'     via \code{tryCatch(chol(Sigma_con), ...)}. This is the canonical test for
#'     positive definiteness: it succeeds if and only if the matrix is PD. A
#'     failed decomposition sets an \code{is_pd = FALSE} flag without crashing
#'     the call.
#'   }
#'   \item{GLM fallback for non-PD GEE or GLMM covariance}{
#'     If the GEE naive covariance or the GLMM fixed-effect \code{vcov} is still
#'     not positive definite (GLMM failures are common with boundary variance
#'     estimates or convergence issues), a standard
#'     \code{glm(..., family = binomial)} is re-fitted on the same concordant
#'     data and its Fisher-information-based \code{vcov} is used instead. That
#'     matrix is always PD for a non-degenerate design matrix, providing a
#'     graceful degradation from GEE/GLMM to GLM.
#'   }
#'   \item{Symmetrization}{
#'     \code{Sigma_con <- (Sigma_con + t(Sigma_con)) / 2} is applied after
#'     the GLM fallback re-extraction and again unconditionally before the matrix
#'     is passed to Stan. Name-indexed covariance matching can introduce tiny
#'     floating-point asymmetries; averaging with the transpose enforces exact
#'     symmetry as required by Stan's Cholesky-based samplers.
#'   }
#'   \item{NA sanitization}{
#'     Aliased model terms (collinear columns in the concordant data) cause
#'     \code{coef()} and \code{vcov()} to return \code{NA}. These are caught
#'     before Stan sees them: \code{NA} entries in the prior mean vector default
#'     to \code{0} (a neutral prior), and any \code{NA} in the covariance matrix
#'     triggers replacement of the entire matrix with a wide diagonal
#'     (\code{diag(100, p)}), an uninformative but valid prior.
#'   }
#'   \item{Diffuse prior fallback}{
#'     Multiple failure paths (too few concordant pairs, persistent non-PD
#'     matrix, \code{NA}-contaminated covariance) fall back to
#'     \code{diag(100, p)}, which is positive definite by construction and
#'     encodes a broad, independent prior that lets the discordant-pair
#'     likelihood dominate.
#'   }
#'   \item{Treatment prior reset}{
#'     After extracting the concordant-model covariance, the entire first row
#'     and column (corresponding to the treatment coefficient) are zeroed and
#'     the diagonal entry is set to \code{prior_variance_treatment} (default
#'     100). This decouples the treatment prior from the nuisance-covariate
#'     prior: concordant pairs inform covariate shrinkage but are not allowed to
#'     shrink the treatment coefficient, which remains independently diffuse.
#'     Zeroing the off-diagonals also eliminates any spurious prior correlation
#'     between treatment and covariates induced by the concordant-model fit.
#'   }
#' }
#'
#' @seealso \code{\link{summary.bclogit}}, \code{\link{confint.bclogit}},
#' \code{\link{vcov.bclogit}}, \code{\link{coef.bclogit}}
#' @examples
#' \donttest{
#' # Example usage
#' data("fhs")
#' fit <- bclogit(PREVHYP ~ TOTCHOL + CIGPDAY + BMI + HEARTRTE, 
#'   data = fhs, treatment = PERIOD, strata = RANDID)
#' summary(fit)
#' }
#' @export
bclogit <- function(formula, data, treatment = NULL, strata = NULL,
                    subset = NULL, na.action = NULL,
                    concordant_method = "GLM", prior_type = "Naive",
                    chains = 4, return_raw_stan_output = FALSE,
                    prior_variance_treatment = 100, stan_refresh = 0, ...) {
    UseMethod("bclogit")
}
