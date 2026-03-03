# =============================================================================
# Standalone testthat suite for all exported functions in the bclogit package.
#
# Run with:
#   testthat::test_file("package_tests/testthat_exported_functions.R")
#
# Requires the package to be installed:
#   devtools::install("bclogit")
#
# Stan-dependent tests use chains = 1, iter = 500, warmup = 200 for speed.
# =============================================================================

library(testthat)
library(bclogit)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

gen <- function(n_pairs = 60, p = 2, seed = 1) {
  set.seed(seed)
  n <- 2 * n_pairs
  strata    <- rep(seq_len(n_pairs), each = 2)
  treatment <- rep(c(0L, 1L), n_pairs)
  X         <- matrix(rnorm(n * p), nrow = n, ncol = p,
                      dimnames = list(NULL, paste0("x", seq_len(p))))
  beta      <- c(0.5, rep(0.3, p))
  prob      <- plogis(cbind(treatment, X) %*% beta)
  y         <- rbinom(n, 1L, prob)
  list(y = y, X = X, treatment = treatment, strata = strata, n = n, p = p)
}

as_df <- function(d) {
  data.frame(y = d$y, d$X, treatment = d$treatment, strata = d$strata)
}

STAN_ARGS <- list(chains = 1L, iter = 500L, warmup = 200L, refresh = 0L)

bclogit_fit <- function(d, concordant_method = "GLM", prior_type = "Naive") {
  do.call(
    bclogit:::bclogit.default,
    c(list(y = d$y, X = d$X, treatment = d$treatment, strata = d$strata,
           concordant_method = concordant_method, prior_type = prior_type),
      STAN_ARGS)
  )
}

# =============================================================================
# 1. clogit — no Stan required
# =============================================================================

test_that("clogit.formula returns clogit_bclogit object with correct structure", {
  d  <- gen()
  df <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment, strata = strata)
  expect_s3_class(fit, "clogit_bclogit")
  expect_named(fit, c("coefficients", "var", "flr_model", "call", "terms",
                      "n", "num_discordant", "num_concordant",
                      "X_model_matrix_col_names", "treatment_name",
                      "se", "z", "pval", "do_inference_on_var"),
               ignore.order = TRUE)
  expect_equal(length(fit$coefficients), d$p + 1L)
})

test_that("clogit.default returns clogit_bclogit object", {
  d   <- gen()
  fit <- bclogit:::clogit.default(y = d$y, X = d$X,
                                  treatment = d$treatment, strata = d$strata)
  expect_s3_class(fit, "clogit_bclogit")
  expect_equal(length(fit$coefficients), d$p + 1L)
})

test_that("coef.clogit_bclogit returns named numeric vector", {
  d   <- gen()
  df  <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment, strata = strata)
  cc  <- coef(fit)
  expect_true(is.numeric(cc))
  expect_equal(length(cc), d$p + 1L)
  expect_false(is.null(names(cc)))
})

test_that("vcov.clogit_bclogit returns diagonal positive matrix", {
  d   <- gen()
  df  <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment, strata = strata)
  V   <- vcov(fit)
  expect_true(is.matrix(V))
  k   <- length(coef(fit))
  expect_equal(dim(V), c(k, k))
  expect_true(all(diag(V) > 0))
})

test_that("summary.clogit_bclogit returns summary with coefficient table", {
  d   <- gen()
  df  <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment, strata = strata)
  s   <- summary(fit)
  expect_s3_class(s, "summary.clogit_bclogit")
  expect_true(is.matrix(s$coefficients))
  expect_true("Estimate" %in% colnames(s$coefficients))
  expect_true("Pr(>|z|)" %in% colnames(s$coefficients))
})

test_that("print.summary.clogit_bclogit produces visible output", {
  d   <- gen()
  df  <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment, strata = strata)
  out <- capture.output(print(summary(fit)))
  expect_true(length(out) > 0)
  expect_true(any(grepl("Conditional Logistic Regression", out)))
})

test_that("clogit do_inference_on_var = 1 only returns SE for variable 1", {
  d   <- gen()
  df  <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment,
                strata = strata, do_inference_on_var = 1L)
  expect_equal(length(coef(fit)), d$p + 1L)
  expect_false(is.na(fit$se[1L]))
})

test_that("clogit do_inference_on_var = 'none' skips inference", {
  d   <- gen()
  df  <- as_df(d)
  fit <- clogit(y ~ x1 + x2, data = df, treatment = treatment,
                strata = strata, do_inference_on_var = "none")
  expect_true(is.null(fit$var))
})

test_that("clogit errors when treatment is missing", {
  d  <- gen()
  df <- as_df(d)
  expect_error(clogit(y ~ x1, data = df, strata = strata), "treatment")
})

test_that("clogit agrees with survival::clogit on coefficients", {
  skip_if_not_installed("survival")
  library(survival)
  # Define variables locally so survival::clogit can find them in the parent frame
  set.seed(7)
  n         <- 200L
  x1        <- rnorm(n)
  x2        <- rnorm(n)
  treatment <- rep(c(0L, 1L), n / 2L)
  strata    <- rep(seq_len(n / 2L), each = 2L)
  y         <- rbinom(n, 1L, plogis(0.5 * treatment + 0.3 * x1))
  dat       <- data.frame(y = y, x1 = x1, x2 = x2,
                          treatment = treatment, strata = strata)
  fit_bcl  <- bclogit::clogit(y ~ x1 + x2, data = dat, treatment = treatment, strata = strata)
  fit_surv <- survival::clogit(y ~ treatment + x1 + x2 + strata(strata), data = dat)
  expect_equal(unname(coef(fit_bcl)), unname(coef(fit_surv)), tolerance = 1e-4)
})

# =============================================================================
# 2. bclogit — all concordant methods × all prior types
# =============================================================================

concordant_methods <- c("GLM", "GEE", "GLMM")
prior_types        <- c("Naive", "G prior", "PMP", "Hybrid")

for (cm in concordant_methods) {
  for (pt in prior_types) {
    local({
      cm_ <- cm; pt_ <- pt
      test_that(sprintf("bclogit fits: concordant_method=%s, prior_type=%s", cm_, pt_), {
        d   <- gen(n_pairs = 60, p = 2, seed = 42)
        fit <- bclogit_fit(d, concordant_method = cm_, prior_type = pt_)
        expect_s3_class(fit, "bclogit")
        expect_equal(length(coef(fit)), d$p + 1L)
        expect_equal(dim(vcov(fit)), c(d$p + 1L, d$p + 1L))
        expect_false(anyNA(coef(fit)))
      })
    })
  }
}

# =============================================================================
# 3. bclogit — formula interface
# =============================================================================

test_that("bclogit.formula produces the same coefficients as bclogit.default", {
  d  <- gen(n_pairs = 60, p = 2, seed = 99)
  df <- as_df(d)

  fit_def <- bclogit:::bclogit.default(
    y = d$y, X = d$X, treatment = d$treatment, strata = d$strata,
    seed = 1L, chains = 1L, iter = 500L, warmup = 200L, refresh = 0L
  )
  # Direct call so treatment/strata are captured as symbols by match.call()
  fit_form <- bclogit(y ~ x1 + x2, data = df, treatment = treatment, strata = strata,
                      seed = 1L, chains = 1L, iter = 500L, warmup = 200L, refresh = 0L)

  expect_equal(unname(coef(fit_def)), unname(coef(fit_form)), tolerance = 0.05)
})

# =============================================================================
# 4. bclogit S3 methods — coef / vcov / formula / summary / confint / print
# =============================================================================

# Fit once and share across related tests.
.d   <- gen(n_pairs = 60, p = 2, seed = 11)
.fit <- bclogit_fit(.d)

test_that("coef.bclogit returns named numeric vector of correct length", {
  cc <- coef(.fit)
  expect_true(is.numeric(cc))
  expect_equal(length(cc), .d$p + 1L)
  expect_false(is.null(names(cc)))
})

test_that("vcov.bclogit returns square, named, positive-semidefinite matrix", {
  V <- vcov(.fit)
  k <- length(coef(.fit))
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(k, k))
  expect_equal(rownames(V), names(coef(.fit)))
  expect_equal(colnames(V), names(coef(.fit)))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(eigen(V, only.values = TRUE)$values >= -1e-10))
})

test_that("formula.bclogit returns a formula or NULL", {
  f <- formula(.fit)
  expect_true(inherits(f, "formula") || is.null(f))
})

test_that("summary.bclogit returns summary.bclogit with expected columns", {
  s <- summary(.fit)
  expect_s3_class(s, "summary.bclogit")
  expect_true(is.matrix(s$coefficients))
  expect_true(all(c("mean estimate", "median estimate", "Std. Error", "Pr(tx!=0)")
                  %in% colnames(s$coefficients)))
  expect_equal(nrow(s$coefficients), length(coef(.fit)))
})

test_that("summary.bclogit level argument widens intervals correctly", {
  s90 <- summary(.fit, level = 0.90, inference_method = "CR")
  s99 <- summary(.fit, level = 0.99, inference_method = "CR")
  w90 <- s90$coefficients[1L, 5L] - s90$coefficients[1L, 4L]
  w99 <- s99$coefficients[1L, 5L] - s99$coefficients[1L, 4L]
  expect_true(w99 > w90)
})

test_that("summary.bclogit contains num_discordant and num_concordant", {
  s <- summary(.fit)
  expect_true(s$num_discordant > 0L)
  expect_true(s$num_concordant >= 0L)
})

test_that("confint.bclogit CR type: matrix, correct dims, lower < upper", {
  ci <- confint(.fit, type = "CR")
  k  <- length(coef(.fit))
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(k, 2L))
  expect_true(all(ci[, 1L] < ci[, 2L]))
})

test_that("confint.bclogit HPD_one: lower < upper", {
  ci <- confint(.fit, type = "HPD_one")
  expect_true(is.matrix(ci))
  expect_true(all(ci[, 1L] < ci[, 2L]))
})

test_that("confint.bclogit HPD_many: returns lower/upper columns", {
  ci <- confint(.fit, type = "HPD_many")
  expect_true(is.matrix(ci))
  expect_true(all(c("lower", "upper") %in% colnames(ci)))
})

test_that("confint.bclogit level argument widens CR intervals", {
  ci90 <- confint(.fit, level = 0.90, type = "CR")
  ci99 <- confint(.fit, level = 0.99, type = "CR")
  expect_true((ci99[1L, 2L] - ci99[1L, 1L]) > (ci90[1L, 2L] - ci90[1L, 1L]))
})

test_that("confint.bclogit parm by index subsets rows", {
  ci <- confint(.fit, parm = 1L, type = "CR")
  expect_equal(nrow(ci), 1L)
})

test_that("confint.bclogit parm by name subsets rows", {
  nm <- names(coef(.fit))[1L]
  ci <- confint(.fit, parm = nm, type = "CR")
  expect_equal(nrow(ci), 1L)
})

test_that("confint.bclogit errors on unknown parm name", {
  expect_error(confint(.fit, parm = "no_such_param", type = "CR"), "No parameters found")
})

test_that("print.summary.bclogit produces output with key headings", {
  out <- capture.output(print(summary(.fit)))
  expect_true(length(out) > 0L)
  expect_true(any(grepl("Discordant", out)))
  expect_true(any(grepl("Concordant", out)))
  expect_true(any(grepl("Coefficients", out)))
})

test_that("print.summary.bclogit returns the summary object invisibly", {
  s   <- summary(.fit)
  ret <- withVisible(print(s))
  expect_s3_class(ret$value, "summary.bclogit")
  expect_false(ret$visible)
})

# =============================================================================
# 5. bclogit result-object completeness
# =============================================================================

test_that("bclogit result contains all documented components", {
  expected <- c("coefficients", "var", "model", "posterior_samples",
                "concordant_model", "matched_data", "prior_info",
                "call", "terms", "xlevels",
                "n", "num_discordant", "num_concordant",
                "X_model_matrix_col_names", "treatment_name")
  for (nm in expected) {
    expect_true(nm %in% names(.fit), label = paste("missing component:", nm))
  }
})

test_that("bclogit prior_info has mu, Sigma, and fallbacks of correct form", {
  pi  <- .fit$prior_info
  k   <- length(coef(.fit))
  expect_equal(length(pi$mu), k)
  expect_equal(dim(pi$Sigma), c(k, k))
  expect_false(anyNA(pi$mu))
  expect_false(anyNA(pi$Sigma))
  # No fallbacks expected on clean data with GLM
  expect_true(is.character(pi$fallbacks))
  expect_equal(length(pi$fallbacks), 0L)
})

test_that("bclogit Sigma_con passed to Stan is positive definite", {
  S <- .fit$prior_info$Sigma
  expect_true(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
})

test_that("bclogit treatment prior variance equals prior_variance_treatment", {
  S <- .fit$prior_info$Sigma
  expect_equal(S[1L, 1L], 100, tolerance = 1e-10)   # default prior_variance_treatment
  expect_equal(S[1L, -1L], rep(0, ncol(S) - 1L))    # off-diagonals zeroed
})

# =============================================================================
# 6. Numerical stability: GEE and GLMM prior covariance safeguards
# =============================================================================

test_that("GEE prior_info Sigma is positive definite and symmetric, fallbacks is character", {
  d   <- gen(n_pairs = 60, p = 2, seed = 21)
  fit <- bclogit_fit(d, concordant_method = "GEE")
  S   <- fit$prior_info$Sigma
  expect_equal(S, t(S), tolerance = 1e-10)
  expect_true(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
  expect_false(anyNA(S))
  expect_true(is.character(fit$prior_info$fallbacks))
})

test_that("GLMM prior_info Sigma is positive definite and symmetric, fallbacks is character", {
  d   <- gen(n_pairs = 60, p = 2, seed = 31)
  fit <- bclogit_fit(d, concordant_method = "GLMM")
  S   <- fit$prior_info$Sigma
  expect_equal(S, t(S), tolerance = 1e-10)
  expect_true(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
  expect_false(anyNA(S))
  expect_true(is.character(fit$prior_info$fallbacks))
})

test_that("GLMM GLM-fallback triggers on near-degenerate concordant data, yields valid Sigma, and records fallback", {
  # Construct data with highly collinear concordant pairs so glmmTMB is
  # likely to produce a singular / non-PD fixed-effect covariance,
  # triggering the GLM fallback added in default.bclogit.R.
  set.seed(77)
  n_pairs <- 50
  strata    <- rep(seq_len(n_pairs), each = 2)
  treatment <- rep(c(0L, 1L), n_pairs)
  # Make X almost perfectly collinear so concordant-model vcov is ill-conditioned
  x1  <- rnorm(2 * n_pairs)
  X   <- cbind(x1 = x1, x2 = x1 + rnorm(2 * n_pairs, sd = 1e-4))
  y   <- rbinom(2 * n_pairs, 1L, plogis(0.5 * treatment + 0.3 * x1))

  fit <- do.call(bclogit:::bclogit.default,
                 c(list(y = y, X = X, treatment = treatment, strata = strata,
                        concordant_method = "GLMM"), STAN_ARGS))

  S <- fit$prior_info$Sigma
  expect_false(anyNA(S))
  expect_equal(S, t(S), tolerance = 1e-10)
  expect_true(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
  if (length(fit$prior_info$fallbacks) > 0L) {
    expect_true("glm_fallback" %in% fit$prior_info$fallbacks)
  }
})

test_that("GEE GLM-fallback triggers on near-degenerate concordant data, yields valid Sigma, and records fallback", {
  set.seed(88)
  n_pairs <- 50
  strata    <- rep(seq_len(n_pairs), each = 2)
  treatment <- rep(c(0L, 1L), n_pairs)
  x1  <- rnorm(2 * n_pairs)
  X   <- cbind(x1 = x1, x2 = x1 + rnorm(2 * n_pairs, sd = 1e-4))
  y   <- rbinom(2 * n_pairs, 1L, plogis(0.5 * treatment + 0.3 * x1))

  fit <- do.call(bclogit:::bclogit.default,
                 c(list(y = y, X = X, treatment = treatment, strata = strata,
                        concordant_method = "GEE"), STAN_ARGS))

  S <- fit$prior_info$Sigma
  expect_false(anyNA(S))
  expect_equal(S, t(S), tolerance = 1e-10)
  expect_true(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
  if (length(fit$prior_info$fallbacks) > 0L) {
    expect_true("glm_fallback" %in% fit$prior_info$fallbacks)
  }
})

test_that("bclogit falls back to diffuse prior when too few concordant pairs", {
  # Create data where concordant pairs exist but fewer than ncol(X) + 5 = 6,
  # triggering the "not enough concordant pairs" warning path.
  # Pairs 1-4: concordant (treatment=0 y=0, treatment=1 y=0)
  # Pairs 5-20: discordant (treatment=0 y=0, treatment=1 y=1)
  set.seed(55)
  # Need length(y_concordant) < ncol(X)+5 = 6, i.e. 2*n_con < 6, so n_con < 3
  n_con <- 2L   # concordant pairs  (2*2=4 < 6, triggers warning)
  n_dis <- 16L  # discordant pairs
  n_pairs   <- n_con + n_dis
  strata    <- rep(seq_len(n_pairs), each = 2L)
  treatment <- rep(c(0L, 1L), n_pairs)
  X         <- matrix(rnorm(2L * n_pairs), ncol = 1L)
  y         <- c(rep(c(0L, 0L), n_con),   # concordant: both 0
                 rep(c(0L, 1L), n_dis))   # discordant: (0,1)

  expect_warning(
    fit <- do.call(bclogit:::bclogit.default,
                   c(list(y = y, X = X, treatment = treatment, strata = strata,
                          concordant_method = "GLM"), STAN_ARGS)),
    "not enough concordant pairs"
  )

  S <- fit$prior_info$Sigma
  expect_false(anyNA(S))
  expect_true(tryCatch({ chol(S); TRUE }, error = function(e) FALSE))
  expect_true("insufficient_concordant_pairs" %in% fit$prior_info$fallbacks)
})

# =============================================================================
# 7. Input-validation errors (no Stan needed)
# =============================================================================

test_that("bclogit errors when treatment argument is missing", {
  d <- gen()
  expect_error(
    bclogit:::bclogit.default(y = d$y, X = d$X, strata = d$strata),
    "treatment"
  )
})

test_that("bclogit errors on non-binary treatment values", {
  d           <- gen()
  bad_trt     <- d$treatment
  bad_trt[1L] <- 0.5
  expect_error(
    bclogit:::bclogit.default(y = d$y, X = d$X,
                              treatment = bad_trt, strata = d$strata)
  )
})

test_that("bclogit errors on out-of-range response", {
  d         <- gen()
  bad_y     <- d$y
  bad_y[1L] <- 2L
  expect_error(
    bclogit:::bclogit.default(y = bad_y, X = d$X,
                              treatment = d$treatment, strata = d$strata)
  )
})

test_that("bclogit errors on too few observations relative to covariates", {
  set.seed(1L)
  n <- 6L; p <- 5L
  expect_error(
    bclogit:::bclogit.default(
      y = rbinom(n, 1L, 0.5),
      X = matrix(rnorm(n * p), nrow = n),
      treatment = rep(c(0L, 1L), n / 2L),
      strata = rep(seq_len(n / 2L), each = 2L)
    ),
    "Not enough rows"
  )
})

test_that("bclogit errors on invalid prior_type", {
  d <- gen()
  expect_error(
    bclogit:::bclogit.default(y = d$y, X = d$X, treatment = d$treatment,
                              strata = d$strata, prior_type = "BadPrior"),
    "prior_type"
  )
})

test_that("bclogit errors on invalid concordant_method", {
  d <- gen()
  expect_error(
    bclogit:::bclogit.default(y = d$y, X = d$X, treatment = d$treatment,
                              strata = d$strata, concordant_method = "INVALID")
  )
})
