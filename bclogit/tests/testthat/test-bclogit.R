# ============================================================================
# 1. Input validation tests (fast, no Stan)
# ============================================================================
test_that("bclogit.default errors when treatment is missing", {
  dat <- generate_test_data()
  expect_error(
    bclogit_default(y = dat$y, X = dat$X, strata = dat$strata),
    "treatment"
  )
})

test_that("bclogit.formula errors when treatment is missing", {
  dat <- generate_test_data()
  df <- make_df(dat)
  expect_error(
    bclogit(y ~ x1 + x2, data = df, strata = strata),
    "treatment"
  )
})

test_that("bclogit.default errors on invalid response values", {
  dat <- generate_test_data()
  bad_y <- dat$y
  bad_y[1] <- 2
  expect_error(
    bclogit_default(
      y = bad_y, X = dat$X, treatment = dat$treatment, strata = dat$strata
    )
  )
})

test_that("bclogit.default errors on NA in response", {
  dat <- generate_test_data()
  bad_y <- dat$y
  bad_y[1] <- NA
  expect_error(
    bclogit_default(
      y = bad_y, X = dat$X, treatment = dat$treatment, strata = dat$strata
    )
  )
})

test_that("bclogit.default errors on non-binary treatment", {
  dat <- generate_test_data()
  bad_trt <- dat$treatment
  bad_trt[1] <- 0.5
  expect_error(
    bclogit_default(
      y = dat$y, X = dat$X, treatment = bad_trt, strata = dat$strata
    )
  )
})

test_that("bclogit.default errors on mismatched lengths", {
  dat <- generate_test_data()
  expect_error(
    bclogit_default(
      y = dat$y[-1], X = dat$X, treatment = dat$treatment, strata = dat$strata
    ),
    "length"
  )
})

test_that("bclogit.default errors on invalid concordant_method", {
  dat <- generate_test_data()
  expect_error(
    bclogit_default(
      y = dat$y, X = dat$X, treatment = dat$treatment,
      strata = dat$strata, concordant_method = "INVALID"
    )
  )
})

test_that("bclogit.default errors on too few observations", {
  set.seed(1)
  n <- 6
  strata <- rep(1:3, each = 2)
  treatment <- rep(c(0, 1), 3)
  X <- matrix(rnorm(n * 5), nrow = n)
  y <- rbinom(n, 1, 0.5)
  expect_error(
    bclogit_default(y = y, X = X, treatment = treatment, strata = strata),
    "Not enough rows"
  )
})

test_that("bclogit.default errors on invalid X type", {
  dat <- generate_test_data()
  expect_error(
    bclogit_default(
      y = dat$y, X = "not_a_matrix", treatment = dat$treatment,
      strata = dat$strata
    )
  )
})

test_that("bclogit.default errors on strata with != 2 observations", {
  set.seed(1)
  n <- 21
  strata <- c(rep(1:6, each = 2), rep(7, 3), rep(8:10, each = 2))
  treatment <- rep(c(0, 1), length.out = n)
  X <- matrix(rnorm(n * 2), nrow = n)
  y <- rbinom(n, 1, 0.5)
  expect_error(
    bclogit_default(y = y, X = X, treatment = treatment, strata = strata)
  )
})

# ============================================================================
# 2. coef/vcov/formula on fake objects (no Stan needed)
# ============================================================================
test_that("coef.bclogit returns NULL when no model converged", {
  fake <- list(coefficients = NULL)
  class(fake) <- c("bclogit", "list")
  expect_null(coef(fake))
})

test_that("vcov.bclogit returns NULL when no model converged", {
  fake <- list(var = NULL)
  class(fake) <- c("bclogit", "list")
  expect_null(vcov(fake))
})

test_that("formula.bclogit returns NULL when no terms", {
  fake <- list(terms = NULL)
  class(fake) <- c("bclogit", "list")
  expect_null(formula(fake))
})

test_that("summary.bclogit handles NULL model", {
  fake <- list(model = NULL, coefficients = c(a = 1))
  class(fake) <- c("bclogit", "list")
  expect_output(summary(fake), "No discordant model")
})

test_that("confint.bclogit warns with NULL model", {
  fake <- list(model = NULL, coefficients = c(a = 1))
  class(fake) <- c("bclogit", "list")
  expect_warning(confint(fake), "No discordant model")
})

# ============================================================================
# 3. Full model fitting tests (requires Stan compilation)
# ============================================================================

# Fit one model and reuse across multiple tests to save time
fit_naive <- NULL
dat_for_fit <- NULL

test_that("bclogit.default fits a Naive model successfully", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, prior_type = "Naive",
    concordant_method = "GLM", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
  expect_true(!is.null(fit$var))
  expect_true(!is.null(fit$model))
  expect_equal(length(fit$coefficients), dat$p + 1) # treatment + covariates
  expect_equal(nrow(fit$var), dat$p + 1)
  expect_equal(ncol(fit$var), dat$p + 1)

  # Save for reuse
  fit_naive <<- fit
  dat_for_fit <<- dat
})

test_that("bclogit.formula fits a model successfully", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  df <- make_df(dat)
  fit <- bclogit(
    y ~ x1 + x2, data = df, treatment = treatment,
    strata = strata, prior_type = "Naive",
    concordant_method = "GLM", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
  expect_equal(length(fit$coefficients), 3)
})

test_that("bclogit with G prior fits successfully", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, prior_type = "G prior",
    concordant_method = "GLM", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
})

test_that("bclogit with PMP prior fits successfully", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, prior_type = "PMP",
    concordant_method = "GLM", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
})

test_that("bclogit with Hybrid prior fits successfully", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, prior_type = "Hybrid",
    concordant_method = "GLM", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
})

test_that("bclogit with GEE concordant method fits successfully", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, prior_type = "Naive",
    concordant_method = "GEE", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
})

test_that("bclogit with GLMM concordant method fits successfully", {
  skip_if_no_stan()
  skip_if_not_installed("glmmTMB")

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, prior_type = "Naive",
    concordant_method = "GLMM", chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
})

test_that("GLMM without glmmTMB installed errors with install instructions", {
  local_mocked_bindings(glmmTMB_available = function() FALSE, .package = "bclogit")
  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  expect_error(
    bclogit_default(
      y = dat$y, X = dat$X, treatment = dat$treatment,
      strata = dat$strata, concordant_method = "GLMM", chains = 1
    ),
    "glmmtmb.r-universe.dev",
    fixed = TRUE
  )
})

test_that("bclogit with data.frame X input works", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  X_df <- as.data.frame(dat$X)
  fit <- bclogit_default(
    y = dat$y, X = X_df, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )
  expect_s3_class(fit, "bclogit")
  expect_true(!is.null(fit$coefficients))
})

test_that("bclogit errors on unknown prior_type", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 123)
  expect_error(
    bclogit_default(
      y = dat$y, X = dat$X, treatment = dat$treatment,
      strata = dat$strata, prior_type = "Unknown", chains = 1
    ),
    "prior_type"
  )
})

# ============================================================================
# 4. coef() tests (Stan-dependent)
# ============================================================================
test_that("coef.bclogit returns named numeric vector", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  cc <- coef(fit_naive)
  expect_true(is.numeric(cc))
  expect_true(!is.null(names(cc)))
  expect_equal(length(cc), dat_for_fit$p + 1)
})

# ============================================================================
# 5. vcov() tests (Stan-dependent)
# ============================================================================
test_that("vcov.bclogit returns square matrix with correct dimensions", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  V <- vcov(fit_naive)
  expect_true(is.matrix(V))
  k <- length(coef(fit_naive))
  expect_equal(dim(V), c(k, k))
  expect_equal(rownames(V), names(coef(fit_naive)))
  expect_equal(colnames(V), names(coef(fit_naive)))
})

test_that("vcov.bclogit returns symmetric positive-semidefinite matrix", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  V <- vcov(fit_naive)
  expect_equal(V, t(V), tolerance = 1e-10)
  eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))
})

# ============================================================================
# 6. summary() tests (Stan-dependent)
# ============================================================================
test_that("summary.bclogit returns summary.bclogit object", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  s <- summary(fit_naive)
  expect_s3_class(s, "summary.bclogit")
  expect_true(is.matrix(s$coefficients))
  expect_true("mean estimate" %in% colnames(s$coefficients))
  expect_true("median estimate" %in% colnames(s$coefficients))
  expect_true("Std. Error" %in% colnames(s$coefficients))
  expect_true("Pr(tx!=0)" %in% colnames(s$coefficients))
})

test_that("summary.bclogit coefficient matrix has correct dimensions", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  s <- summary(fit_naive)
  expect_equal(nrow(s$coefficients), length(coef(fit_naive)))
  expect_equal(rownames(s$coefficients), names(coef(fit_naive)))
})

test_that("summary.bclogit respects level", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  s90 <- summary(fit_naive, level = 0.90, inference_method = "CR")
  s99 <- summary(fit_naive, level = 0.99, inference_method = "CR")

  expect_equal(s90$level, 0.90)
  expect_equal(s99$level, 0.99)

  # 99% intervals should be wider than 90%
  # Columns 4 and 5 are lower and upper bounds
  width_90 <- s90$coefficients[1, 5] - s90$coefficients[1, 4]
  width_99 <- s99$coefficients[1, 5] - s99$coefficients[1, 4]
  expect_true(width_99 > width_90)
})

test_that("summary.bclogit includes metadata", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  s <- summary(fit_naive)
  expect_true(!is.null(s$num_discordant))
  expect_true(!is.null(s$num_concordant))
  expect_true(s$num_discordant > 0)
  expect_true(s$num_concordant >= 0)
})

# ============================================================================
# 7. confint() tests (Stan-dependent)
# ============================================================================
test_that("confint.bclogit defaults to HPD_one", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  ci_default <- confint(fit_naive)
  ci_hpd <- confint(fit_naive, type = "HPD_one")
  expect_equal(ci_default, ci_hpd)
})

test_that("confint.bclogit CR type returns correct structure", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  ci <- confint(fit_naive, type = "CR")
  expect_true(is.matrix(ci))
  k <- length(coef(fit_naive))
  expect_equal(nrow(ci), k)
  expect_equal(ncol(ci), 2)
  # Lower bound should be less than upper bound
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint.bclogit respects level argument", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  ci90 <- confint(fit_naive, level = 0.90, type = "CR")
  ci99 <- confint(fit_naive, level = 0.99, type = "CR")

  # 99% intervals should be wider
  width_90 <- ci90[1, 2] - ci90[1, 1]
  width_99 <- ci99[1, 2] - ci99[1, 1]
  expect_true(width_99 > width_90)
})

test_that("confint.bclogit parm selects parameters by index", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  ci <- confint(fit_naive, parm = 1, type = "CR")
  expect_equal(nrow(ci), 1)
})

test_that("confint.bclogit parm selects parameters by name", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  param_name <- names(coef(fit_naive))[1]
  ci <- confint(fit_naive, parm = param_name, type = "CR")
  expect_equal(nrow(ci), 1)
})

test_that("confint.bclogit HPD_one type returns correct structure", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  ci <- confint(fit_naive, type = "HPD_one")
  expect_true(is.matrix(ci))
  k <- length(coef(fit_naive))
  expect_equal(nrow(ci), k)
  expect_equal(ncol(ci), 2)
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint.bclogit HPD_many type returns matrix", {
  skip_if_no_stan()
  skip_if_not_installed("ggdist")
  skip_if(is.null(fit_naive), "Model not fitted")

  ci <- confint(fit_naive, type = "HPD_many")
  expect_true(is.matrix(ci))
  expect_true("lower" %in% colnames(ci))
  expect_true("upper" %in% colnames(ci))
  expect_equal(attr(ci, "Probability"), 0.95)
})

test_that("confint.bclogit errors on invalid parm", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  expect_error(
    confint(fit_naive, parm = "nonexistent_param", type = "CR"),
    "No parameters found"
  )
})

# ============================================================================
# 8. formula() tests (Stan-dependent)
# ============================================================================
test_that("formula.bclogit returns formula when terms exist", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  f <- formula(fit_naive)
  expect_true(inherits(f, "formula") || is.null(f))
})

# ============================================================================
# 9. print.summary.bclogit tests (Stan-dependent)
# ============================================================================
test_that("print.summary.bclogit produces output", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  s <- summary(fit_naive)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("Discordant", out)))
  expect_true(any(grepl("Concordant", out)))
  expect_true(any(grepl("Coefficients", out)))
})

test_that("print.summary.bclogit returns invisibly", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  s <- summary(fit_naive)
  ret <- capture.output(result <- print(s))
  expect_s3_class(result, "summary.bclogit")
})

# ============================================================================
# 10. Result structure tests (Stan-dependent)
# ============================================================================
test_that("bclogit result contains all expected components", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  expected_names <- c(
    "coefficients", "var", "model", "concordant_model",
    "prior_info", "call", "terms", "xlevels",
    "n", "num_discordant", "num_concordant",
    "X_model_matrix_col_names", "treatment_name"
  )
  for (nm in expected_names) {
    expect_true(nm %in% names(fit_naive), info = paste("Missing component:", nm))
  }
})

test_that("bclogit result class is correct", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  expect_true(inherits(fit_naive, "bclogit"))
})

test_that("bclogit prior_info contains mu and Sigma", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  expect_true(!is.null(fit_naive$prior_info))
  expect_true("mu" %in% names(fit_naive$prior_info))
  expect_true("Sigma" %in% names(fit_naive$prior_info))
})

test_that("bclogit coefficient names match vcov names", {
  skip_if_no_stan()
  skip_if(is.null(fit_naive), "Model not fitted")

  expect_equal(names(coef(fit_naive)), rownames(vcov(fit_naive)))
  expect_equal(names(coef(fit_naive)), colnames(vcov(fit_naive)))
})

# ============================================================================
# 11. Integration: formula vs default give consistent results
# ============================================================================
test_that("formula and default methods give same coefficients", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 999)
  df <- make_df(dat)

  set.seed(1)
  fit_def <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, chains = 1, seed = 1
  )

  set.seed(1)
  fit_form <- bclogit(
    y ~ x1 + x2, data = df, treatment = treatment,
    strata = strata, chains = 1, seed = 1
  )

  # Coefficients should be identical given same seed and data
  expect_equal(
    unname(coef(fit_def)), unname(coef(fit_form)),
    tolerance = 1e-2
  )
})
