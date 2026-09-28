library(testthat)
library(bclogit)

context("Test all combinations of concordant_method and prior_type")

# Use a small number of pairs to keep Stan sampling fast for tests
N_PAIRS <- 40
P <- 2

test_that("All combinations of concordant_method and prior_type work", {
  skip_if_no_stan()

  concordant_methods <- c("GLM", "GEE", if (requireNamespace("glmmTMB", quietly = TRUE)) "GLMM")
  prior_types <- c("Naive", "G prior", "PMP", "Hybrid")

  dat <- generate_test_data(n_pairs = N_PAIRS, p = P, seed = 123)

  for (cm in concordant_methods) {
    for (pt in prior_types) {
      msg <- paste0("Testing concordant_method = ", cm, ", prior_type = ", pt)
      # message(msg) # testthat captures message

      # We use chains = 1, iter = 500 to keep it relatively fast
      # but enough to see if it crashes or fails to converge
      fit <- expect_error(
        bclogit_default(
          y = dat$y,
          X = dat$X,
          treatment = dat$treatment,
          strata = dat$strata,
          concordant_method = cm,
          prior_type = pt,
          chains = 1,
          iter = 500,
          refresh = 0
        ),
        NA, # expect no error
        label = msg
      )

      expect_s3_class(fit, "bclogit")
      expect_true(!is.null(fit$coefficients), info = msg)
      expect_equal(length(fit$coefficients), P + 1, info = msg)
      expect_true(!is.null(fit$var), info = msg)
      expect_equal(dim(fit$var), c(P + 1, P + 1), info = msg)
      expect_true(!is.null(fit$matched_data), info = msg)
      expect_true(is.list(fit$matched_data), info = msg)
    }
  }
})

test_that("bclogit.formula works for all combinations", {
  skip_if_no_stan()

  # Just test a subset to keep it from being too slow, 
  # but enough to ensure formula interface works.
  # Or test all if it's not too bad.
  concordant_methods <- c("GLM", "GEE", if (requireNamespace("glmmTMB", quietly = TRUE)) "GLMM")
  prior_types <- c("Naive", "G prior", "PMP", "Hybrid")

  dat <- generate_test_data(n_pairs = N_PAIRS, p = P, seed = 456)
  df <- make_df(dat)

  # To save time, we might only test a few, but the request was "all"
  for (cm in concordant_methods) {
    for (pt in prior_types) {
      msg <- paste0("Testing bclogit formula: cm = ", cm, ", pt = ", pt)

      fit <- expect_error(
        bclogit(
          y ~ x1 + x2,
          data = df,
          treatment = treatment,
          strata = strata,
          concordant_method = cm,
          prior_type = pt,
          chains = 1,
          iter = 500,
          refresh = 0
        ),
        NA,
        label = msg
      )

      expect_s3_class(fit, "bclogit")
      expect_true(!is.null(fit$coefficients), info = msg)
      expect_true(!is.null(fit$matched_data), info = msg)
    }
  }
})
