#' @param subset An optional vector specifying a subset of observations.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param X A data.frame, data.table, or model.matrix containing the variables.
#' @return A list of class \code{"bclogit"} containing:
#'   \item{coefficients}{Estimated coefficients (posterior means).}
#'   \item{var}{Variance-covariance matrix of the coefficients (posterior covariance).}
#'   \item{model}{The fitted Stan model object for the discordant pairs.}
#'   \item{posterior_samples}{Raw posterior samples as a 3D array
#'     (iterations x chains x parameters) from \code{rstan::extract(model, permuted = FALSE)}.
#'     Only populated when \code{return_raw_stan_output = TRUE}; \code{NULL} otherwise.}
#'   \item{concordant_model}{The fitted model object for the concordant pairs (GLM, GEE, or GLMM).}
#'   \item{matched_data}{The processed matched pairs data from the C++ pre-modeling step.}
#'   \item{prior_info}{A list with elements \code{mu} (prior mean vector), \code{Sigma}
#'     (prior covariance matrix), and \code{fallbacks} (character vector of robustness
#'     fallbacks triggered; \code{character(0)} if none).}
#'   \item{call}{The function call.}
#'   \item{terms}{The model terms.}
#'   \item{xlevels}{Factor level information (always \code{NULL} for this method).}
#'   \item{n}{Total number of observations.}
#'   \item{num_discordant}{Number of discordant pairs used for fitting.}
#'   \item{num_concordant}{Number of concordant pairs used for the prior.}
#'   \item{X_model_matrix_col_names}{Column names of the covariate model matrix.}
#'   \item{treatment_name}{Name of the treatment variable.}
#' @describeIn bclogit Default method for matrix/data input.
#' @export
bclogit.default <- function(formula = NULL,
                            data = NULL,
                            treatment = NULL,
                            strata = NULL,
                            subset = NULL,
                            na.action = NULL,
                            concordant_method = "GLM",
                            prior_type = "Naive",
                            chains = 4,
                            return_raw_stan_output = FALSE,
                            prior_variance_treatment = 100,
                            stan_refresh = 0,
                            ...,
                            y = NULL,
                            X = NULL,
                            treatment_name = NULL,
                            call = NULL) {
  # ------------------------------------------------------------------------
  # 1. Input Validation and Pre-processing
  # ------------------------------------------------------------------------
  if (missing(treatment)) stop("The 'treatment' argument is required.")
  assertVector(subset, null.ok = TRUE)
  assertFunction(na.action, null.ok = TRUE)
  assertChoice(prior_type, c("Naive", "G prior", "PMP", "Hybrid"))
  assertCount(chains, positive = TRUE)
  assertString(treatment_name, null.ok = TRUE)

  # Use copy of X if needed or just use X directly
  if (test_data_frame(X, types = c("numeric", "integer", "factor", "logical"))) {
    data_mat <- model.matrix(~., data = X)
    if ("(Intercept)" %in% colnames(data_mat)) {
      data_mat <- data_mat[, colnames(data_mat) != "(Intercept)", drop = FALSE]
    }
  } else if (test_matrix(X, mode = "numeric")) {
    if (sd(X[, 1]) == 0) {
      X <- X[, -1, drop = FALSE]
    }
    data_mat <- as.matrix(X)
  } else {
    assert(
      check_data_frame(X),
      check_matrix(X, mode = "numeric")
    )
  }

  n <- nrow(data_mat)
  assertChoice(concordant_method, c("GLM", "GEE", "GLMM"))

  # Relaxed assertions to allow factors/characters, then convert
  assertVector(y, any.missing = FALSE, len = n)
  if (is.factor(y)) {
    y <- as.numeric(as.character(y))
  } else {
    y <- as.numeric(y)
  }
  assertNumeric(y, lower = 0, upper = 1, any.missing = FALSE, len = n)

  assertVector(treatment, any.missing = FALSE, len = n)
  if (is.factor(treatment)) {
    treatment <- as.numeric(as.character(treatment))
  } else {
    treatment <- as.numeric(treatment)
  }
  assertNumeric(treatment, lower = 0, upper = 1, any.missing = FALSE, len = n)

  if (is.factor(strata) || is.character(strata)) {
    strata <- as.numeric(as.factor(strata))
  }
  assertNumeric(strata, any.missing = FALSE, len = n)

  if (!all(treatment %in% c(0, 1))) {
    stop("Treatment must be binary 0 or 1.")
  }
  if (is.null(treatment_name)) {
    orig_trt <- substitute(treatment)
    if (is.name(orig_trt)) {
      treatment_name <- as.character(orig_trt)
    } else {
      treatment_name <- paste(deparse(orig_trt), collapse = "")
      if (nchar(treatment_name) > 30 || grepl("^c\\(", treatment_name) || grepl("^[0-9\\.\\-]+", treatment_name)) {
        treatment_name <- "treatment"
      }
    }
  }

  # Filter to keep only strata with exactly 2 observations
  strata_counts <- table(strata)
  valid_strata <- names(strata_counts)[strata_counts == 2]
  keep_idx <- strata %in% as.numeric(valid_strata)

  if (!all(keep_idx)) {
    y <- y[keep_idx]
    data_mat <- data_mat[keep_idx, , drop = FALSE]
    treatment <- treatment[keep_idx]
    strata <- strata[keep_idx]
    n <- length(y)
  }

  # Capture terms for summary/printing usage if possible (though we used matrix for fitting)
  model_terms <- tryCatch(terms(y ~ ., data = as.data.frame(data_mat)), error = function(e) NULL)

  # ------------------------------------------------------------------------
  # 2. Data Preparation
  # ------------------------------------------------------------------------

  X_concordant <- NULL
  y_concordant <- NULL
  treatment_concordant <- NULL

  X_diffs_discordant <- NULL
  y_diffs_discordant <- NULL
  treatment_diffs_discordant <- NULL

  # Internal variables for processing
  X <- data_mat
  X_model_matrix_col_names <- colnames(X)
  if (is.null(X_model_matrix_col_names)) {
    X_model_matrix_col_names <- paste0("X", 1:ncol(X))
  }

  if (length(y) <= ncol(X) + 2) {
    stop("Not enough rows. Must be at least the number of covariates plus 2")
  }


  matched_data <- process_matched_pairs_cpp(
    strata = strata,
    y = y,
    X = X,
    treatment = treatment
  )

  X_concordant <- matched_data$X_concordant
  y_concordant <- matched_data$y_concordant
  treatment_concordant <- matched_data$treatment_concordant
  strata_concordant <- matched_data$strata_concordant
  X_diffs_discordant <- matched_data$X_diffs_discordant
  y_diffs_discordant <- matched_data$y_diffs_discordant
  treatment_diffs_discordant <- matched_data$treatment_diffs_discordant

  # Ensure column names match for model formula consistency
  if (!is.null(X_concordant) && ncol(X_concordant) == length(X_model_matrix_col_names)) {
    colnames(X_concordant) <- X_model_matrix_col_names
  }

  # Early check: need enough discordant pairs to fit
  num_discordant <- length(y_diffs_discordant)
  if (is.null(X_diffs_discordant) || num_discordant < ncol(X) + 2) {
    stop(sprintf(
      "There are not enough discordant pairs (%d found). Need at least %d to fit the model.",
      num_discordant, ncol(X) + 2
    ))
  }

  # ------------------------------------------------------------------------
  # 3. Model Fitting
  # ------------------------------------------------------------------------


  # --- Concordant Pairs / Reservoir Model ---

  concordant_model <- NULL
  prior_fallbacks <- character(0)

  # Check for concordant pairs
  if (length(y_concordant) < ncol(X_concordant) + 5) {
    # Not enough data for concordant model
    warning("There are not enough concordant pairs or reservoir entries. The prior for the discordant pairs will be non-informative.")
    prior_fallbacks <- c(prior_fallbacks, "insufficient_concordant_pairs")
    b_con <- rep(0, ncol(X_concordant) + 1)
    Sigma_con <- diag(101, ncol(X_concordant) + 1)
  } else {
    if (concordant_method == "GLM") {
      concordant_model <- glm(y_concordant ~ treatment_concordant + X_concordant, family = "binomial")
    }
    if (concordant_method == "GEE") {
      gee_df <- data.frame(y_concordant, treatment_concordant, strata_concordant)
      concordant_model <- geepack::geeglm(
        y_concordant ~ treatment_concordant + X_concordant,
        id = strata_concordant,
        family = binomial(link = "logit"),
        corstr = "exchangeable",
        data = gee_df
      )
    }
    if (concordant_method == "GLMM") {
      concordant_model <- glmmTMB::glmmTMB(
        y_concordant ~ treatment_concordant + X_concordant + (1 | strata_concordant),
        family = binomial(),
        data   = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
      )
    }

    # Standardize prior extraction
    if (!is.null(concordant_model)) {
      # Get full coefficients including NAs for aliased terms
      if (concordant_method == "GLM") {
        full_b_con <- coef(concordant_model)
        full_Sigma_con <- vcov(concordant_model)
      } else if (concordant_method == "GEE") {
        full_b_con <- coef(concordant_model)
        # Use model-based (naive) covariance instead of sandwich estimator,
        # which can be degenerate (negative variances) for small cluster sizes
        full_Sigma_con <- concordant_model$geese$vbeta.naiv
        rownames(full_Sigma_con) <- names(full_b_con)
        colnames(full_Sigma_con) <- names(full_b_con)
      } else if (concordant_method == "GLMM") {
        full_b_con <- glmmTMB::fixef(concordant_model)$cond
        full_Sigma_con <- vcov(concordant_model)$cond
      }

      # Identify target names: treatment + X columns
      if (ncol(X_concordant) == 1) {
        target_names <- c("treatment_concordant", "X_concordant")
      } else {
        target_names <- c("treatment_concordant", paste0("X_concordant", X_model_matrix_col_names))
      }

      # Construct target vector
      K_stan <- ncol(X_diffs_discordant) + 1
      b_con <- numeric(K_stan)
      Sigma_con <- diag(100, K_stan)

      # Match coefficients
      for (i in seq_along(target_names)) {
        t_name <- target_names[i]
        if (t_name %in% names(full_b_con)) {
          val <- full_b_con[t_name]
          if (!is.na(val)) {
            b_con[i] <- val
          }
        }
      }

      # Match covariance
      valid_indices <- which(names(full_b_con) %in% target_names)
      valid_names <- names(full_b_con)[valid_indices]

      if (length(valid_names) > 0) {
        vcov_names <- rownames(full_Sigma_con)
        for (r_name in valid_names) {
          if (r_name %in% vcov_names) {
            target_idx_r <- match(r_name, target_names)
            for (c_name in valid_names) {
              if (c_name %in% vcov_names) {
                target_idx_c <- match(c_name, target_names)
                val <- full_Sigma_con[r_name, c_name]
                if (!is.na(val)) {
                  Sigma_con[target_idx_r, target_idx_c] <- val
                }
              }
            }
          }
        }
      }

      # Override treatment prior
      b_con[1] <- 0
      Sigma_con[1, ] <- 0
      Sigma_con[, 1] <- 0
      Sigma_con[1, 1] <- prior_variance_treatment
    }
  }

  # Sanitize priors to ensure no NA values are passed
  if (exists("b_con") && any(is.na(b_con))) {
    b_con[is.na(b_con)] <- 0
    prior_fallbacks <- c(prior_fallbacks, "na_in_prior_mean")
  }
  if (exists("Sigma_con")) {
    items_fixed <- FALSE
    if (any(is.na(Sigma_con))) {
      warning("Prior covariance contains NAs. Replacing with independent and diffuse prior.")
      Sigma_con <- diag(100, nrow(Sigma_con))
      prior_fallbacks <- c(prior_fallbacks, "na_in_prior_covariance")
      items_fixed <- TRUE
    }

    if (!items_fixed) {
      # Check positive definiteness
      is_pd <- tryCatch(
        {
          chol(Sigma_con)
          TRUE
        },
        error = function(e) FALSE
      )

      if (!is_pd) {
        if (concordant_method == "GEE") {
          warning("GEE prior covariance is not positive definite (common with small clusters). Falling back to GLM covariance for the prior.")
          prior_fallbacks <- c(prior_fallbacks, "glm_fallback")
          glm_fallback <- glm(y_concordant ~ treatment_concordant + X_concordant, family = "binomial")
          full_b_fb <- coef(glm_fallback)
          full_S_fb <- vcov(glm_fallback)
          # Re-extract into b_con / Sigma_con using the same matching logic
          b_con <- numeric(K_stan)
          Sigma_con <- diag(100, K_stan)
          for (i in seq_along(target_names)) {
            t_name <- target_names[i]
            if (t_name %in% names(full_b_fb)) {
              val <- full_b_fb[t_name]
              if (!is.na(val)) b_con[i] <- val
            }
          }
          vcov_names_fb <- rownames(full_S_fb)
          valid_names_fb <- names(full_b_fb)[names(full_b_fb) %in% target_names]
          for (r_name in valid_names_fb) {
            if (r_name %in% vcov_names_fb) {
              target_idx_r <- match(r_name, target_names)
              for (c_name in valid_names_fb) {
                if (c_name %in% vcov_names_fb) {
                  target_idx_c <- match(c_name, target_names)
                  val <- full_S_fb[r_name, c_name]
                  if (!is.na(val)) Sigma_con[target_idx_r, target_idx_c] <- val
                }
              }
            }
          }
          b_con[1] <- 0
          Sigma_con[1, ] <- 0
          Sigma_con[, 1] <- 0
          Sigma_con[1, 1] <- prior_variance_treatment
          Sigma_con <- (Sigma_con + t(Sigma_con)) / 2
        } else if (concordant_method == "GLMM") {
          warning("GLMM prior covariance is not positive definite (common with boundary variance estimates or convergence failures). Falling back to GLM covariance for the prior.")
          prior_fallbacks <- c(prior_fallbacks, "glm_fallback")
          glm_fallback <- glm(y_concordant ~ treatment_concordant + X_concordant, family = "binomial")
          full_b_fb <- coef(glm_fallback)
          full_S_fb <- vcov(glm_fallback)
          b_con <- numeric(K_stan)
          Sigma_con <- diag(100, K_stan)
          for (i in seq_along(target_names)) {
            t_name <- target_names[i]
            if (t_name %in% names(full_b_fb)) {
              val <- full_b_fb[t_name]
              if (!is.na(val)) b_con[i] <- val
            }
          }
          vcov_names_fb <- rownames(full_S_fb)
          valid_names_fb <- names(full_b_fb)[names(full_b_fb) %in% target_names]
          for (r_name in valid_names_fb) {
            if (r_name %in% vcov_names_fb) {
              target_idx_r <- match(r_name, target_names)
              for (c_name in valid_names_fb) {
                if (c_name %in% vcov_names_fb) {
                  target_idx_c <- match(c_name, target_names)
                  val <- full_S_fb[r_name, c_name]
                  if (!is.na(val)) Sigma_con[target_idx_r, target_idx_c] <- val
                }
              }
            }
          }
          b_con[1] <- 0
          Sigma_con[1, ] <- 0
          Sigma_con[, 1] <- 0
          Sigma_con[1, 1] <- prior_variance_treatment
          Sigma_con <- (Sigma_con + t(Sigma_con)) / 2
        } else {
          warning("Prior covariance is not positive definite. Replacing with independent and diffuse prior.")
          Sigma_con <- diag(100, nrow(Sigma_con))
          prior_fallbacks <- c(prior_fallbacks, "diffuse_prior")
        }
      }
    }
  }

  discordant_model <- NULL
  converged_discordant <- NULL

  # --- Discordant Pairs Model ---
  if (length(y_diffs_discordant) < ncol(X_diffs_discordant) + 2) {
    stop("There are not enough discordant pairs. The model will not be fit.")
  } else {
    Sigma_con <- (Sigma_con + t(Sigma_con)) / 2
    y_dis_0_1 <- ifelse(y_diffs_discordant == -1, 0, 1)
    wX_dis <- cbind(treatment_diffs_discordant, X_diffs_discordant)

    stan_model_name <- switch(prior_type,
      "Naive" = "mvn_logistic",
      "G prior" = "mvn_logistic_gprior",
      "PMP" = "mvn_logistic_PMP",
      "Hybrid" = "mvn_logistic_Hybrid",
      stop("Unknown prior type")
    )

    # Get pre-compiled model from rstantools
    stan_mod <- stanmodels[[stan_model_name]]

    if (prior_type %in% c("Naive", "G prior")) {
      data_list <- list(
        N = nrow(wX_dis),
        K = ncol(wX_dis),
        X = wX_dis,
        y = y_dis_0_1,
        mu = b_con,
        Sigma = Sigma_con
      )
    }
    if (prior_type %in% c("PMP", "Hybrid")) {
      # PMP and Hybrid share this setup
      proj_matrix <- X_diffs_discordant %*% solve(t(X_diffs_discordant) %*% X_diffs_discordant) %*% t(X_diffs_discordant)
      w_dis_ortho <- treatment_diffs_discordant - proj_matrix %*% treatment_diffs_discordant

      data_list <- list(
        N = nrow(X_diffs_discordant),
        P = ncol(X_diffs_discordant),
        y = y_dis_0_1,
        xw = as.vector(w_dis_ortho),
        X = X_diffs_discordant,
        mu_A = as.array(b_con[-1]),
        Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
      )
    }

    # Sample from the model
    discordant_model <- tryCatch(
      {
        rstan::sampling(stan_mod, data = data_list, refresh = stan_refresh, chains = chains, ...)
      },
      error = function(e) {
        warning(sprintf("Sampling failed: %s", e$message))
        NULL
      }
    )
  }

  # ------------------------------------------------------------------------
  # 4. Result Construction
  # ------------------------------------------------------------------------

  # Extract coefficients from the Stan model if available
  coefficients <- NULL
  var_cov <- NULL

  posterior_samples <- NULL

  if (!is.null(discordant_model)) {
    # Extract posterior samples
    sims <- rstan::extract(discordant_model)
    if (return_raw_stan_output) {
      # Raw samples array: iterations x chains x parameters
      posterior_samples <- rstan::extract(discordant_model, permuted = FALSE)
    }

    if (prior_type %in% c("Naive", "G prior")) {
      beta_post <- sims$beta
      col_names_final <- c(treatment_name, X_model_matrix_col_names)
    } else {
      # PMP / Hybrid
      beta_w_post <- as.vector(sims$beta_w)
      beta_nuis_post <- sims$beta_nuis

      beta_post <- cbind(beta_w_post, beta_nuis_post)
      col_names_final <- c(treatment_name, X_model_matrix_col_names)
    }

    # Calculate posterior means
    coefficients <- colMeans(beta_post)
    names(coefficients) <- col_names_final

    # Calculate posterior covariance
    var_cov <- cov(beta_post)
    # Assign names
    rownames(var_cov) <- col_names_final
    colnames(var_cov) <- col_names_final
  }

  # Construct result object
  res <- list(
    coefficients = coefficients,
    var = var_cov,
    model = discordant_model,
    posterior_samples = posterior_samples,
    concordant_model = concordant_model,
    matched_data = matched_data,
    prior_info = list(
      mu = if (exists("b_con")) b_con else NULL,
      Sigma = if (exists("Sigma_con")) Sigma_con else NULL,
      fallbacks = prior_fallbacks
    ),

    # Standard S3 components
    call = if (!is.null(call)) call else match.call(),
    terms = model_terms,
    xlevels = NULL,

    # Metadata for summary
    n = n,
    num_discordant = length(y_diffs_discordant),
    num_concordant = length(y_concordant) / 2,
    X_model_matrix_col_names = X_model_matrix_col_names,
    treatment_name = treatment_name
  )

  class(res) <- c("bclogit", "list")
  return(res)
}
