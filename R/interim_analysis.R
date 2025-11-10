#' Interim analysis with sample size re-estimation
#'
#' Performs an interim global test (exact OLS or permutation). If the interim
#' statistic crosses the efficacy boundary \code{C1}, the trial stops. Otherwise,
#' it can do SSR in the promising zone using either \code{SSR-Power} or \code{SSR-CP}.
#'
#' Stage-1 sample sizes are taken directly from \code{nrow(y_trt1)} and
#' \code{nrow(y_ct1)}.
#'
#' @param y_trt1 Stage-1 treatment data (matrix/data.frame: subjects x endpoints).
#' @param y_ct1 Stage-1 control data (matrix/data.frame: subjects x endpoints).
#' @param N_trt Planned total treatment sample size before SSR.
#' @param N_ct Planned total control sample size before SSR.
#' @param N_ct_max Maximum *re-estimated* control sample size allowed.
#' @param timing_actual Interim timing (information fraction, within (0,1))
#' @param alpha0 One-sided alpha for the power-based SSR formula.
#' @param beta0 Type II error rate (\eqn{\beta}).
#' @param alloc_rate Treatment:control allocation.
#' @param n_endpts Number of endpoints.
#' @param SSR_type \code{"SSR-Power"}, \code{"SSR-CP"}, or \code{"No SSR"}.
#' @param global_test_type \code{"exact OLS"} or \code{"Permutation"}. If
#'   \code{"Permutation"} is used, a random seed is needed to ensure reproducibility.
#' @param nPM Number of permutations for the permutation test.
#' @param promising_LL Lower bound of the promising zone. Default: 0.2. (Upper bound: 1-beta0)
#'
#' @details
#' This function calls the following **in-package** helpers:
#'
#' - \code{get_CP(timing, n, mean_cohen_d, K, sum_rho)}
#' - \code{target_CP_get_N_ct(N_ct, timing, mean_cohen_d, K, sum_rho, target_CP)}
#' - \code{get_1_PM_mean_z_score(y_mat, trt)}
#'
#' @return An object of class \code{"interim_analysis"}.
#'
#' @examples
#' rho <- matrix(c(
#'  1.0, 0.1, 0.3, 0.7, 0.1, 0.3,
#'  0.1, 1.0, 0.7, 0.1, 0.3, 0.7,
#'  0.3, 0.7, 1.0, 0.1, 0.3, 0.7,
#'  0.7, 0.1, 0.1, 1.0, 0.1, 0.3,
#'  0.1, 0.3, 0.3, 0.1, 1.0, 0.7,
#'  0.3, 0.7, 0.7, 0.3, 0.7, 1.0
#'  ), nrow = 6, byrow = TRUE)
#'
#' sim <- simulate_example_one_stage_data(
#'   n_trt = 40, n_ct = 40,
#'   n_endpts = 6,
#'   exp_mean_cohen_d = 0.3,
#'   mu_ct = c(4, 6, 4, 5, 7, 6),
#'   sd = c(0.5, 0.5, 1, 1, 2, 2),
#'   rho = rho,
#'   seed = 15
#' )
#' y_ct1  <- sim$y_ct
#' y_trt1 <- sim$y_trt
#'
#' # set.seed(1) # seed needed for "Permutation"
#' ia <- interim_analysis(
#'   y_trt1 = y_trt1,
#'   y_ct1  = y_ct1,
#'   N_trt = 80,
#'   N_ct = 80,
#'   N_ct_max = 160,
#'   timing_actual = 0.5,
#'   alpha0 = 0.025,
#'   beta0 = 0.2,
#'   alloc_rate = 1,
#'   n_endpts = 6,
#'   SSR_type = "SSR-Power",
#'   global_test_type = "exact OLS"
#' )
#' ia
#' summary(ia)
#'
#' @export
interim_analysis <- function(
    y_trt1, y_ct1,
    N_trt, N_ct,
    N_ct_max,
    timing_actual,
    alpha0,
    beta0,
    alloc_rate,
    n_endpts,
    SSR_type,
    global_test_type,
    nPM = 1e2,
    promising_LL = 0.2
) {
  # --- derive stage-1 sizes ---------------------------------------------------
  if (!is.matrix(y_trt1)) y_trt1 <- as.matrix(y_trt1)
  if (!is.matrix(y_ct1))  y_ct1  <- as.matrix(y_ct1)

  n_trt <- nrow(y_trt1)
  n_ct  <- nrow(y_ct1)

  if (ncol(y_trt1) != ncol(y_ct1)) {
    stop("Treatment and control must have the same number of endpoints.", call. = FALSE)
  }

  # --- pooled covariance & correlations --------------------------------------
  S1_pool <- (stats::cov(y_trt1) * (n_trt - 1) + stats::cov(y_ct1) * (n_ct - 1)) /
    (n_trt + n_ct - 2)

  sigma_hat <- sqrt(diag(S1_pool))
  if (any(sigma_hat == 0)) {
    stop("At least one endpoint has zero variance in stage 1.", call. = FALSE)
  }

  rho_ij <- (S1_pool - diag(S1_pool) * diag(1, n_endpts)) /
    (matrix(sigma_hat, n_endpts, 1) %*% sigma_hat)

  sum_rho <- sum(rho_ij)
  A_obs   <- sqrt((n_endpts + sum_rho) / (n_endpts^2))

  # --- effect sizes -----------------------------------------------------------
  cohen_d      <- (colMeans(y_trt1) - colMeans(y_ct1)) / sigma_hat
  mean_cohen_d <- mean(cohen_d)

  # --- global test ------------------------------------------------------------
  if (global_test_type=="exact OLS") {
    bar_t_1     <- (1 / sqrt(1 / n_trt + 1 / n_ct)) * mean_cohen_d
    bar_T_star1 <- bar_t_1 / A_obs
    # df proposed by logan and Tamhane 2024
    new_df      <- 0.5 * (n_ct + n_trt - 2) * (1 + 1 / n_endpts^2)
    p_value_t   <- stats::pt(bar_T_star1, df = new_df, lower.tail = FALSE)
    bar_Z_star1 <- stats::qnorm(1 - p_value_t)
    test_stat_interim <- bar_Z_star1 # test stats is Z-score
    p_value           <- p_value_t
  }else{
    # permutation
    y_mat1   <- rbind(y_trt1, y_ct1)
    trt_ind1 <- rep(c(1, 0), c(n_trt, n_ct))

    mean_z_score <- mean(apply(y_mat1, 2, function(y) {
      -as.numeric(stats::t.test(y ~ trt_ind1)$statistic)
    }))

    perm_mean <- rep(NA, nPM)
    for (pm in seq_len(nPM)) {
      perm_mean[pm] <- get_1_PM_mean_z_score(y_mat = y_mat1, trt = trt_ind1)
    }
    p_value1          <- mean(perm_mean > mean_z_score)
    test_stat_interim <- stats::qnorm(1 - p_value1)
    p_value           <- p_value1
  }

  # --- early stopping ---------------------------------------------------------
  C1 = get_efficacy_boundaries(timing = timing_actual)[1]
  if (test_stat_interim > C1) {
    out <- list(
      decision = "Reject H0 at interim; stop the trial.",
      p_value = p_value,
      test_statistic_Z = test_stat_interim,
      boundary = C1,
      method = list(
        global_test_type = global_test_type,
        SSR_type = SSR_type,
        n_permutations = if (identical(global_test_type, "Permutation")) nPM else NA_integer_
      ),
      estimates = list(
        mean_cohen_d = mean_cohen_d,
        A_obs = A_obs,
        sum_rho = sum_rho
      ),
      ssr = list(
        estimated_CP = NA_real_,
        promising_zone = FALSE,
        n_hat_ct = NA_integer_,
        N_ct_max = N_ct_max
      ),
      sample_sizes = list(
        stage1 = c(control = n_ct, treatment = n_trt),
        stage2 = c(control = 0L, treatment = 0L),
        final  = c(control = n_ct, treatment = n_trt)
      )
    )
    class(out) <- c("interim_analysis", "list")
    return(out)
  }

  # --- not stopped: consider SSR ----------------------------------------------
  mean_cohen_d_adj <- max(mean_cohen_d, 0.001) # # when mean cohen_d <= 0, adjust it to be 0.01 to avoid denominator = 0 in SSR formula
  est_CP           <- NA_real_
  promising_zone   <- FALSE
  n_hat            <- NA_integer_

  if (SSR_type!="No SSR") {
    # conditional power at planned n
    est_CP <- get_CP(
      timing = timing_actual,
      n = N_ct,  # planned control size before SSR
      mean_cohen_d = mean_cohen_d_adj,
      K = n_endpts,
      sum_rho = sum_rho
    )
    if (est_CP > promising_LL && est_CP < (1 - beta0)) {
      promising_zone <- TRUE

      if (SSR_type == "SSR-Power") {
        n_hat <- ceiling(
          A_obs^2 *
            ((stats::qnorm(1 - alpha0) + stats::qnorm(1 - beta0)) *
               sqrt(1 / alloc_rate + 1) / mean_cohen_d_adj)^2
        )
      } else {
        # SSR-CP
        if ((N_trt + N_ct) > 500) {
          LL <- 100
          UL <- 5000
        } else {
          LL <- 10
          UL <- 1000
        }
        root_res <- stats::uniroot(
          target_CP_get_N_ct,
          lower = LL,
          upper = UL,
          timing = timing_actual,
          mean_cohen_d = mean_cohen_d_adj,
          K = n_endpts,
          sum_rho = sum_rho,
          target_CP = 1 - beta0
        )
        n_hat <- ceiling(root_res$root)
      }

      # Stage 2 reestimated sample size, cap within [N_ct, N_ct_max]
      M_ct  <- min(max(n_hat, N_ct), N_ct_max)
      M_trt <- ceiling(alloc_rate * M_ct)
    } else {
      # favorable/unfavorable → keep original planned
      M_ct  <- N_ct
      M_trt <- N_trt
    }
  } else {
    # No SSR
    M_ct  <- N_ct
    M_trt <- N_trt
  }

  m_ct  <- M_ct - n_ct
  m_trt <- M_trt - n_trt

  out <- list(
    decision = "Fail to reject H0 at interim; continue the trial.",
    p_value = p_value,
    test_statistic_Z = test_stat_interim,
    boundary = C1,
    method = list(
      global_test_type = global_test_type,
      SSR_type = SSR_type,
      timing_actual = timing_actual,
      n_permutations = if (identical(global_test_type, "Permutation")) nPM else NA_integer_
    ),
    estimates = list(
      mean_cohen_d = mean_cohen_d,
      A_obs = A_obs,
      sum_rho = sum_rho
    ),
    ssr = list(
      estimated_CP = est_CP,
      promising_zone = promising_zone,
      n_hat_ct = n_hat,
      N_ct_max = N_ct_max
    ),
    sample_sizes = list(
      stage1 = c(control = n_ct, treatment = n_trt),
      stage2 = c(control = m_ct, treatment = m_trt),
      final  = c(control = M_ct, treatment = M_trt)
    )
  )
  class(out) <- c("interim_analysis", "list")
  out
}


#' @export
print.interim_analysis <- function(x, digits = 3, ...) {
  stopifnot(inherits(x, "interim_analysis"))
  fmt <- function(v) if (is.na(v)) "NA" else formatC(v, digits = digits, format = "f")

  cli::cli_rule("Interim analysis")
  cli::cli_text("{.strong Decision}: {x$decision}")

  cli::cli_text("Global test: {x$method$global_test_type}")
  cli::cli_text("SSR type: {x$method$SSR_type}")
  if (!is.na(x$method$n_permutations)) {
    cli::cli_text("Permutations: {x$method$n_permutations}")
  }

  cli::cli_rule("Results")
  cli::cli_text("Interim boundary (C1): {fmt(x$boundary)}")
  cli::cli_text("Test statistic (Z): {fmt(x$test_statistic_Z)}")
  cli::cli_text("p-value: {fmt(x$p_value)}")

  cli::cli_rule("Estimates")
  cli::cli_dl(c(
    "Mean Cohen's d" = fmt(x$estimates$mean_cohen_d),
    "A_obs (SE of bar_t_1)" = fmt(x$estimates$A_obs),
    "Sum of pairwise correlations" = fmt(x$estimates$sum_rho)
  ))

  cli::cli_rule("SSR")
  cli::cli_text("Estimated conditional power: {fmt(x$ssr$estimated_CP)}")
  cli::cli_text("Promising zone? {if (isTRUE(x$ssr$promising_zone)) 'yes' else 'no'}")
  if (!is.na(x$ssr$n_hat_ct)) {
    cli::cli_text("Estimated control N (before capping): {x$ssr$n_hat_ct}")
  }
  cli::cli_text("Max allowed control N: {x$ssr$N_ct_max}")

  cli::cli_rule("Re-estimated Sample sizes")
  s1 <- x$sample_sizes$stage1
  s2 <- x$sample_sizes$stage2
  sf <- x$sample_sizes$final
  cli::cli_text("{.strong Stage 1} | Control: {s1['control']}, Treatment: {s1['treatment']}")
  cli::cli_text("{.strong Stage 2} | Control: {s2['control']}, Treatment: {s2['treatment']}")
  cli::cli_text("{.strong Final}   | Control: {sf['control']}, Treatment: {sf['treatment']}")

  invisible(x)
}

#' @export
summary.interim_analysis <- function(object, ...) {
  stopifnot(inherits(object, "interim_analysis"))
  tibble::tibble(
    decision        = object$decision,
    global_test     = object$method$global_test_type,
    SSR_type        = object$method$SSR_type,
    test_statistic_Z  = object$test_statistic_Z,
    p_value         = object$p_value,
    C1              = object$boundary,
    mean_cohen_d    = object$estimates$mean_cohen_d,
    A_obs           = object$estimates$A_obs,
    sum_rho         = object$estimates$sum_rho,
    est_CP          = object$ssr$estimated_CP,
    promising_zone  = object$ssr$promising_zone,
    n_hat_ct        = object$ssr$n_hat_ct,
    n_ct_stage1     = object$sample_sizes$stage1["control"],
    n_trt_stage1    = object$sample_sizes$stage1["treatment"],
    n_ct_stage2     = object$sample_sizes$stage2["control"],
    n_trt_stage2    = object$sample_sizes$stage2["treatment"],
    n_ct_final      = object$sample_sizes$final["control"],
    n_trt_final     = object$sample_sizes$final["treatment"]
  )
}


