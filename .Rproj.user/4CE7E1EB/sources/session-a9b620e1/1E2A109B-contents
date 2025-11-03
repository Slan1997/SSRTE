#' Final analysis (inverse normal combination) after an interim with SSR
#'
#' Takes the output of \code{interim_analysis()} together with stage-2 data
#' and performs the final combined test using the inverse normal method:
#'
#' \deqn{ Z_{final} = \sqrt{n1/(n1+n2)} Z_1 + \sqrt{n2/(n1+n2)} Z_2 }
#'
#' where \eqn{n1/(n1+n2)} is the actual interim timing (information fraction),
#' \eqn{Z_1} is the interim Z from \code{interim_analysis()}, and \eqn{Z_2}
#' is the stage-2 Z computed from the stage-2 data.
#'
#' @param interim An object of class \code{"interim_analysis"} returned by
#'   \code{interim_analysis()}.
#' @param y_trt2 Stage-2 treatment data (matrix/data.frame: subjects x endpoints).
#' @param y_ct2 Stage-2 control data (matrix/data.frame: subjects x endpoints).
#' @param alloc_rate Treatment:control allocation (should match interim; default 1).
#'
#' @return An object of class \code{"final_analysis"}.
#'
#' @examples
#' \dontrun{
#' # assume `ia` is the object from interim_analysis(...)
#' # and we have simulated stage-2 data y_trt2, y_ct2
#' fa <- final_analysis(
#'   interim = ia,
#'   y_trt2 = y_trt2,
#'   y_ct2  = y_ct2,
#'
#' )
#' fa
#' summary(fa)
#' }
#'
#' @export
final_analysis <- function(
    interim,
    y_trt2,
    y_ct2,
    alloc_rate = 1
) {
  # ----- checks ---------------------------------------------------------------
  if (!inherits(interim, "interim_analysis")) {
    stop("`interim` must be the result of `interim_analysis()`.", call. = FALSE)
  }

  if (!is.matrix(y_trt2)) y_trt2 <- as.matrix(y_trt2)
  if (!is.matrix(y_ct2))  y_ct2  <- as.matrix(y_ct2)

  if (ncol(y_trt2) != ncol(y_ct2)) {
    stop("Stage-2 treatment and control must have the same number of endpoints.", call. = FALSE)
  }

  # ----- pull info from interim -----------------------------------------------
  z1 <- interim$test_statistic_Z
  if (is.null(z1)) {
    stop("Interim object does not contain `test_statistic_Z`.", call. = FALSE)
  }

  timing_actual <- interim$method$timing_actual
  if (is.null(timing_actual)) {
    stop("Interim object does not contain `method$timing_actual`.", call. = FALSE)
  }

  # C1 was used at interim; now get C2 from the **same** function
  bounds <- get_efficacy_boundaries(timing = timing_actual)
  C1 <- bounds[1]
  C2 <- bounds[2]

  # ----- stage-2 sample sizes -------------------------------------------------
  K      <- ncol(y_trt2)
  m_trt  <- nrow(y_trt2)
  m_ct   <- nrow(y_ct2)

  # ----- stage-2 covariance and correlation bits -----------------------------
  S2_pool <- (stats::cov(y_trt2) * (m_trt - 1) + stats::cov(y_ct2) * (m_ct - 1)) /
    (m_trt + m_ct - 2)

  sigma_hat2 <- sqrt(diag(S2_pool))
  if (any(sigma_hat2 == 0)) {
    stop("At least one endpoint has zero variance in stage 2.", call. = FALSE)
  }

  rho_ij2 <- (S2_pool - diag(S2_pool) * diag(1, K)) /
    (matrix(sigma_hat2, K, 1) %*% sigma_hat2)
  sum_rho2 <- sum(rho_ij2)

  A_obs2 <- sqrt((K + sum_rho2) / (K^2))

  cohen_d2      <- (colMeans(y_trt2) - colMeans(y_ct2)) / sigma_hat2
  mean_cohen_d2 <- mean(cohen_d2)

  # stage-2 t-like stat
  bar_T_star2 <- sqrt(m_ct / (1 / alloc_rate + 1)) * mean_cohen_d2 / A_obs2

  # df proposed by logan and Tamhane 2024
  new_df2 <- 0.5 * (m_ct * (1 + alloc_rate) - 2) * (1 + 1 / K^2)
  p_value_t2 <- stats::pt(bar_T_star2, df = new_df2, lower.tail = FALSE)
  z2 <- stats::qnorm(1 - p_value_t2)

  # ----- Combine via inverse normal ------------------------------------------
  z_final <- sqrt(timing_actual) * z1 + sqrt(1 - timing_actual) * z2
  p_final <- 1 - stats::pnorm(z_final)
  reject_final <- z_final > C2

  out <- list(
    decision = if (reject_final) "Reject H0 at final analysis." else "Fail to reject H0 at final analysis.",
    p_value_final = p_final,
    test_statistic_final_Z = z_final,
    boundary_final = C2,
    interim_boundary = C1,
    stage1 = list(
      Z = z1,
      timing = timing_actual
    ),
    stage2 = list(
      Z = z2,
      p_value_t = p_value_t2,
      m_ct = m_ct,
      m_trt = m_trt,
      A_obs2 = A_obs2,
      mean_cohen_d2 = mean_cohen_d2,
      sum_rho2 = sum_rho2
    ),
    method = list(
      combination = "inverse-normal",
      timing_actual = timing_actual,
      alloc_rate = alloc_rate
    ),
    sample_sizes = list(
      stage1 = interim$sample_sizes$stage1,
      stage2 = c(control = m_ct, treatment = m_trt),
      planned_final = interim$sample_sizes$final
    )
  )

  class(out) <- c("final_analysis", "list")
  out
}

#' @export
print.final_analysis <- function(x, digits = 3, ...) {
  stopifnot(inherits(x, "final_analysis"))
  fmt <- function(v) if (is.na(v)) "NA" else formatC(v, digits = digits, format = "f")

  cli::cli_rule("Final analysis")
  cli::cli_text("{.strong Decision}: {x$decision}")

  cli::cli_rule("Combined result")
  cli::cli_text("Final Z: {fmt(x$test_statistic_final_Z)}")
  cli::cli_text("Final p-value: {fmt(x$p_value_final)}")
  cli::cli_text("Final boundary (C2): {fmt(x$boundary_final)}")

  cli::cli_rule("Interim component")
  cli::cli_text("Interim Z (Z1): {fmt(x$stage1$Z)}")
  cli::cli_text("Interim boundary (C1): {fmt(x$interim_boundary)}")
  cli::cli_text("Timing: {fmt(x$stage1$timing)}")

  cli::cli_rule("Stage 2 component")
  cli::cli_text("Stage 2 Z (Z2): {fmt(x$stage2$Z)}")
  cli::cli_text("Stage 2 p-value: {fmt(x$stage2$p_value_t)}")
  cli::cli_text("Stage 2 n: control = {x$stage2$m_ct}, treatment = {x$stage2$m_trt}")

  cli::cli_rule("Stage 2 estimates")
  cli::cli_dl(c(
    "Mean Cohen's d (stage 2)" = fmt(x$stage2$mean_cohen_d2),
    "A_obs2 (SE of bar_t_2)" = fmt(x$stage2$A_obs2),
    "Sum rho (stage 2)" = fmt(x$stage2$sum_rho2)
  ))

  invisible(x)
}

#' @export
summary.final_analysis <- function(object, ...) {
  stopifnot(inherits(object, "final_analysis"))
  tibble::tibble(
    decision      = object$decision,
    Z_final       = object$test_statistic_final_Z,
    p_final       = object$p_value_final,
    C1            = object$interim_boundary,
    C2            = object$boundary_final,
    Z_stage1      = object$stage1$Z,
    timing        = object$stage1$timing,
    Z_stage2      = object$stage2$Z,
    p_stage2    = object$stage2$p_value_t,
    m_ct_stage2   = object$stage2$m_ct,
    m_trt_stage2  = object$stage2$m_trt
  )
}
