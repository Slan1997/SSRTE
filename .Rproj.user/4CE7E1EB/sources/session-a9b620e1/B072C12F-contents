#' Get the initial planned total sample size for a multivariate trial
#'
#' This computes the **planned** total sample size under a
#' multivariate setting with \eqn{K} endpoints, assuming a common allocation
#' ratio across endpoints.
#'
#' The formula implements Eq. (2.11) in Shi et al. (2025):
#'
#' \deqn{
#' N_{tot,0}
#'   = \frac{K + 2\sum_{p<q} \rho_{pq}}{K^2}
#'     (1/r + 1)(1 + r)
#'     \left(
#'       \frac{z_{\alpha} + z_{\beta}}{\bar{\theta}*}
#'     \right)^2,
#' }
#'
#' where \eqn{K =} \code{n_endpts}, \eqn{r} is the allocation ratio
#' (treatment/control), and \eqn{\bar{\theta}*} is the expected mean Cohen’s d across endpoints
#'
#' @param beta0 Type II error rate (\eqn{\beta}). Power is \code{1 - beta0}.
#' @param r Allocation ratio = treatment / control. For 1:1 allocation, use 1.
#' @param theta_k Expected treatment effect size, i.e., expected mean Cohen’s d across all endpoints (\eqn{\bar{\theta}*}).
#' @param n_endpts Integer. Number of endpoints.
#' @param alpha0 One-sided type I error rate. Default \code{0.025}.
#' @param rho Either a single correlation value (assumed common across all
#'   endpoint pairs) or a \code{n_endpts x n_endpts} correlation matrix.
#' @param timing0 Planned interim timing (information fraction) in \code{(0, 1]},
#'   e.g. \code{0.5} for an interim at half of the total sample.
#'
#' @return A list with elements:
#' \describe{
#'   \item{N_tot0}{Planned total sample size under the formula (a crude result, can be non-integer).}
#'   \item{N_ct}{Planned **control** total sample size (ceiling).}
#'   \item{N_trt}{Planned **treatment** total sample size (ceiling).}
#'   \item{n_ct}{Control sample size at interim (rounded).}
#'   \item{n_trt}{Treatment sample size at interim (rounded).}
#'   \item{N_interim_actual}{Actual total sample size at interim = \code{n_ct + n_trt}.}
#'   \item{N_total_actual}{Final total sample size = \code{N_ct + N_trt}.}
#'   \item{timing_actual}{Actual interim fraction = \code{N_interim_actual / N_total_actual}.}
#' }
#'
#' @examples
#' # 1:1 allocation, 6 endpoints, common correlation 0.3, interim at 50%
#' get_initial_tot_ss(
#'   beta0  = 0.2,
#'   r      = 1,
#'   theta_k = 0.386,
#'   n_endpts = 6,
#'   alpha0 = 0.025,
#'   rho    = 0.3,
#'   timing0 = 0.5
#' )
#'
#' # With a correlation matrix
#' rho <- matrix(c(
#'   1.0, 0.1, 0.3, 0.7, 0.1, 0.3,
#'   0.1, 1.0, 0.7, 0.1, 0.3, 0.7,
#'   0.3, 0.7, 1.0, 0.1, 0.3, 0.7,
#'   0.7, 0.1, 0.1, 1.0, 0.1, 0.3,
#'   0.1, 0.3, 0.3, 0.1, 1.0, 0.7,
#'   0.3, 0.7, 0.7, 0.3, 0.7, 1.0
#' ), nrow = 6, byrow = TRUE)
#'
#' get_initial_tot_ss(
#'   beta0  = 0.2,
#'   r      = 1,
#'   theta_k = 0.386,
#'   n_endpts = 6,
#'   alpha0 = 0.025,
#'   rho    = rho,
#'   timing0 = 0.5
#' )
#'
#' @export
get_initial_tot_ss <- function(beta0,
                               r,
                               theta_k,
                               n_endpts,
                               alpha0 = 0.025,
                               rho,
                               timing0) {
  # --- correlation summary ---------------------------------------------------
  if (length(rho) == 1L) {
    # common correlation across all pairs
    rho_sum <- n_endpts * (n_endpts - 1) * rho
  } else {
    # rho is a matrix: sum of off-diagonals
    if (!is.matrix(rho) || nrow(rho) != n_endpts || ncol(rho) != n_endpts) {
      stop("`rho` must be either a single number or an n_endpts x n_endpts matrix.",
           call. = FALSE)
    }
    rho_sum <- sum(rho) - n_endpts  # subtract diagonal 1's
  }

  # --- total planned sample size (may be non-integer) ------------------------
  z_alpha <- stats::qnorm(1 - alpha0)
  z_beta  <- stats::qnorm(1 - beta0)

  N_tot0 <- ((n_endpts + rho_sum) / n_endpts^2) *
    (1 / r + 1) * (1 + r) *
    ((z_alpha + z_beta) / theta_k)^2

  # --- split into arms -------------------------------------------------------
  # control total
  N_ct <- ceiling(N_tot0 / (1 + r))
  # interim control
  n_ct <- round(N_ct * timing0)
  # treatment totals (ceiling to keep allocation rate)
  N_trt <- ceiling(N_ct * r)
  n_trt <- ceiling(n_ct * r)

  N_interim_actual <- n_ct + n_trt
  N_total_actual   <- N_ct + N_trt
  timing_actual    <- N_interim_actual / N_total_actual

  list(
    N_tot0           = N_tot0,
    N_ct             = N_ct,
    N_trt            = N_trt,
    n_ct             = n_ct,
    n_trt            = n_trt,
    N_interim_actual = N_interim_actual,
    N_total_actual   = N_total_actual,
    timing_actual    = timing_actual
  )
}
