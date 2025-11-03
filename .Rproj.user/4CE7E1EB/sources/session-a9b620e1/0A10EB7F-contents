#' Conditional power at interim (Eq. 2.10)
#'
#' @param timing Numeric in (0, 1). Interim information fraction.
#' @param n Numeric, planned final **control group** sample size (before SSR).
#' @param mean_cohen_d Numeric, observed mean Cohen's d at interim
#'   (i.e. \eqn{\bar d_1}).
#' @param K Integer, number of endpoints.
#' @param sum_rho Numeric, sum of all off-diagonal interim correlations
#'   \eqn{\sum_{p \ne q} \hat\rho_{pq,1}}.
#' @param alloc_rate Numeric > 0, treatment:control allocation ratio. Default \code{1}.
#'
#' @return Numeric scalar: the conditional power.
#'
#' @details
#' This function calls \code{get_efficacy_boundaries(timing)} and takes the
#' final-stage boundary \eqn{z_{\alpha_2}}
#'
#' @examples
#' # 50% look, planned 80 controls total, 1:1 allocation, 6 endpoints
#' get_CP(
#'   timing = 0.5,
#'   n = 80,
#'   mean_cohen_d = 0.3,
#'   K = 6,
#'   sum_rho = 9,
#'   alloc_rate = 1
#' )
#'
#' @export
get_CP <- function(timing,
                   n,
                   mean_cohen_d,
                   K,
                   sum_rho,
                   alloc_rate = 1) {

  if (timing <= 0 || timing >= 1) {
    stop("`timing` must be in (0, 1).", call. = FALSE)
  }
  if (n <= 0) {
    stop("`n` must be > 0 (planned final control sample size).", call. = FALSE)
  }
  if (alloc_rate <= 0) {
    stop("`alloc_rate` must be > 0 (treatment:control ratio).", call. = FALSE)
  }

  # stage-1 and stage-2 CONTROL sizes
  n1 <- round(n * timing)
  n2 <- n - n1
  if (n1 <= 0 || n2 <= 0) {
    stop("`timing * n` must give positive stage-1 and stage-2 sizes.", call. = FALSE)
  }

  # final-stage z-boundary (z_{alpha2})
  C2 <- get_efficacy_boundaries(timing = timing)[2]

  # denominator: A_obs = sqrt( (1 / K^2) * (K + sum_rho) )
  A_obs <- sqrt((K + sum_rho) / (K^2))

  # stage-1 average t-statistic with allocation r
  bar_t_1 <- mean_cohen_d / sqrt(1 / (n1 * alloc_rate) + 1 / n1)

  # Eq. (2.10) numerator pieces
  term1 <- (sqrt(n1) * bar_t_1 - sqrt(n) * C2) / sqrt(n2)
  term2 <- (sqrt(n2) * bar_t_1) / sqrt(n1)

  CP <- stats::pnorm((term1 + term2) / A_obs)
  CP
}
