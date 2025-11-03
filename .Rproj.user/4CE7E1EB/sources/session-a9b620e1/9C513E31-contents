#' Root function for CP-based SSR
#'
#' This is a small helper meant to be used **inside** \code{uniroot()} when
#' doing sample size re-estimation based on a *target conditional power*.
#'
#' Given a candidate control sample size \code{n}, it:
#'
#' 1. computes the conditional power via \code{get_CP(...)} using
#'    Eq. (2.10), and
#' 2. returns \code{CP(n) - target_CP}.
#'
#' Then \code{uniroot()} can be used to find the \code{n} such that
#' \code{CP(n) - target_CP = 0}.
#'
#' @param n Numeric, **candidate** final control group sample size (this is the
#'   variable over which \code{uniroot()} searches).
#' @param timing Numeric in (0, 1), interim information fraction.
#' @param mean_cohen_d Numeric, observed mean Cohen’s d at interim.
#' @param K Integer, number of endpoints.
#' @param sum_rho Numeric, sum of all off-diagonal endpoint correlations at
#'   interim, i.e. \eqn{\sum_{p \ne q} \hat\rho_{pq,1}}.
#' @param target_CP Numeric in (0, 1), the conditional power you want to hit
#' @param alloc_rate Numeric > 0, treatment:control allocation ratio; default \code{1}.
#'
#' @return A single numeric \code{get_CP(...) - target_CP}. Suitable for
#'   \code{uniroot()}.
#'
#' @examples
#' \dontrun{
#' root_res <- stats::uniroot(
#'   target_CP_get_N_ct,
#'   lower = 10,
#'   upper = 1000,
#'   timing = 0.5,
#'   mean_cohen_d = 0.3,,
#'   K = 6,
#'   sum_rho = 9,
#'   target_CP = 0.8,
#'   alloc_rate = 1
#'   )
#' n_hat <- ceiling(root_res$root)
#' }
#'
#' @export
target_CP_get_N_ct <- function(timing,
                               n,
                               mean_cohen_d,
                               K,
                               sum_rho,
                               target_CP,
                               alloc_rate = 1) {

  CP <- get_CP(
    timing       = timing,
    n            = n,
    mean_cohen_d = mean_cohen_d,
    K            = K,
    sum_rho      = sum_rho,
    alloc_rate   = alloc_rate
  )

  CP - target_CP
}
