#' Simulate one-stage multivariate continuous endpoint data
#'
#' This function generates simulated treatment and control data
#' with the **same covariance structure** but different means based on
#' a target (expected) Cohen’s d value.
#'
#' @param n_trt Integer. Number of subjects in the treatment group (one-stage).
#' @param n_ct Integer. Number of subjects in the control group (one-stage).
#' @param n_endpts Integer. Number of endpoints.
#' @param exp_mean_cohen_d Numeric scalar. Expected mean Cohen’s d across all endpoints.
#'   The mean shift is computed as \code{exp_mean_cohen_d * sd}.
#' @param mu_ct Numeric vector of length \code{n_endpts}. Mean vector for the control group.
#' @param sd Numeric vector of length \code{n_endpts}. Standard deviations for each endpoint.
#' @param rho Numeric \code{n_endpts x n_endpts} correlation matrix.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list with two matrices:
#' \describe{
#'   \item{y_trt}{Matrix of treatment group data (\code{n_trt x n_endpts}).}
#'   \item{y_ct}{Matrix of control group data (\code{n_ct x n_endpts}).}
#' }
#'
#' @details
#' The covariance matrix is constructed as:
#' \deqn{\Sigma = \text{diag}(sd) \, \rho \, \text{diag}(sd)}
#'
#' The treatment mean vector is computed as:
#' \deqn{\mu_{\text{trt}} = \mu_{\text{ct}} + (\text{exp\_mean\_cohen\_d} \times sd)}
#'
#' @examples
#'   rho <- matrix(c(
#'     1.0, 0.1, 0.3, 0.7, 0.1, 0.3,
#'     0.1, 1.0, 0.7, 0.1, 0.3, 0.7,
#'     0.3, 0.7, 1.0, 0.1, 0.3, 0.7,
#'     0.7, 0.1, 0.1, 1.0, 0.1, 0.3,
#'     0.1, 0.3, 0.3, 0.1, 1.0, 0.7,
#'     0.3, 0.7, 0.7, 0.3, 0.7, 1.0
#'   ), nrow = 6, byrow = TRUE)
#'
#'   sim <- simulate_example_one_stage_data(
#'     n_trt = 40, n_ct = 40,
#'     n_endpts = 6,
#'     exp_mean_cohen_d = 0.3,
#'     mu_ct = c(4, 6, 4, 5, 7, 6),
#'     sd = c(0.5, 0.5, 1, 1, 2, 2),
#'     rho = rho,
#'     seed = 15
#'   )
#'   str(sim)
#'
#' @importFrom mvtnorm rmvnorm
#' @export
simulate_example_one_stage_data <- function(n_trt, n_ct,
                                            n_endpts,
                                            exp_mean_cohen_d,
                                            mu_ct, sd,
                                            rho,
                                            seed) {
  # --- argument checks ------------------------------------------------------
  if (!is.numeric(n_trt) || length(n_trt) != 1L || n_trt <= 0) {
    stop("`n_trt` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n_ct) || length(n_ct) != 1L || n_ct <= 0) {
    stop("`n_ct` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n_endpts) || length(n_endpts) != 1L || n_endpts <= 0) {
    stop("`n_endpts` must be a positive integer.", call. = FALSE)
  }
  if (length(mu_ct) != n_endpts)
    stop("`mu_ct` must have length `n_endpts`.", call. = FALSE)
  if (length(sd) != n_endpts)
    stop("`sd` must have length `n_endpts`.", call. = FALSE)
  if (!is.matrix(rho) || any(dim(rho) != n_endpts))
    stop("`rho` must be an n_endpts x n_endpts matrix.", call. = FALSE)
  if (any(sd <= 0))
    stop("All entries in `sd` must be positive.", call. = FALSE)

  # --- compute mean shift and covariance ------------------------------------
  exp_cohen_d <- rep(exp_mean_cohen_d, n_endpts)
  exp_delta   <- exp_cohen_d * sd
  mu_trt      <- mu_ct + exp_delta
  cov_mat     <- outer(sd, sd) * rho

  # --- simulate -------------------------------------------------------------
  set.seed(seed)
  y_trt_mat <- mvtnorm::rmvnorm(n_trt, mean = mu_trt, sigma = cov_mat)
  y_ct_mat  <- mvtnorm::rmvnorm(n_ct,  mean = mu_ct,  sigma = cov_mat)

  list(y_trt = y_trt_mat, y_ct = y_ct_mat)
}
