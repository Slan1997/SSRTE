#' One permutation mean z-score
#'
#' This helper performs **one** permutation step for the permutation test
#' described in Li et al. (2020).
#'
#' 1. randomly permutes the treatment indicator,
#' 2. for **each endpoint** runs a two-sample \code{t.test()},
#' 3. takes the \emph{negative} of the test statistic (to match the direction in
#'    the interim code),
#' 4. returns the **mean** of those endpointwise test statistics.
#'
#'
#' @param y_mat Numeric matrix (or coercible to matrix) of outcomes with
#'   \code{n} rows (subjects) and \code{K} columns (endpoints).
#' @param trt Integer or numeric vector of length \code{n} giving the original
#'   treatment indicator (e.g. 1 = treatment, 0 = control).
#'
#' @return A single numeric value: the permuted mean z-score.
#'
#' @examples
#' \dontrun{
#' y_mat1   <- rbind(y_trt1, y_ct1)
#' trt_ind1 <- rep(c(1, 0), c(n_trt1, n_ct1))
#' mean_z_score <- mean(apply(y_mat1, 2, function(y) {
#'   -as.numeric(stats::t.test(y ~ trt_ind1)$statistic)}))
#' perm_mean <- rep(NA, nPM)
#' for (pm in seq_len(nPM)) {
#'   perm_mean[pm] <- get_1_PM_mean_z_score(y_mat = y_mat1, trt = trt_ind1)
#' }
#' p_value1          <- mean(perm_mean > mean_z_score)
#' test_stat_interim <- stats::qnorm(1 - p_value1)
#' }
#'
#' @export
get_1_PM_mean_z_score <- function(y_mat, trt) {
  if (!is.matrix(y_mat)) {
    y_mat <- as.matrix(y_mat)
  }
  n <- nrow(y_mat)
  if (length(trt) != n) {
    stop("`trt` must have the same length as the number of rows in `y_mat`.", call. = FALSE)
  }

  # permute treatment labels
  pm_trt <- sample(trt)

  # endpointwise t-tests, then average the test statistics
  mean(apply(
    y_mat,
    2,
    function(y) -as.numeric(stats::t.test(y ~ pm_trt)$statistic)
  ))
}
