#' Get one-sided Lan-DeMets O'Brien-Fleming efficacy boundaries (2-look design)
#'
#' Compute one-sided group-sequential **efficacy** boundaries for a
#' 2-stage (1 interim + final) trial using the Lan-DeMets
#' O'Brien-Fleming-type alpha-spending function from \pkg{gsDesign}.
#'
#' This is a small convenience wrapper around \code{gsDesign()} for the
#' very common case where one has **one interim at time t** and then
#' the final at time 1.
#'
#' @param timing Numeric scalar in \code{(0, 1)} giving the **interim**
#'   information fraction (e.g. \code{0.5} for 50 percent information at interim).
#'   The final look is taken to be at 1.
#' @param alpha0 One-sided type I error to spend. Default \code{0.025}.
#' @param sfu Upper (efficacy) spending function. Default
#'   \code{gsDesign::sfLDOF}, which gives a Lan-DeMets O'Brien-Fleming–type
#'   spending.
#'
#' @return A numeric vector of length 2: \code{c(C1, C2)}, where
#'   \code{C1} is the interim efficacy boundary and \code{C2} is the final
#'   efficacy boundary (both on the Z-scale).
#'
#' @examples
#' efficacy_bounds <- get_efficacy_boundaries(timing = 0.5)
#' efficacy_bounds
#'
#' @seealso \code{\link[gsDesign]{gsDesign}}, \code{\link[gsDesign]{sfLDOF}}
#'
#' @importFrom gsDesign gsDesign sfLDOF
#' @export
get_efficacy_boundaries = function(timing, # this should be the actual timing
                                   alpha0 = .025,
                                   sfu = sfLDOF){
  # check timing: single number, (0,1)
  if (length(timing) != 1L) {
    stop("`timing` must be a single number in (0, 1).", call. = FALSE)
  }
  if (timing <= 0 || timing >= 1) {
    stop("`timing` must be in (0, 1).", call. = FALSE)
  }

  x <- gsDesign(
    k = 2,
    test.type = 1,
    alpha = alpha0,
    timing = timing,
    sfu = sfu
  )
  C1 = x$upper$bound[1]
  C2 = x$upper$bound[2]
  return(c(C1,C2))
}


