#' Simulation study for sample size re-estimation using a totality of evidence approach
#'
#' This runs \code{nsim} simulated two-stage trials. For each replicate it:
#' 1) optionally computes the *initial* (planned) sample size using
#'    \code{get_initial_tot_ss()},
#' 2) generates stage-1 data,
#' 3) runs \code{interim_analysis()},
#' 4) if the trial stops for efficacy, records that and moves on,
#' 5) otherwise generates stage-2 data according to the re-estimated sizes and
#'    runs \code{final_analysis()},
#' 6) records whether H0 was rejected (can be rejected at either the interim or
#'    final analysis).
#'
#' The returned object contains per-replicate interim and final summaries and
#' the overall proportion of rejections.
#'
#' @details
#' This function calls the following in-package functions:
#' - \code{simulate_example_one_stage_data()}
#' - \code{get_initial_tot_ss()}
#' - \code{interim_analysis()}
#' - \code{final_analysis()}
#'
#' @param nsim Number of simulated trials.
#' @param alpha One-sided alpha
#' @param beta Type II error rate (\eqn{\beta}). (for power / CP).
#' @param n_endpts Number of endpoints.
#' @param allocate_rate Treatment:control allocation ratio.
#' @param mean_cohen_d_truth True mean Cohen's d across endpoints used to simulate data.
#' @param mu_ct Mean vector for the control group (length = \code{n_endpts}).
#' @param sd Standard deviations for each endpoint (length = \code{n_endpts}).
#' @param rho Correlation matrix for the endpoints.
#' @param timing_interim Planned interim timing (information fraction).
#' @param initial_tot_ss Logical. If \code{TRUE}, compute initial sample sizes
#'   via \code{get_initial_tot_ss()}; if \code{FALSE}, use \code{tot_ss_ct}.
#' @param exp_mean_cohen_d_at_planning Expected mean Cohen's d at planning;
#'   required when \code{initial_tot_ss = TRUE}.
#' @param tot_ss_ct Optional fixed total control sample size
#'   (used only when \code{initial_tot_ss = FALSE}).
#' @param SSR_type \code{"SSR-Power"}, \code{"SSR-CP"}, or \code{"No SSR"}.
#' @param global_test_type \code{"exact OLS"} or \code{"Permutation"}.
#' @param nPM Number of permutations for the permutation test. Default: 100.
#' @param promising_LL Lower bound of the promising zone. Default: 0.2. (Upper bound: 1-beta0)
#' @param max_ss_index Maximum inflation factor for SSR
#'   (e.g., 2 means the sample size can at most be doubled relative to the
#'   original plan).
#' @param seed1 Base seed for stage-1 simulations.
#' @param seed2 Base seed for stage-2 simulations.
#' @param seed_PM Base seed for Permutation. Default seed: 1.
#'
#' @return A list with
#' \itemize{
#'   \item \code{interim_analysis}: tibble with one row per simulation
#'   \item \code{final_analysis}: tibble with one row per simulation that got to final
#'   \item \code{proportion_of_rejection}: overall rejection proportion (i.e., type I error or power)
#' }
#'
#' @examples
#' nsim <- 1e3
#' alpha <- 0.025
#' beta  <- 0.2
#' n_endpts <- 6
#' allocate_rate <- 1
#' mean_cohen_d_truth <- 0.25
#' timing_interim <- 2/3
#' max_ss_index <- 2
#'
#' mu_ct <- c(4, 6, 4, 5, 7, 6)
#' sd <- c(0.5, 0.5, 1, 1, 2, 2)
#' rho <- matrix(c(
#'   1.0, 0.1, 0.3, 0.7, 0.1, 0.3,
#'   0.1, 1.0, 0.7, 0.1, 0.3, 0.7,
#'   0.3, 0.7, 1.0, 0.1, 0.3, 0.7,
#'   0.7, 0.1, 0.1, 1.0, 0.1, 0.3,
#'   0.1, 0.3, 0.3, 0.1, 1.0, 0.7,
#'   0.3, 0.7, 0.7, 0.3, 0.7, 1.0
#' ), nrow = 6, byrow = TRUE)
#'
#' SSR_type <- "SSR-Power"
#' global_test_type <- "exact OLS"
#'
#' seed1 <- 51
#' seed2 <- 82
#'
#' # at planning, we thought exp_mean_cohen_d = 0.3
#' exp_mean_cohen_d_at_planning <- 0.3
#'
#' res <- SSRTE_simstudy(
#'   nsim = nsim,
#'   alpha = alpha,
#'   beta = beta,
#'   n_endpts = n_endpts,
#'   allocate_rate = allocate_rate,
#'   mean_cohen_d_truth = mean_cohen_d_truth,
#'   mu_ct = mu_ct,
#'   sd = sd,
#'   rho = rho,
#'   timing_interim = timing_interim,
#'   initial_tot_ss = TRUE,
#'   exp_mean_cohen_d_at_planning = exp_mean_cohen_d_at_planning,
#'   SSR_type = SSR_type,
#'   global_test_type = global_test_type,
#'   max_ss_index = max_ss_index,
#'   seed1 = seed1,
#'   seed2 = seed2
#' )
#'
#' res
#'
#' @export
SSRTE_simstudy <- function(
    nsim,
    alpha, beta,
    n_endpts,
    allocate_rate,
    mean_cohen_d_truth,
    mu_ct,
    sd,
    rho,
    timing_interim,
    initial_tot_ss = TRUE,
    exp_mean_cohen_d_at_planning = NULL,
    tot_ss_ct = NULL,
    SSR_type,
    global_test_type,
    max_ss_index,
    seed1,
    seed2,
    seed_PM = 1,
    nPM = 1e2,
    promising_LL = 0.2
) {
  ia_out <- tibble::tibble()
  fa_out <- tibble::tibble()
  rej_ind <- rep(NA_integer_, nsim)

  for (si in seq_len(nsim)) {

    # ------------------- 1. get planned sample sizes -------------------------
    if (isTRUE(initial_tot_ss)) {
      if (is.null(exp_mean_cohen_d_at_planning)) {
        stop("Must provide exp_mean_cohen_d_at_planning when initial_tot_ss = TRUE.",
             call. = FALSE)
      }

      init_ss <- get_initial_tot_ss(
        beta0    = beta,
        r        = allocate_rate,
        theta_k  = exp_mean_cohen_d_at_planning,
        n_endpts = n_endpts,
        alpha0   = alpha,
        rho      = rho,
        timing0  = timing_interim
      )

      N_ct <- init_ss$N_ct
      n_ct <- init_ss$n_ct

      N_trt <- init_ss$N_trt
      n_trt <- init_ss$n_trt

      timing_actual <- init_ss$timing_actual

    } else {
      # user supplies total control ss
      if (is.null(tot_ss_ct)) {
        stop("initial_tot_ss = FALSE, but tot_ss_ct is NULL.", call. = FALSE)
      }
      N_ct <- as.integer(tot_ss_ct)
      n_ct <- round(N_ct * timing_interim)

      N_trt <- ceiling(N_ct * allocate_rate)
      n_trt <- ceiling(n_ct * allocate_rate)

      timing_actual <- (n_ct + n_trt) / (N_ct + N_trt)
    }

    # ------------------- 2. stage 1 data -------------------------------------
    sim1 <- simulate_example_one_stage_data(
      n_trt = n_trt,
      n_ct  = n_ct,
      n_endpts = n_endpts,
      exp_mean_cohen_d = mean_cohen_d_truth,
      mu_ct = mu_ct,
      sd = sd,
      rho = rho,
      seed = seed1 + si
    )

    y_ct1  <- sim1$y_ct
    y_trt1 <- sim1$y_trt

    # ------------------- 3. interim ------------------------------------------
    if (global_test_type=="Permutation"){
      set.seed(seed_PM)
    }
    ia <- interim_analysis(
      y_trt1 = y_trt1,
      y_ct1  = y_ct1,
      N_trt = N_trt,
      N_ct = N_ct,
      N_ct_max = max_ss_index * N_ct,
      timing_actual = timing_actual,
      alpha0 = alpha,
      beta0 = beta,
      alloc_rate = allocate_rate,
      n_endpts = n_endpts,
      SSR_type = SSR_type,
      global_test_type = global_test_type,
      nPM = nPM,
      promising_LL = promising_LL
    )

    sum_ia <- dplyr::mutate(summary(ia), sim_id = si)
    ia_out <- dplyr::bind_rows(ia_out, sum_ia)

    # stop at interim ---------------------------------------------------------
    if (ia$decision == "Reject H0 at interim; stop the trial.") {
      rej_ind[si] <- 1L
      next
    }

    # ------------------- 4. stage 2 sample sizes from interim ----------------
    n_ct2  <- ia$sample_sizes$stage2["control"]
    n_trt2 <- ia$sample_sizes$stage2["treatment"]

    # ------------------- 5. stage 2 data -------------------------------------
    sim2 <- simulate_example_one_stage_data(
      n_trt = n_trt2,
      n_ct  = n_ct2,
      n_endpts = n_endpts,
      exp_mean_cohen_d = mean_cohen_d_truth,
      mu_ct = mu_ct,
      sd = sd,
      rho = rho,
      seed = seed2 + si
    )
    y_ct2  <- sim2$y_ct
    y_trt2 <- sim2$y_trt

    # ------------------- 6. final --------------------------------------------
    fa <- final_analysis(
      interim = ia,
      y_trt2  = y_trt2,
      y_ct2   = y_ct2,
      alloc_rate = allocate_rate
    )

    sum_fa <- dplyr::mutate(summary(fa), sim_id = si)
    fa_out <- dplyr::bind_rows(fa_out, sum_fa)

    if (fa$decision == "Reject H0 at final analysis.") {
      rej_ind[si] <- 1L
    } else {
      rej_ind[si] <- 0L
    }
  }

  proportion_of_rejection <- mean(rej_ind)

  list(
    interim_analysis = dplyr::relocate(ia_out, sim_id, .before = 1),
    final_analysis   = dplyr::relocate(fa_out, sim_id, .before = 1),
    proportion_of_rejection = proportion_of_rejection
  )
}
