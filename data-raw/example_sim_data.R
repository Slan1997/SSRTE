## code to prepare `sim1` and `sim2` datasets goes here
# data-raw/example_sim_data.R

sim1 <- simulate_example_one_stage_data(
  n_trt = 29, n_ct = 29, n_endpts = 6,
  exp_mean_cohen_d = 0.3,
  mu_ct = c(4,6,4,5,7,6),
  sd  = c(1.5,0.5,1.6,1,0.6,1.2),
  rho = matrix(c(
    1.0,0.1,0.3,0.5,0.1,0.3,
    0.1,1.0,0.5,0.1,0.3,0.5,
    0.3,0.5,1.0,0.1,0.3,0.5,
    0.5,0.1,0.1,1.0,0.1,0.3,
    0.1,0.3,0.3,0.1,1.0,0.5,
    0.3,0.5,0.5,0.3,0.5,1.0
  ), 6, byrow = TRUE),
  seed = 42
)

sim2 <- simulate_example_one_stage_data(
  n_trt = 61, n_ct = 61, n_endpts = 6,
  exp_mean_cohen_d = 0.3,
  mu_ct = c(4,6,4,5,7,6),
  sd  = c(1.5,0.5,1.6,1,0.6,1.2),
  rho = matrix(c(
    1.0,0.1,0.3,0.5,0.1,0.3,
    0.1,1.0,0.5,0.1,0.3,0.5,
    0.3,0.5,1.0,0.1,0.3,0.5,
    0.5,0.1,0.1,1.0,0.1,0.3,
    0.1,0.3,0.3,0.1,1.0,0.5,
    0.3,0.5,0.5,0.3,0.5,1.0
  ), 6, byrow = TRUE),
  seed = 58
)

usethis::use_data(sim1, sim2, overwrite = TRUE, compress = "xz")

