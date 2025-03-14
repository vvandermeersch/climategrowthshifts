data {
  int<lower=1> N;                 // Number of observations
  array[N] real<lower=1919> year; // Observation year
  vector[N] log_rw_obs;           // Log of Ring Width Per 1 mm
}

parameters {
  real alpha;          // Log ring width baseline
  real<lower=0> rho;   // Lifetime growth time scale
  real<lower=0> gamma; // Lifetime proportional growth variation
  real<lower=0> sigma; // Proportional measurement error
}

model {
  matrix[N, N] cov =  gp_exp_quad_cov(year, gamma, rho)
                    + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);

  alpha ~ normal(0, 0.69);            // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  rho ~ lognormal(3.55, 0.24);       // log(20) <~ rho <~ log(60)
  gamma ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)
  sigma ~ normal(0, 0.095 / 2.57);   // -log(1.1) <~ sigma <~ +log(1.1)
}

generated quantities {
  vector[N] rw;
  array[N] real log_rw_pred;
  {
    matrix[N, N] cov =  gp_exp_quad_cov(year, gamma, rho)
                      + diag_matrix(rep_vector(1e-10, N));
    matrix[N, N] L_cov = cholesky_decompose(cov);
    vector[N] log_rw
      = multi_normal_cholesky_rng(rep_vector(alpha, N), L_cov);
    rw = exp(log_rw);
    log_rw_pred = normal_rng(log_rw, sigma);
  }
}
