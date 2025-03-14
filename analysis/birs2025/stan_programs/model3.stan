functions {
  vector gp_pred_rng(array[] real x2,
                     vector y1, array[] real x1,
                     vector mu,
                     real alpha1, real rho1,
                     real alpha2, real rho2,
                     real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   gp_exp_quad_cov(x1, alpha1, rho1)
                         + gp_exp_quad_cov(x1, alpha2, rho2)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1 - mu);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 =   gp_exp_quad_cov(x1, x2, alpha1, rho1)
                               + gp_exp_quad_cov(x1, x2, alpha2, rho2);
      vector[N2] f2_mu = mu + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   gp_exp_quad_cov(x2, alpha1, rho1)
                              + gp_exp_quad_cov(x2, alpha2, rho2)
                              - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N;                 // Number of observations
  array[N] real<lower=1919> year; // Observation year
  vector[N] log_rw_obs;           // Log of Ring Width Per 1 mm
  vector[N] gsl;                  // Growing season length (days)
}

transformed data {
  real gsl0 = 150; // GSL Baseline (days)
}

parameters {
  real alpha;             // Log ring width baseline
  real<lower=0> beta_gsl; // GSL slope (1 / days)

  real<lower=0> rho;   // Lifetime growth time scale
  real<lower=0> gamma; // Lifetime proportional growth variation

  real<lower=0> rho_sh;   // Short-term growth time scale
  real<lower=0> gamma_sh; // Short-term proportional growth variation

  real<lower=0> sigma; // Proportional measurement error
}

model {
  vector[N] mu = alpha + beta_gsl * (gsl - gsl0);

  matrix[N, N] cov =  gp_exp_quad_cov(year, gamma, rho)
                    + gp_exp_quad_cov(year, gamma_sh, rho_sh)
                    + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);

  alpha ~ normal(0, 0.69); // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  beta_gsl ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)

  rho ~ lognormal(3.55, 0.24);       // log(20) <~ rho <~ log(60)
  gamma ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)

  rho_sh ~ lognormal(1.7, 0.26);       // log(3) <~ rho_sh <~ log(10)
  gamma_sh ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)

  sigma ~ normal(0, 0.095 / 2.57);   // -log(1.1) <~ sigma <~ +log(1.1)

  log_rw_obs ~ multi_normal_cholesky(mu, L_cov);
}

generated quantities {
  vector[N] mu0 = rep_vector(alpha, N);
  vector[N] mu1 = beta_gsl * (gsl - gsl0);
  vector[N] rw;
  array[N] real log_rw_pred;
  {
    vector[N] log_rw = gp_pred_rng(year, log_rw_obs, year,
                                   mu0 + mu1,
                                   gamma, rho, gamma_sh, rho_sh,
                                   sigma, 1e-10);
    rw = exp(log_rw);
    log_rw_pred = normal_rng(log_rw, sigma);
  }
}
