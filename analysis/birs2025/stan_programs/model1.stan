functions {
  vector gp_pred_rng(array[] real x2,
                     vector y1, array[] real x1,
                     vector mu,
                     real alpha, real rho,
                     real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   gp_exp_quad_cov(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1 - mu);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = gp_exp_quad_cov(x1, x2, alpha, rho);
      vector[N2] f2_mu = mu + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   gp_exp_quad_cov(x2, alpha, rho)
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

  log_rw_obs ~ multi_normal_cholesky(rep_vector(alpha, N), L_cov);
}

generated quantities {
  vector[N] rw;
  array[N] real log_rw_pred;
  {
    vector[N] log_rw = gp_pred_rng(year, log_rw_obs, year,
                                   mu, gamma, rho, sigma, 1e-10);
    rw = exp(log_rw);
    log_rw_pred = normal_rng(log_rw, sigma);
  }
}
