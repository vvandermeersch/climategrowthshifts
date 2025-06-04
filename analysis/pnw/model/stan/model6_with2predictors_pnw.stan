functions {
  vector gp_pred_rng(array[] real x2,
                     vector y1, array[] real x1,
                     vector mu, real alpha, real rho,
                     real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K = gp_exp_quad_cov(x1, alpha, rho)
      + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);
      
      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1 - mu);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = gp_exp_quad_cov(x1, x2, alpha, rho);
      vector[N2] f2_mu = mu + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 = gp_exp_quad_cov(x2, alpha, rho)
                              - v_pred' * v_pred + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
  
  vector gp2_pred_rng(array[] real x2,
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
  int<lower=1> N;           // Total number of observations
  int<lower=1> N_trees;     // Number of trees
  int<lower=1> N_all_years; // Number of years spanned by all observations
  int<lower=1> N_stands;    // NUmber of stands
  
  // Indices of tree stands
  array[N_trees] int<lower=1, upper=N_stands> stand_idxs;
  
  // Number of observations per tree
  array[N_trees] int<lower=1, upper=N> N_years;
  
  // Years spanned by all observations
  array[N_all_years] real all_years;
  
  // Tree observation years
  array[N] real years;
  
  // Indices of tree observation years
  array[N] int<lower=1, upper=N_all_years> all_years_idxs;
  
  // Ragged array indexing
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;
  
  vector[N] log_rw_obs; // Log of Observed Ring Width Per 1 mm
  vector[N] gdd_obs;    // Observed gdd (during GS) (x10 degC)
  vector[N] sm_obs;    // Observed soil moisture (during GS) (m3.m-3)
  
  
}

transformed data {
  real gdd0 = 100;
  real sm0 = 25;
}

parameters {
  real alpha; // Log ring width baseline
  real<lower=0> beta_gdd; // GDD slope 
  real<lower=0> beta_sm; // SM slope 
  
  real<lower=0> rho;   // Lifetime growth time scale
  real<lower=0> gamma; // Lifetime proportional growth variation
  
  // Short-term proportional growth functional behavior
  // for species Abam
  array[N_stands] vector[N_all_years] f_tilde_sh; // Non-centered functional behavior
  real<lower=3> rho_sh;   // Time scale
  real<lower=0> gamma_sh; // Marginal variation

  real<lower=0> sigma; // Proportional measurement error
}

transformed parameters {
  
  array[N_stands] vector[N_all_years] f_sh;
  {
    matrix[N_all_years, N_all_years] cov
    =   gp_exp_quad_cov(all_years, gamma_sh, rho_sh)
    + diag_matrix(rep_vector(1e-10, N_all_years));
    matrix[N_all_years, N_all_years] L_cov = cholesky_decompose(cov);
    
    for (s in 1:N_stands) {
      f_sh[s] = L_cov * f_tilde_sh[s];
    }
  }
}

model {
  matrix[N_all_years, N_all_years] L_cov;
  matrix[N_all_years, N_all_years] cov
    =  gp_exp_quad_cov(all_years, gamma, rho)
    + diag_matrix(rep_vector(square(sigma), N_all_years));
  L_cov = cholesky_decompose(cov);
  
  
  alpha ~ normal(0, 0.69); // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  beta_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  beta_sm ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  
  rho ~ lognormal(3.55, 0.24);       // 20 <~ rho <~ 60
  gamma ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)
  
  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26);       // 3 <~ rho_sh <~ 10
  gamma_sh ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)
  
  sigma ~ normal(0, 0.095 / 2.57);   // -log(1.1) <~ sigma <~ +log(1.1)
  
  for (t in 1:N_trees) {
    array[N_years[t]] int tree_idxs
    = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);
    
    array[N_years[t]] int all_years_idxs_tree
    = all_years_idxs[tree_idxs];
    array[N_years[t]] real years_tree = years[tree_idxs];
    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    
    int stand_idx = stand_idxs[t];
    
    vector[N_years[t]] mu =  alpha
    + beta_gdd * (gdd_obs_tree - gdd0)
    + beta_sm * (sm_obs_tree - sm0)
    + f_sh[stand_idx, all_years_idxs_tree];
    
    // ***** WARNING *****
      // Making the assumption that N_years[t] = N_all_years for all t
      // so that L_cov can be used for all trees.  For ragged tree
      // observations L_cov needs to be marginalized to the specific
      // observation years for each tree.
      log_rw_obs[tree_idxs]
      ~ multi_normal_cholesky(mu, L_cov);
  }
}

generated quantities {
  vector[N] mu0 = rep_vector(alpha, N);
  vector[N] mu1;
  vector[N] rw;
  array[N] real log_rw_pred;
  array[N] real log_rw_ltonly_pred;
  array[N] real log_rw_woclim_pred;
  
  for (t in 1:N_trees) {
    array[N_years[t]] int tree_idxs
    = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);
    
    array[N_years[t]] int all_years_idxs_tree
    = all_years_idxs[tree_idxs];
    array[N_years[t]] real years_tree = years[tree_idxs];
    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    
    int stand_idx = stand_idxs[t];
    
     vector[N_years[t]] mu =  alpha
    + beta_gdd * (gdd_obs_tree - gdd0)
    + beta_sm * (sm_obs_tree - sm0)
    + f_sh[stand_idx, all_years_idxs_tree];
    
    vector[N_years[t]] mu_ltonly;
    
    for(i in 1:N_years[t]) {
      
      mu_ltonly[i] = alpha;
      
    }
    
    vector[N_years[t]] mu_woclim = alpha
    + f_sh[stand_idx, all_years_idxs_tree];
    
    vector[N_years[t]] log_rw = gp_pred_rng(years_tree,
                                            log_rw_obs[tree_idxs],
                                            years_tree,
                                            mu,
                                            gamma,
                                            rho,
                                            sigma,
                                            1e-10);
                                            
    vector[N_years[t]] log_rw_ltonly = gp_pred_rng(years_tree,
                                            log_rw_obs[tree_idxs],
                                            years_tree,
                                            mu_ltonly,
                                            gamma,
                                            rho,
                                            sigma,
                                            1e-10);
                                            
    vector[N_years[t]] log_rw_woclim = gp_pred_rng(years_tree,
                                            log_rw_obs[tree_idxs],
                                            years_tree,
                                            mu_woclim,
                                            gamma,
                                            rho,
                                            sigma,
                                            1e-10);                                           
                                            
    rw[tree_idxs] = exp(log_rw);
    log_rw_pred[tree_idxs] = normal_rng(log_rw, sigma);
    log_rw_ltonly_pred[tree_idxs] = normal_rng(log_rw_ltonly, sigma);
    log_rw_woclim_pred[tree_idxs] = normal_rng(log_rw_woclim, sigma);
    
    mu1[tree_idxs] = beta_gdd * (gdd_obs_tree - gdd0)
    + beta_sm * (sm_obs_tree - sm0);
  }
}
