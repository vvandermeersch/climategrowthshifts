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
  int<lower=1> N_all_years; // Max number of years
  int<lower=1> N_stands;    // Number of stands
  
  // Indices of tree stands
  array[N_trees] int<lower=1, upper=N_stands> stand_idxs;
  
  // Number of observations per tree
  array[N_trees] int<lower=1, upper=N> N_years;
  
  // Years spanned by all observations
  array[N_all_years] real all_years;
  
  // Number of year observed in each stand
  // array[N_stands] int<lower=1, upper=N_all_years> N_years_stand;
  
  // Tree observation years
  array[N] real years;
  
  // Indices of tree observation years
  array[N] int<lower=1, upper=N_all_years> all_years_idxs;
  
  // Ragged array indexing
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;
  
  vector[N] log_rw_obs; // Log of Observed Ring Width Per 1 mm
  vector[N] gdd_obs; // Observed gdd (all year) (x100 degC)
  vector[N] sm_obs; // Observed soil moisture (MJJ) (%)
  vector[N] vpd_obs; // Observed VPD (MJJ) (hPa)
  vector[N] pre_obs; // Observed summer precipitation (JJA) (cm)
  
  // corr_matrix[N_species] Cphy; // phylogenetic relationship matrix (fixed)
  
}

transformed data {
  real gdd0 = 10;
  real vpd0 = 8;
  real pre0 = 0;
}

parameters {
  real alpha; // Log ring width baseline
  
  // GDD slope (1/100degC)
  real beta_gdd;
  
  // VPD slope (1/hPa)
  real beta_vpd;
  
  // Precipitation slope (1/cm)
  real beta_pre;
  
  // Lifetime proportional growth scale (here I implement the hard contraint on both clade and species parameters?)
  real<lower=0> rho_sp; 
  
  // Lifetime proportional growth variation
  real<lower=0> gamma_sp; 
  
  real<lower=0> sigma; // Proportional measurement error
}

model {
  
  matrix[N_all_years, N_all_years] L_cov;
  matrix[N_all_years, N_all_years] cov
    =  gp_exp_quad_cov(all_years, gamma_sp, rho_sp)
    + diag_matrix(rep_vector(square(sigma), N_all_years));
  L_cov = cholesky_decompose(cov);
  
  
  alpha ~ normal(0, log(5)/2.32); // 0.2 mm < exp(alpha) * 1 mm < 5 mm
  
  beta_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) < beta_gsl < log(1.8)
  beta_pre ~ normal(0, log(1.8) / 2.57); // -log(1.8) < beta_pre < log(1.8)
  beta_vpd ~ normal(0, log(1.8) / 2.57); // -log(1.8) < beta_vpd < log(1.8)
  
  rho_sp ~ lognormal(3.7, 0.35); // 20 < rho < 90
  gamma_sp ~ normal(0, log(10) / 2.57); // 0 < gamma < log(10)
  
  sigma ~ normal(0, log(1.1) / 2.57); // 0 < sigma < log(1.1)
  
  for (t in 1:N_trees) {
    array[N_years[t]] int tree_idxs
    = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);
    
    array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
    array[N_years[t]] real years_tree = years[tree_idxs];
    
    int stand_idx = stand_idxs[t];
    
    // L_cov needs to be marginalized to the specific
    // observation years for each tree.
    // matrix[N_years[t], N_years[t]] cov
    // =  gp_exp_quad_cov(years_tree, gamma_sp[species_idx], rho_sp[species_idx])
    //  + diag_matrix(rep_vector(square(sigma), N_years[t]));
    // matrix[N_years[t], N_years[t]] L_cov = cholesky_decompose(cov);
    matrix[N_years[t], N_years[t]] L_cov_tree = block(L_cov, 1, 1, N_years[t], N_years[t]);
    
    
    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];
    vector[N_years[t]] pre_obs_tree = pre_obs[tree_idxs];
    
    vector[N_years[t]] mu =  alpha
    + beta_gdd * (gdd_obs_tree - gdd0)
    + beta_vpd * (vpd_obs_tree - vpd0)
    + beta_pre * (pre_obs_tree - pre0);
    
    log_rw_obs[tree_idxs] ~ multi_normal_cholesky(mu, L_cov_tree);
    
  }
}

generated quantities {

  array[N] real log_rw_pred;
  
  for (t in 1:N_trees) {
    
    array[N_years[t]] int tree_idxs = linspaced_int_array(N_years[t], tree_start_idxs[t], tree_end_idxs[t]);
    array[N_years[t]] real years_tree = years[tree_idxs];
    array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
    
    int stand_idx = stand_idxs[t];
    
    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];
    vector[N_years[t]] pre_obs_tree = pre_obs[tree_idxs];
    
    vector[N_years[t]] mu1 = alpha 
    + beta_gdd * (gdd_obs_tree - gdd0) 
    + beta_vpd * (vpd_obs_tree - vpd0) 
    + beta_pre * (pre_obs_tree - pre0);
    
    vector[N_years[t]] log_rw = gp_pred_rng(years_tree, log_rw_obs[tree_idxs], years_tree,
                                mu1, gamma_sp, rho_sp, sigma, 1e-10);
                                
    // rw[tree_idxs] = exp(log_rw[tree_idxs]);
    log_rw_pred[tree_idxs] = normal_rng(log_rw, sigma);
                                 
  }
}
