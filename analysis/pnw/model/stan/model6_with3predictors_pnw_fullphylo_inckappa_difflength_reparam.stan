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
  int<lower=1> N_species;   // Number of species
  int<lower=1> N_trees;     // Number of trees
  int<lower=1> N_all_years; // Max number of years
  int<lower=1> N_stands;    // Number of stands
  
  // Indices of tree stands
  array[N_trees] int<lower=1, upper=N_stands> stand_idxs;
  
  // Indices of tree species
  array[N_trees] int<lower=1, upper=N_species> species_idxs;
  
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
  vector[N] gdd_obs; // Observed gdd (all year) (x100 degC)
  vector[N] sm_obs; // Observed soil moisture (MJJ) (%)
  vector[N] vpd_obs; // Observed VPD (JJMJJA) (hPa)
  
  corr_matrix[N_species] Cphy; // phylogenetic relationship matrix (fixed)
  
}

transformed data {
  real gdd0 = 10;
  real sm0 = 25;
  real vpd0 = 8;
}

parameters {
  real alpha; // Log ring width baseline
  
  // GDD slope (1/100degC)
  vector[N_species] beta_tilde_gdd_sp; 
  real beta_gdd; // root value
  real<lower=0, upper=1> lambda_beta_gdd;  // phylogenetic structure        
  real<lower=0> sigma_beta_gdd; // overall rate of change (brownian motion?)
  
  // SM slope (1/%)
  vector[N_species] beta_tilde_sm_sp; 
  real beta_sm; // root value
  real<lower=0, upper=1> lambda_beta_sm;  // phylogenetic structure        
  real<lower=0> sigma_beta_sm; // overall rate of change (brownian motion?)
  
  // VPD slope (1/hPa)
  vector[N_species] beta_tilde_vpd_sp; 
  real beta_vpd; // root value
  real<lower=0, upper=1> lambda_beta_vpd;  // phylogenetic structure        
  real<lower=0> sigma_beta_vpd; // overall rate of change (brownian motion?)

  // Lifetime proportional growth scale
  // array[N_species] real<lower=0> rho_sp;
  vector<lower=0>[N_species] rho_sp; 
  real<lower=0> rho; // root value
  real<lower=0, upper=1> lambda_rho;  // phylogenetic structure        
  real<lower=0> sigma_rho; // overall rate of change (brownian motion?)
  
  // Lifetime proportional growth variation
  vector<lower=0>[N_species] gamma_tilde_sp; 
  real<lower=0> gamma; // root value
  real<lower=0, upper=1> lambda_gamma;  // phylogenetic structure        
  real<lower=0> sigma_gamma; // overall rate of change (brownian motion?)
  
  
  // Short-term proportional growth functional behavior
  // for species JUOC
  array[N_stands] vector[N_all_years] f_tilde_sh; // Non-centered functional behavior
  real<lower=0> rho_sh;   // Time scale
  // real<lower=0> gamma_sh; // Marginal variation - NOW FIXED TO 1!
  
  // Short-term proportional growth species scaling
  vector<lower=0>[N_species] kappa_sh_sp; 
  real<lower=0> kappa_sh; // root value
  real<lower=0, upper=1> lambda_kappa_sh;  // phylogenetic structure        
  real<lower=0> sigma_kappa_sh; // overall rate of change (brownian motion?)
  
  real<lower=0> sigma; // Proportional measurement error
}

transformed parameters {
  // array[N_species] real kappa_sh = append_array({1}, kappa_sh_free);
  
  array[N_stands] vector[N_all_years] f_sh;
  {
    matrix[N_all_years, N_all_years] cov
    =   gp_exp_quad_cov(all_years, 1, rho_sh)
    + diag_matrix(rep_vector(1e-10, N_all_years));
    matrix[N_all_years, N_all_years] L_cov = cholesky_decompose(cov);
    
    for (s in 1:N_stands) {
      f_sh[s] = L_cov * f_tilde_sh[s];
    }
  }
  
  corr_matrix[N_species] C_beta_gdd = lambda_beta_gdd * Cphy;
  C_beta_gdd = C_beta_gdd - diag_matrix(diagonal(C_beta_gdd)) + diag_matrix(diagonal(Cphy));
  corr_matrix[N_species] C_beta_sm = lambda_beta_sm * Cphy;
  C_beta_sm = C_beta_sm - diag_matrix(diagonal(C_beta_sm)) + diag_matrix(diagonal(Cphy));
  corr_matrix[N_species] C_beta_vpd = lambda_beta_vpd * Cphy;
  C_beta_vpd = C_beta_vpd - diag_matrix(diagonal(C_beta_vpd)) + diag_matrix(diagonal(Cphy));
  
  corr_matrix[N_species] C_rho = lambda_rho * Cphy;
  C_rho = C_rho - diag_matrix(diagonal(C_rho)) + diag_matrix(diagonal(Cphy));
  corr_matrix[N_species] C_gamma = lambda_gamma * Cphy;
  C_gamma = C_gamma - diag_matrix(diagonal(C_gamma)) + diag_matrix(diagonal(Cphy));
  
  corr_matrix[N_species] C_kappa_sh = lambda_kappa_sh * Cphy;
  C_kappa_sh = C_kappa_sh - diag_matrix(diagonal(C_kappa_sh)) + diag_matrix(diagonal(Cphy));
  
  // more numerically stable and more efficient to use pre-factored covariance matrices (i.e. multi_normal_cholesky in the following)
  matrix[N_species,N_species] L_beta_gdd =  cholesky_decompose(sigma_beta_gdd^2*C_beta_gdd); 
  matrix[N_species,N_species] L_beta_sm =  cholesky_decompose(sigma_beta_sm^2*C_beta_sm); 
  matrix[N_species,N_species] L_beta_vpd =  cholesky_decompose(sigma_beta_vpd^2*C_beta_vpd); 
  
  matrix[N_species,N_species] L_rho =  cholesky_decompose(sigma_rho^2*C_rho); 
  matrix[N_species,N_species] L_gamma =  cholesky_decompose(sigma_gamma^2*C_gamma); 
  
  matrix[N_species,N_species] L_kappa_sh =  cholesky_decompose(sigma_kappa_sh^2*C_kappa_sh); 
  
  // Reparametrization
  vector[N_species] beta_gdd_sp = beta_gdd + L_beta_gdd*beta_tilde_gdd_sp;
  vector[N_species] beta_sm_sp = beta_sm + L_beta_sm*beta_tilde_sm_sp;
  vector[N_species] beta_vpd_sp = beta_vpd + L_beta_vpd*beta_tilde_vpd_sp;
  vector<lower=0>[N_species] gamma_sp = gamma + L_gamma*gamma_tilde_sp;; 
}

model {
  
  array[N_species] matrix[N_all_years, N_all_years] L_cov;
  for (sp in 1:N_species) {
    matrix[N_all_years, N_all_years] cov
    =  gp_exp_quad_cov(all_years, gamma_sp[sp], rho_sp[sp])
    + diag_matrix(rep_vector(square(sigma), N_all_years));
    L_cov[sp] = cholesky_decompose(cov);
  }
  
  alpha ~ normal(0, 0.69); // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  
  beta_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  beta_sm ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  beta_vpd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  
  // Temporary priors (to discuss)
  lambda_beta_gdd ~ beta(1.5, 1.5);
  sigma_beta_gdd  ~ normal(0, 1);
  lambda_beta_sm ~ beta(1.5, 1.5);
  sigma_beta_sm  ~ normal(0, 1);
  lambda_beta_vpd ~ beta(1.5, 1.5);
  sigma_beta_vpd  ~ normal(0, 1);
  lambda_rho ~ beta(1.5, 1.5);
  sigma_rho  ~ normal(0, 1);
  lambda_gamma ~ beta(1.5, 1.5);
  sigma_gamma  ~ normal(0, 1);
  lambda_kappa_sh ~ beta(1.5, 1.5);
  sigma_kappa_sh  ~ normal(0, 1);
  
  // beta_gdd_sp ~ multi_normal_cholesky(rep_vector(beta_gdd,N_species), L_beta_gdd); 
  // beta_sm_sp ~ multi_normal_cholesky(rep_vector(beta_sm,N_species), L_beta_sm); 
  // beta_vpd_sp ~ multi_normal_cholesky(rep_vector(beta_vpd,N_species), L_beta_vpd); 
  beta_tilde_gdd_sp ~ std_normal();
  beta_tilde_sm_sp ~ std_normal();
  beta_tilde_vpd_sp ~ std_normal();
  
  rho_sp ~ multi_normal_cholesky(rep_vector(rho,N_species), L_rho); 
  // gamma_sp ~ multi_normal_cholesky(rep_vector(gamma,N_species), L_gamma); 
  gamma_tilde_sp ~ std_normal();
  
  rho ~ lognormal(3.55, 0.24);       // 20 <~ rho <~ 60
  gamma ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)
  
  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26);       // 3 <~ rho_sh <~ 10
  // gamma_sh ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)
  
  kappa_sh ~ lognormal(1, 0.41 / 2.32); // 2/3 <~ kappa_sh <~ 3/2
  kappa_sh_sp ~ multi_normal_cholesky(rep_vector(kappa_sh,N_species), L_kappa_sh); 
  
  sigma ~ normal(0, 0.095 / 2.57);   // -log(1.1) <~ sigma <~ +log(1.1)
  
  for (t in 1:N_trees) {
    array[N_years[t]] int tree_idxs
    = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);
    
    array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
    array[N_years[t]] real years_tree = years[tree_idxs];
    
    int stand_idx = stand_idxs[t];
    int species_idx = species_idxs[t];
    
    // L_cov needs to be marginalized to the specific
    // observation years for each tree.
    // matrix[N_years[t], N_years[t]] cov
    // =  gp_exp_quad_cov(years_tree, gamma_sp[species_idx], rho_sp[species_idx])
    //  + diag_matrix(rep_vector(square(sigma), N_years[t]));
    // matrix[N_years[t], N_years[t]] L_cov = cholesky_decompose(cov);
    matrix[N_years[t], N_years[t]] L_cov_tree = block(L_cov[species_idx], 1, 1, N_years[t], N_years[t]);
    
    
    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];
    
    vector[N_years[t]] mu =  alpha
    + beta_gdd_sp[species_idx] * (gdd_obs_tree - gdd0)
    + beta_sm_sp[species_idx] * (sm_obs_tree - sm0)
    + beta_vpd_sp[species_idx] * (vpd_obs_tree - vpd0)
    + kappa_sh_sp[species_idx]
    * f_sh[stand_idx, all_years_idxs_tree];
    
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
    int species_idx = species_idxs[t];
    
    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];
    
    vector[N_years[t]] mu1 = alpha + beta_gdd_sp[species_idx] * (gdd_obs_tree - gdd0) + beta_sm_sp[species_idx] * (sm_obs_tree - sm0)
    + beta_vpd_sp[species_idx] * (vpd_obs_tree - sm0) + kappa_sh_sp[species_idx] * f_sh[stand_idx, all_years_idxs_tree];
    
    vector[N_years[t]] log_rw = gp_pred_rng(years_tree, log_rw_obs[tree_idxs], years_tree,
                                mu1, gamma_sp[species_idx], rho_sp[species_idx], sigma, 1e-10);
                                
    // rw[tree_idxs] = exp(log_rw[tree_idxs]);
    log_rw_pred[tree_idxs] = normal_rng(log_rw, sigma);
                                 
  }
}
