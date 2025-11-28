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
  // int<lower=1> N_clades;    // Number of clades - for now, Gymnosperms=1 or Angiosperms=2
  
  array[N_stands] int<lower=1, upper=N> N_stand_trees; // Number of trees per stand
  array[N_stands] int<lower=1, upper=N> N_stand_years; // Max. number of observed years per stand
  array[N_stands] int<lower=1, upper=N> stand_start_years_idxs; // Indice of first year observed per stand
  
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
  
  // Ragged array indexing for trees
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;
  
  // Ragged array indexing for stands
  array[N_stands] int<lower=1, upper=N> stand_trees_start_idxs;
  array[N_stands] int<lower=1, upper=N> stand_trees_end_idxs;
  
  vector[N] log_rw_obs; // Log of Observed Ring Width Per 1 mm
  vector[N] gdd_obs; // Observed gdd (all year) (x100 degC)
  vector[N] sm_obs; // Observed soil moisture (MJJ) (%)
  vector[N] vpd_obs; // Observed VPD (JJMJJA) (hPa)
  
  // Indices of species clades (Gymnosperms=1 or Angiosperms=2)
  // array[N_species] int<lower=1, upper=N_clades> clade_idxs;
  
  // Partial centering parameters for growth shocks
  // vector<lower=0, upper=1>[N] w_sck;
  
}

transformed data {
  real gdd0 = 10;
  real sm0 = 25;
  real vpd0 = 8;
}

parameters {
  real alpha; // Log ring width baseline
  
  vector[N] f_tilde;
  
  // GDD slope (1/100degC)
  // vector[N_clades] mu_gdd;
  // vector<lower=0>[N_clades] tau_gdd;
  vector[N_species] beta_gdd;
  
  // SM slope (1/%)
  // vector[N_clades] mu_sm;
  //vector<lower=0>[N_clades] tau_sm;
  vector[N_species] beta_sm;
  
  // VPD slope (1/hPa)
  // vector[N_clades] mu_vpd;
  // vector<lower=0>[N_clades] tau_vpd;
  vector[N_species] beta_vpd;
  
  // Short-term proportional growth functional behavior
  array[N_stands] vector[N_all_years] f_tilde_sh; // Non-centered functional behavior
  real<lower=1> rho_sh;   // Time scale
  real<lower=0> gamma_sh; // Marginal variation - now fixed to 1! (and scaled by kappa)
  
  // Lifetime proportional growth scale (here I implement the hard contraint on both clade and species parameters?)
  // vector<lower=rho_sh>[N_clades] mu_rho;
  // vector<lower=0>[N_clades] tau_rho;
  vector<lower=rho_sh>[N_species] rho_sp; 
  
  // Lifetime proportional growth variation
  // vector<lower=0>[N_clades] mu_gamma;
  // vector<lower=0>[N_clades] tau_gamma;
  vector<lower=0>[N_species] gamma_sp; 
  
  // Short-term proportional growth species scaling
  // vector<lower=0>[N_clades] mu_kappa;
  // vector<lower=0>[N_clades] tau_kappa;
  // vector<lower=0>[N_species] kappa_sh; 
  
  // Growth shocks
  // real<lower=0> inner_tau_sck; // Inner yearly log variation scale
  real<lower=0> outer_tau_sck; // Outer yearly log variation scale (the shocks!) Now it's a Cauchy scale!
  //real<lower=0> eta; // eta corresponds to sqrt(sigma^2 + tau_inner^2)
  // real<lower=0> nu; // nu corresponds to tau_outer^2 - tau_inner^2
  
  vector<lower=0>[N] outer_tau2_aux; // latent variances for Cauchy
  
  // Probability of shocks
  vector<lower=0, upper=1>[N_stands] phi_sck; // Probability of stand-level shock
  real<lower=0, upper=1> omega_conc_sck; // Probability of tree-level shock given stand in shock (concordant shock)
  real<lower=0, upper=omega_conc_sck> omega_nonconc_sck; // Probability of tree-level shock given stand NOT in shock (nonconcordant shock)
  
  // Growth shock species scaling
  // vector<lower=0>[N_clades] mu_kappa_sck;
  // vector<lower=0>[N_clades] tau_kappa_sck;
  // vector<lower=0>[N_species] kappa_sck; 
  
  // Proportional measurement error
  real<lower=0> sigma; 
}

transformed parameters {
  // array[N_species] real kappa_sh = append_array({1}, kappa_sh_free);
  
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
  
  vector[N] f;
  {
    array[N_species] matrix[N_all_years, N_all_years] L_cov;
    for (sp in 1:N_species) {
      matrix[N_all_years, N_all_years] cov
      =  gp_exp_quad_cov(all_years, gamma_sp[sp], rho_sp[sp])
      + diag_matrix(rep_vector(1e-10, N_all_years));
      
      L_cov[sp] = cholesky_decompose(cov);
    }
    for (t in 1:N_trees) {
      array[N_years[t]] int tree_idxs
      = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);
                          
      int species_idx = species_idxs[t];
      
      // L_cov needs to be marginalized to the specific
      // observation years for each tree.
      matrix[N_years[t], N_years[t]] L_cov_tree = block(L_cov[species_idx], 1, 1, N_years[t], N_years[t]);
      
      f[tree_idxs] = L_cov_tree * f_tilde[tree_idxs];
    }
    
  }
  
  // vector[N_species] beta_gdd;
  // vector[N_species] beta_sm;
  // vector[N_species] beta_vpd;
  // 
  // for (sp in 1:N_species) {
  //   beta_gdd[sp] = mu_gdd[clade_idxs[sp]] + tau_gdd[clade_idxs[sp]] * beta_gdd_tilde[sp];
  //   beta_sm[sp] = mu_sm[clade_idxs[sp]] + tau_sm[clade_idxs[sp]] * beta_sm_tilde[sp];
  //   beta_vpd[sp] = mu_vpd[clade_idxs[sp]] + tau_vpd[clade_idxs[sp]] * beta_vpd_tilde[sp];
  // }

}

model {
  
  alpha ~ normal(0, 0.69); // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  
  beta_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  beta_sm ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  beta_vpd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  
  rho_sp ~ lognormal(2.65, 0.135); // 10 <~ rho <~ 20
  
  gamma_sp ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)
  
  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26);       // 3 <~ rho_sh <~ 10
  gamma_sh ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)
  
  phi_sck ~ beta(2, 20); 
  omega_conc_sck ~ beta(230, 14); // 0.9 <~ omega_conc_sck <~ 0.97 (most trees, but not ALL trees)
  omega_nonconc_sck ~ beta(1, 20); // 0 <~ omega_conc_sck <~ 0.15 (should be rare, but... who knows?)
  
  sigma ~ normal(0, log(1.1) / 2.57);   // 0 <~ sigma <~ +log(1.1)
  
  vector[N] mu;
  for (t in 1:N_trees) {
    array[N_years[t]] int tree_idxs
    = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);

    array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
    array[N_years[t]] real years_tree = years[tree_idxs];

    int stand_idx = stand_idxs[t];
    int species_idx = species_idxs[t];

    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];

    mu[tree_idxs] =  alpha
    + beta_gdd[species_idx] * (gdd_obs_tree - gdd0)
    + beta_sm[species_idx] * (sm_obs_tree - sm0)
    + beta_vpd[species_idx] * (vpd_obs_tree - vpd0)
    + f_sh[stand_idx, all_years_idxs_tree];
    
    f_tilde[tree_idxs] ~ normal(0,1);
  }
  
  // Cauchy is a scale mixture of Gaussian density functions!
  outer_tau_sck ~ normal(0, 1/ 2.57); // outer_tau_sck is now the scale of Cauchy. With outer_tau_sck = 1, we get a Cauchy between -12 and 12
  outer_tau2_aux ~ inv_gamma(0.5, 0.5 * square(outer_tau_sck)); 
  
  // Observational model with marginalized shocks
  for (s in 1:N_stands) {
    
    array[N_stand_trees[s]] int stand_trees_idxs = linspaced_int_array(N_stand_trees[s], 
      stand_trees_start_idxs[s], stand_trees_end_idxs[s]);
    vector[N_stand_years[s]] log_p0 = rep_vector(0, N_stand_years[s]);
    vector[N_stand_years[s]] log_p1 = rep_vector(0, N_stand_years[s]);
    
    for(ts in 1:N_stand_trees[s]){
      int t = stand_trees_idxs[ts];
      array[N_years[t]] int tree_idxs = linspaced_int_array(N_years[t], tree_start_idxs[t], tree_end_idxs[t]);
      array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
      
      for(y in 1:N_years[t]) {
        int ys = all_years_idxs_tree[y]-stand_start_years_idxs[s]+1;
        
        log_p0[ys] += log_mix(omega_nonconc_sck,
                          normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[tree_idxs[y]] 
                          + f[tree_idxs[y]], sqrt(outer_tau2_aux[tree_idxs[y]] + sigma^2)),
                          normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[tree_idxs[y]] 
                          + f[tree_idxs[y]], sigma));
        log_p1[ys] += log_mix(omega_conc_sck,
                          normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[tree_idxs[y]] 
                          + f[tree_idxs[y]], sqrt(outer_tau2_aux[tree_idxs[y]] + sigma^2)),
                          normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[tree_idxs[y]] 
                          + f[tree_idxs[y]], sigma));
        
      }
    }
    
    for(y in 1:N_stand_years[s]) {
      target += log_mix(phi_sck[s], log_p1[y], log_p0[y]);
    }
  }
  
}

// generated quantities {
// 
//   array[N] real log_rw_pred;
// 
//   for (t in 1:N_trees) {
// 
//     array[N_years[t]] int tree_idxs = linspaced_int_array(N_years[t], tree_start_idxs[t], tree_end_idxs[t]);
//     array[N_years[t]] real years_tree = years[tree_idxs];
//     array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
// 
//     int stand_idx = stand_idxs[t];
//     int species_idx = species_idxs[t];
// 
//     vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
//     vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
//     vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];
// 
//     vector[N_years[t]] delta_sck_tree = delta_sck[tree_idxs];
// 
//     vector[N_years[t]] mu1 = alpha
//     + beta_gdd[species_idx] * (gdd_obs_tree - gdd0)
//     + beta_sm[species_idx] * (sm_obs_tree - sm0)
//     + beta_vpd[species_idx] * (vpd_obs_tree - vpd0)
//     + f_sh[stand_idx, all_years_idxs_tree]
//     + delta_sck_tree;
// 
//     vector[N_years[t]] log_rw = gp_pred_rng(years_tree, log_rw_obs[tree_idxs], years_tree,
//                                 mu1, gamma_sp[species_idx], rho_sp[species_idx], sigma, 1e-10);
// 
//     // rw[tree_idxs] = exp(log_rw[tree_idxs]);
//     log_rw_pred[tree_idxs] = normal_rng(log_rw, sigma);
// 
//   }
// }

generated quantities {
  
  vector[N] delta_sck = rep_vector(0,N); // latent amplitude of shock
  array[N] int sck_state; // latent state, zt = 0 or zt = 1
   array[N] real log_rw_pred;

  for (t in 1:N_trees) {
    
    array[N_years[t]] int tree_idxs
    = linspaced_int_array(N_years[t],
                          tree_start_idxs[t],
                          tree_end_idxs[t]);

    array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
    array[N_years[t]] real years_tree = years[tree_idxs];

    int stand_idx = stand_idxs[t];
    int species_idx = species_idxs[t];

    vector[N_years[t]] gdd_obs_tree = gdd_obs[tree_idxs];
    vector[N_years[t]] sm_obs_tree = sm_obs[tree_idxs];
    vector[N_years[t]] vpd_obs_tree = vpd_obs[tree_idxs];

    vector[N_years[t]] mu =  alpha
    + beta_gdd[species_idx] * (gdd_obs_tree - gdd0)
    + beta_sm[species_idx] * (sm_obs_tree - sm0)
    + beta_vpd[species_idx] * (vpd_obs_tree - vpd0)
    + f_sh[stand_idx, all_years_idxs_tree];
    
    // mixture weight for shock
    real mw_shock = phi_sck[stand_idx]*omega_conc_sck + (1-phi_sck[stand_idx])*omega_nonconc_sck;
    
    for(y in 1:N_years[t]){
      
      // mixture components
      real log_pshock = log(mw_shock) + normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[y] + f[tree_idxs[y]], sqrt(outer_tau2_aux[tree_idxs[y]] + sigma^2));
      real log_pshock_plus_pnonshock = log_mix(mw_shock, 
        normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[y] + f[tree_idxs[y]], sqrt(outer_tau2_aux[tree_idxs[y]] + sigma^2)),
        normal_lpdf(log_rw_obs[tree_idxs[y]] | mu[y] + f[tree_idxs[y]], sigma));
      
      // probability of shock state
      real lambda_shock = exp(log_pshock - log_pshock_plus_pnonshock);
      
      sck_state[tree_idxs[y]] = bernoulli_rng(lambda_shock); // or something like categorical_rng(lambda_shock);?

      
      if(sck_state[tree_idxs[y]] == 1) {
        // we can reconstruct shock posterior using the normal-normal conjugancy
        real residual = log_rw_obs[tree_idxs[y]] - mu[y] - f[tree_idxs[y]];
        real conjugate_mean = (outer_tau2_aux[tree_idxs[y]] / (outer_tau2_aux[tree_idxs[y]] + sigma^2)) * residual;
        real conjugate_sd   = sqrt((outer_tau2_aux[tree_idxs[y]] * sigma^2) / (outer_tau2_aux[tree_idxs[y]] + sigma^2));
        delta_sck[tree_idxs[y]] = normal_rng(conjugate_mean, conjugate_sd);
        log_rw_pred[tree_idxs[y]] = normal_rng(mu[y] + f[tree_idxs[y]] + delta_sck[tree_idxs[y]], sigma);
      }else{
        log_rw_pred[tree_idxs[y]] = normal_rng(mu[y] + f[tree_idxs[y]], sigma);
      }
      
    }
  }
}
