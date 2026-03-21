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
  
  // see truncated-random-number-generation on Stan's user guide
  real normal_ub_rng(real mu, real sigma, real ub) {
    real p_ub = normal_cdf(ub | mu, sigma);
    // deal with the limit case of the truncated normal
    // to avoid error: "Exception: uniform_rng: Upper bound parameter is 0"
    if (p_ub == 0)
      return ub;
    real u = uniform_rng(0, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
  
  // create a function for reduce_sum
  real loglikelihood_partial_sum(array[] int stand_ids_slice,
                              int start, int end,
                              array[] int N_stand_trees,
                              array[] int N_stand_years,
                              array[] int stand_tree_idxs,
                              array[] int stand_trees_start_idxs,
                              array[] int stand_trees_end_idxs,
                              array[] int tree_start_idxs,
                              array[] int tree_end_idxs,
                              int N_all_years,
                              array[] int N_years,
                              array[] int all_years_idxs,
                              array[] int stand_start_years_idxs,
                              array[] int stand_idxs,
                              array[] int species_idxs,
                              array[] int stand_species_idxs,
                              vector rw_obs,
                              vector gdd_obs,
                              vector pre_obs,
                              vector vpd_obs,
                              real gdd0,
                              real pre0,
                              real vpd0,
                              vector alpha,
                              vector beta_gdd, 
                              vector beta_pre,
                              vector beta_vpd,
                              vector kappa_sh,
                              array[] vector f_sh,
                              vector f_tilde,
                              vector f_ind_tilde,
                              array[] matrix L_cov,
                              matrix L_cov_ind,
                              real epsilon,
                              real sigma,
                              vector tau_sck,
                              vector omega_conc_sck,
                              vector omega_nonconc_sck,
                              vector phi_sck) {

    real lp = 0;
    
    profile("slice_loop") {
      for (i in 1:(end-start+1)) {
  
        int s = stand_ids_slice[i];
        
        array[N_all_years] int stand_clim_idxs = linspaced_int_array(N_all_years, 
          1+(s-1)*N_all_years, s*N_all_years);
        
        vector[N_stand_years[s]] log_p0 = rep_vector(0, N_stand_years[s]);
        vector[N_stand_years[s]] log_p1 = rep_vector(0, N_stand_years[s]);
        
        profile("trees_loop") {
          for(t in stand_tree_idxs[stand_trees_start_idxs[s]:stand_trees_end_idxs[s]]){
            
            int sp = species_idxs[t];
            int stsp = stand_species_idxs[t];
            
            int tree_start = tree_start_idxs[t];
            int tree_end  = tree_end_idxs[t];
          
            vector[N_years[t]] f;
            profile("compute_f") {
              // f = block(L_cov[sp], 1, 1, N_years[t], N_years[t]) * f_tilde[tree_start:tree_end];
              f = L_cov[sp][1:N_years[t], 1:N_years[t]] * f_tilde[tree_start:tree_end]; // 7% faster
            }
            
            vector[N_years[t]] f_ind;
            profile("compute_f_ind") {
              f_ind = L_cov_ind[1:N_years[t], 1:N_years[t]] * f_ind_tilde[tree_start:tree_end]; 
            }
            
            array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_start:tree_end];
            array[N_years[t]] int tree_clim_idxs = stand_clim_idxs[all_years_idxs_tree];
            
            vector[N_years[t]] mu;
            profile("compute_mu") {
              mu = alpha[sp]
              + beta_gdd[sp] * (gdd_obs[tree_clim_idxs] - gdd0)
              + beta_pre[sp] * (pre_obs[tree_clim_idxs] - pre0)
              + beta_vpd[sp] * (vpd_obs[tree_clim_idxs] - vpd0)
              + kappa_sh[sp] * f_sh[s,all_years_idxs_tree];
            }
            
            profile("store_trees") {
              for(y in 1:N_years[t]) {
                
                int idx = tree_start + y - 1;
                int ys =  all_years_idxs[idx]-stand_start_years_idxs[s]+1;
                
                // real log_rw;
                real mu_f = mu[y] + f[y] + f_ind[y];
                
                if(rw_obs[idx] >= epsilon){
                  real log_rw = log(rw_obs[idx]);
                  // log_p0[ys] += normal_lpdf(log_rw | mu_f, sigma);
                  log_p0[ys] += log_mix(omega_nonconc_sck[stsp],
                                  normal_lpdf(log_rw | mu_f, sqrt(tau_sck[sp]^2 + sigma^2)),
                                  normal_lpdf(log_rw | mu_f, sigma));
                  log_p1[ys] += log_mix(omega_conc_sck[stsp],
                                  normal_lpdf(log_rw | mu_f, sqrt(tau_sck[sp]^2 + sigma^2)),
                                  normal_lpdf(log_rw | mu_f, sigma));
                }else{
                  // log_p0[ys] += normal_lcdf(log(epsilon) | mu_f, sigma);
                  log_p0[ys] += log_mix(omega_nonconc_sck[stsp],
                                  normal_lcdf(log(epsilon) | mu_f, sqrt(tau_sck[sp]^2 + sigma^2)),
                                  normal_lcdf(log(epsilon) | mu_f, sigma));
                  log_p1[ys] += log_mix(omega_conc_sck[stsp],
                                  normal_lcdf(log(epsilon) | mu_f, sqrt(tau_sck[sp]^2 + sigma^2)),
                                  normal_lcdf(log(epsilon) | mu_f, sigma));
                }
              }
            }
          }
        }
        
        profile("compute_logmix") {
          for(y in 1:N_stand_years[s]) {
            lp += log_mix(phi_sck[s], log_p1[y], log_p0[y]);
          }
        }
        
      }
    }

    return lp;
  }
  
}


data {
  
  int<lower=1> N;           // Total number of observations
  int<lower=1> N_species;   // Number of species
  int<lower=1> N_trees;     // Number of trees
  int<lower=1> N_all_years; // Max number of years
  int<lower=1> N_stands;    // Number of stands
  int<lower=1> N_clades;    // Number of clades - for now, Gymnosperms=1 or Angiosperms=2
  
  array[N_stands] int<lower=1, upper=N> N_stand_trees; // Number of trees per stand
  array[N_stands] int<lower=1, upper=N> N_stand_years; // Max. number of observed years per stand
  array[N_stands] int<lower=1, upper=N> stand_start_years_idxs; // Indice of first year observed per stand
  
  // Indices of tree stands
  array[N_trees] int<lower=1, upper=N_stands> stand_idxs;
  
  // Indices of tree stand_species (species within a stand)
  int<lower=1> N_stand_species;
  array[N_trees] int<lower=1, upper=N_stand_species> stand_species_idxs;
  
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
  array[N_stands] int<lower=1, upper=N_trees> stand_trees_start_idxs;
  array[N_stands] int<lower=1, upper=N_trees> stand_trees_end_idxs;
  array[N_trees] int<lower=1, upper=N_trees> stand_tree_idxs;
  
  vector[N] rw_obs; // Log of Observed Ring Width Per 1 mm
  
  // Climate (at stand level)
  vector[N_stands*N_all_years] gdd_obs; // Observed gdd (all year) (x100 degC)
  vector[N_stands*N_all_years] pre_obs; // Observed winter precipitation (NDJFMA) (dm)
  vector[N_stands*N_all_years] vpd_obs; // Observed VPD (JJMJJA) (hPa)
  
  // Indices of species clades (Gymnosperms=1 or Angiosperms=2)
  array[N_species] int<lower=1, upper=N_clades> clade_idxs;
  
  // Partial centering parameters for growth shocks
  // vector<lower=0, upper=1>[N] w_sck;
  
  int<lower=1> grainsize;
  
}

transformed data {
  real gdd0 = 10;
  real pre0 = 5;
  real vpd0 = 23;
  
  real epsilon = 1e-3;
  
  // for reduce_sum 
  array[N_stands] int stand_ids;
  for (s in 1:N_stands) stand_ids[s] = s;

}

parameters {
  
   // Log ring width baseline
  vector[N_clades] mu_alpha;
  vector<lower=0>[N_clades] tau_alpha;
  vector[N_species] alpha;
  
  vector[N] f_tilde;
  vector[N] f_ind_tilde;
  
  // GDD slope (1/100degC)
  vector[N_clades] mu_gdd;
  vector<lower=0>[N_clades] tau_gdd;
  vector[N_species] beta_gdd;
  
  // SM slope (1/%)
  vector[N_clades] mu_pre;
  vector<lower=0>[N_clades] tau_pre;
  vector[N_species] beta_pre;
  
  // VPD slope (1/hPa)
  vector[N_clades] mu_vpd;
  vector<lower=0>[N_clades] tau_vpd;
  vector[N_species] beta_vpd;
  
  // Short-term proportional growth functional behavior
  array[N_stands] vector[N_all_years] f_tilde_sh; // Non-centered functional behavior
  real<lower=1> rho_sh;   // Time scale
  // real<lower=0> gamma_sh; // Marginal variation - now fixed to 1! (and scaled by kappa)
  
  // Individual short-term proportional growth scale 
  real<lower=1> rho_ind;
  real<lower=0> gamma_ind;
  
  // Lifetime proportional growth scale (here I implement the hard contraint on both clade and species parameters?)
  vector<lower=rho_ind>[N_clades] mu_rho;
  vector<lower=0>[N_clades] tau_rho;
  vector<lower=rho_ind>[N_species] rho_sp; 
  
  // Lifetime proportional growth variation
  vector<lower=0>[N_clades] mu_gamma;
  vector<lower=0>[N_clades] tau_gamma;
  vector<lower=0>[N_species] gamma_sp; 
  
  // Short-term proportional growth species scaling
  vector<lower=0>[N_clades] mu_kappa;
  vector<lower=0>[N_clades] tau_kappa;
  vector<lower=0>[N_species] kappa_sh; 
  
  // Growth shocks
  vector<lower=0>[N_clades] mu_tau_sck;
  vector<lower=0>[N_clades] tau_tau_sck;
  vector<lower=0>[N_species] tau_sck; // Outer yearly log variation scale (the shocks!)
  
  // Probability of stand-level shock 
  real<lower=0, upper=1> phi_sck0; // probability
  real<lower=0> tau_phi_sck; // log-odds
  vector[N_stands] alpha_phi_sck; // log-odds
  
  // Probability of tree-level shock given stand in shock (concordant shock)
  real<lower=0, upper=1> omega_conc_sck0; // probability
  real<lower=0> tau_omega_conc_sck; // log-odds
  // vector[N_stand_species] alpha_tilde_omega_conc_sck; // log-odds
  vector[N_stand_species] alpha_omega_conc_sck; // log-odds
  
  // Probability of tree-level shock given stand NOT in shock (nonconcordant shock)
  // We want the upper bound to be omega_conc_shock...
  // It's a bound that varies! 
  // For now we decide that there is no particular bound (the prior should be enough!)
  // real<lower=0, upper=1> omega_nonconc_sck0; 
  // real<lower=0> tau_omega_nonconc_sck; // log-odds
  // vector[N_stand_species] alpha_tilde_omega_nonconc_sck; // log-odds
  real mu_logdelta_omega_nonconc_sck;
  real<lower=0> tau_logdelta_omega_nonconc_sck;
  // vector[N_stand_species] delta_tilde_omega_nonconc_sck;
  vector[N_stand_species] logdelta_omega_nonconc_sck;
  
  // Proportional measurement error
  real<lower=0> sigma; 
}

transformed parameters {
  // array[N_species] real kappa_sh = append_array({1}, kappa_sh_free);
  
  array[N_stands] vector[N_all_years] f_sh;
  matrix[N_all_years, N_all_years] L_cov_sh;
  profile("L_cov_fsh") {
    {
      matrix[N_all_years, N_all_years] cov
      =   gp_exp_quad_cov(all_years, 1, rho_sh)
      + diag_matrix(rep_vector(1e-10, N_all_years));
      L_cov_sh = cholesky_decompose(cov);
      
      for (s in 1:N_stands) {
        f_sh[s] = L_cov_sh * f_tilde_sh[s];
      }
    }
  }
  
  array[N_species] matrix[N_all_years, N_all_years] L_cov;
  profile("f_lg") {
      for (sp in 1:N_species) {
        matrix[N_all_years, N_all_years] cov
        =  gp_exp_quad_cov(all_years, gamma_sp[sp], rho_sp[sp])
        + diag_matrix(rep_vector(1e-10, N_all_years));
        
        L_cov[sp] = cholesky_decompose(cov);
      }
  }
  
  matrix[N_all_years, N_all_years] L_cov_ind;
  matrix[N_all_years, N_all_years] cov_ind
        =  gp_exp_quad_cov(all_years, gamma_ind, rho_ind)
        + diag_matrix(rep_vector(1e-10, N_all_years));
  L_cov_ind = cholesky_decompose(cov_ind);
  
  real mu_phi_sck = logit(phi_sck0); // log-odds
  // vector[N_stands] alpha_phi_sck = mu_phi_sck + tau_phi_sck*alpha_tilde_phi_sck; // log-odds
  vector<lower=0, upper=1>[N_stands] phi_sck = inv_logit(alpha_phi_sck); // probabilities
  
  real mu_omega_conc_sck = logit(omega_conc_sck0); // log-odds
  // vector[N_stand_species] alpha_omega_conc_sck = mu_omega_conc_sck + tau_omega_conc_sck*alpha_tilde_omega_conc_sck; // log-odds
  vector<lower=0, upper=1>[N_stand_species] omega_conc_sck = inv_logit(alpha_omega_conc_sck); // probabilities

  // real mu_omega_nonconc_sck = logit(omega_nonconc_sck0); // log-odds
  // vector<upper=alpha_omega_conc_sck>[N_stand_species] alpha_omega_nonconc_sck = mu_omega_nonconc_sck + tau_omega_nonconc_sck*alpha_tilde_omega_nonconc_sck; // log-odds
  // vector<lower=0, upper=1>[N_stand_species] omega_nonconc_sck = inv_logit(alpha_omega_nonconc_sck); // probabilities
  
  // positive shift
  // vector[N_stand_species] delta_omega_nonconc_sck = exp(mu_delta_omega_nonconc_sck 
  //   + tau_delta_omega_nonconc_sck * delta_tilde_omega_nonconc_sck);
  vector[N_stand_species] delta_omega_nonconc_sck = exp(logdelta_omega_nonconc_sck);
  vector[N_stand_species] alpha_omega_nonconc_sck = alpha_omega_conc_sck-delta_omega_nonconc_sck;
  vector<lower=0, upper=1>[N_stand_species] omega_nonconc_sck = inv_logit(alpha_omega_nonconc_sck); // probabilities
}

model {
  
  mu_alpha ~ normal(0, log(5)/2.32); // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  tau_alpha ~ normal(0, log(5^0.5)/2.32); // before  normal(0, log(5^0.25)/2.32)
  
  mu_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  mu_pre ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_pre <~ log(1.8)
  mu_vpd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_vpd <~ log(1.8)
  tau_gdd ~ normal(0, log(1.8^0.25) / 2.57); // variation of the order of 25%?
  tau_pre ~ normal(0, log(1.8^0.25) / 2.57); // variation of the order of 25%?
  tau_vpd ~ normal(0, log(1.8^0.25) / 2.57); // variation of the order of 25%?
  
  mu_rho ~ lognormal(3.55, 0.24); // 20 <~ rho <~ 60 // MODIFIED THIS
  tau_rho ~ normal(0, 20 / 2.57); // max. variation of 20 years across species? // MODIFIED THIS
  
  mu_gamma ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)
  tau_gamma ~ normal(0, 0.23 / 2.57); // max. variation of the order of 10% for max. gamma = log(10)? 
  
  mu_kappa ~ lognormal(0, 0.41 / 2.32); // 2/3 <~ kappa_sh <~ 3/2
  tau_kappa ~ normal(0, 0.2 / 2.57); // variation of the order of 10% for max. kappa = 2?
  
  mu_tau_sck ~ normal(0, log(20) / 2.57); // 2/3 <~ kappa_sh <~ 3/2
  tau_tau_sck ~ normal(0, log(2) / 2.57); // variation of the order of 10% for max. kappa = 2?
  
  for (sp in 1:N_species) {
    alpha[sp] ~ normal(mu_alpha[clade_idxs[sp]], tau_alpha[clade_idxs[sp]]);
    
    beta_gdd[sp] ~ normal(mu_gdd[clade_idxs[sp]], tau_gdd[clade_idxs[sp]]);
    beta_pre[sp] ~ normal(mu_pre[clade_idxs[sp]], tau_pre[clade_idxs[sp]]);
    beta_vpd[sp] ~ normal(mu_vpd[clade_idxs[sp]], tau_vpd[clade_idxs[sp]]);

    rho_sp[sp] ~ normal(mu_rho[clade_idxs[sp]], tau_rho[clade_idxs[sp]]);
    gamma_sp[sp] ~ normal(mu_gamma[clade_idxs[sp]] , tau_gamma[clade_idxs[sp]]);

    kappa_sh[sp] ~ normal(mu_kappa[clade_idxs[sp]], tau_kappa[clade_idxs[sp]]);

    tau_sck[sp] ~ normal(mu_tau_sck[clade_idxs[sp]], tau_tau_sck[clade_idxs[sp]]);
  }

  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26);       // 3 <~ rho_sh <~ 10
  //gamma_sh ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)
  
  rho_ind ~ lognormal(0.80, 0.40);   // 1 <~ rho_sh <~ 5
  gamma_ind ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)
  
  phi_sck0 ~ beta(2, 20); 
  tau_phi_sck ~ normal(0, 1); // MODIFIED THIS
  alpha_phi_sck ~ normal(mu_phi_sck, tau_phi_sck);
  
  omega_conc_sck0 ~ beta(30, 20); // MODIFIED THIS... concordance is something 50-75% trees that shock on average...
  tau_omega_conc_sck ~ normal(0, 1); // MODIFIED THIS
  // alpha_tilde_omega_conc_sck ~ normal(0, 1);
  alpha_omega_conc_sck ~ normal(mu_omega_conc_sck, tau_omega_conc_sck); 
  
  // omega_nonconc_sck0 ~ beta(1, 72); // 0 <~ omega_conc_sck <~ 0.05 (should be rare, but... who knows?)
  // tau_omega_nonconc_sck ~ normal(0, 1/2.57); 
  // alpha_tilde_omega_nonconc_sck ~ normal(0, 1);
  
  mu_logdelta_omega_nonconc_sck ~ normal(log(5), log(3)/2.57); // before:  normal(log(8), log(2)/2.57)
  tau_logdelta_omega_nonconc_sck ~ normal(0, 1/2.57); // before: normal(0, 0.3/2.57)
  // delta_tilde_omega_nonconc_sck ~ normal(0, 1);
  logdelta_omega_nonconc_sck ~ normal(mu_logdelta_omega_nonconc_sck, tau_logdelta_omega_nonconc_sck);
  
  sigma ~ normal(0, log(1.1) / 2.57);   // 0 <~ sigma <~ +log(1.1)
  
  profile("prepare") {
    f_tilde ~ normal(0,1);
    f_ind_tilde  ~ normal(0,1);
  }
  
  // Observational model with marginalized shocks
   profile("likelihood") {
     target += reduce_sum(
      loglikelihood_partial_sum,
      stand_ids,
      grainsize, // grain size
      N_stand_trees,
      N_stand_years,
      stand_tree_idxs,
      stand_trees_start_idxs,
      stand_trees_end_idxs,
      tree_start_idxs,
      tree_end_idxs,
      N_all_years,
      N_years,
      all_years_idxs,
      stand_start_years_idxs,
      stand_idxs,
      species_idxs,
      stand_species_idxs,
      rw_obs,
      gdd_obs,
      pre_obs,
      vpd_obs,
      gdd0,
      pre0,
      vpd0,
      alpha,
      beta_gdd, 
      beta_pre,
      beta_vpd,
      kappa_sh,
      f_sh,
      f_tilde,
      f_ind_tilde,
      L_cov,
      L_cov_ind,
      epsilon,
      sigma,
      tau_sck,
      omega_conc_sck,
      omega_nonconc_sck,
      phi_sck);
   }
  
}
// 
// generated quantities {
// 
//   vector[N] delta_sck = rep_vector(0,N); // latent amplitude of shock
//   array[N] int sck_state; // latent state, zt = 0 or zt = 1
//   array[N] real log_rw_pred;
// 
//   for (t in 1:N_trees) {
// 
//     int stand_idx = stand_idxs[t];
//     int species_idx = species_idxs[t];
//     int stand_species_idx = stand_species_idxs[t];
// 
//     array[N_all_years] int stand_clim_idxs = linspaced_int_array(N_all_years,
//           1+(stand_idx-1)*N_all_years, stand_idx*N_all_years);
// 
//     int tree_start = tree_start_idxs[t];
//     int tree_end  = tree_end_idxs[t];
// 
//     vector[N_years[t]] f;
// 
//     // f = block(L_cov[species_idx], 1, 1, N_years[t], N_years[t]) * f_tilde[tree_start:tree_end];
//     f = L_cov[species_idx][1:N_years[t], 1:N_years[t]] * f_tilde[tree_start:tree_end]; // 7% faster
// 
//     array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_start:tree_end];
//     array[N_years[t]] int tree_clim_idxs = stand_clim_idxs[all_years_idxs_tree];
// 
// 
//     vector[N_years[t]] mu;
//     mu = alpha
//     + beta_gdd[species_idx] * (gdd_obs[tree_clim_idxs] - gdd0)
//     + beta_pre[species_idx] * (pre_obs[tree_clim_idxs] - pre0)
//     + beta_vpd[species_idx] * (vpd_obs[tree_clim_idxs] - vpd0)
//     + kappa_sh[species_idx] * f_sh[stand_idx, all_years_idxs_tree];
// 
//     // mixture weight for shock
//     real mw_shock = phi_sck[stand_idx]*omega_conc_sck[stand_species_idx] + (1-phi_sck[stand_idx])*omega_nonconc_sck[stand_species_idx];
//     // real mw_shock = phi_sck[stand_idx]*omega_conc_sck[stand_species_idx];
//     real log_pshock;
//     real log_pshock_plus_pnonshock;
// 
//     for(y in 1:N_years[t]){
// 
//       int idx = tree_start + y - 1;
//       real mu_f = mu[y] + f[y];
// 
//       if(rw_obs[idx] >= epsilon){
//         real log_rw = log(rw_obs[idx]);
// 
//         log_pshock = log(mw_shock) + normal_lpdf(log_rw | mu_f,
//           sqrt(tau_sck[species_idx]^2 + sigma^2));
//         log_pshock_plus_pnonshock = log_mix(mw_shock,
//           normal_lpdf(log_rw | mu_f, sqrt(tau_sck[species_idx]^2 + sigma^2)),
//           normal_lpdf(log_rw | mu_f, sigma));
//       }else{
//         log_pshock = log(mw_shock) + normal_lcdf(log(epsilon) | mu_f,
//           sqrt(tau_sck[species_idx]^2 + sigma^2));
//         log_pshock_plus_pnonshock = log_mix(mw_shock,
//           normal_lcdf(log(epsilon)| mu_f, sqrt(tau_sck[species_idx]^2 + sigma^2)),
//           normal_lcdf(log(epsilon) | mu_f, sigma));
//       }
// 
//       // probability of shock state
//       real lambda_shock = exp(log_pshock - log_pshock_plus_pnonshock);
// 
//       sck_state[idx] = bernoulli_rng(lambda_shock); // or something like categorical_rng(lambda_shock);?
// 
//       if(sck_state[idx] == 0){
//         log_rw_pred[idx] = normal_rng(mu_f, sigma);
//       }else if(rw_obs[idx] >= epsilon){
//         real log_rw = log(rw_obs[idx]);
//         // we can reconstruct shock posterior using the normal-normal conjugancy
//         real residual = log_rw - mu_f;
//         real conjugate_mean = (tau_sck[species_idx]^2 / (tau_sck[species_idx]^2 + sigma^2)) * residual;
//         real conjugate_sd   = sqrt((tau_sck[species_idx]^2 * sigma^2) / (tau_sck[species_idx]^2 + sigma^2));
//         delta_sck[idx] = normal_rng(conjugate_mean, conjugate_sd);
//         log_rw_pred[idx] = normal_rng(mu_f + delta_sck[idx], sigma);
//       }else{
//         // sample from a truncated normal distribution? between -inf and log(epsilon)
//         real log_y = normal_ub_rng(mu_f, sqrt(tau_sck[species_idx]^2 + sigma^2), log(epsilon));
//         real residual = log_y - mu_f;
//         real conjugate_mean = (tau_sck[species_idx]^2 / (tau_sck[species_idx]^2 + sigma^2)) * residual;
//         real conjugate_sd   = sqrt((tau_sck[species_idx]^2 * sigma^2) / (tau_sck[species_idx]^2 + sigma^2));
//         delta_sck[idx] = normal_rng(conjugate_mean, conjugate_sd);
//         log_rw_pred[idx] = normal_rng(mu_f + delta_sck[idx], sigma);
//       }
// 
//     }
//   }
// }
