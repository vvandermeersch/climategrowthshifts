functions {
  
  // GP helpers
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
      matrix[N1, N1] K = gp_exp_quad_cov(x1, alpha1, rho1)
      + gp_exp_quad_cov(x1, alpha2, rho2)
      + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);
      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1 - mu);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = gp_exp_quad_cov(x1, x2, alpha1, rho1)
                              + gp_exp_quad_cov(x1, x2, alpha2, rho2);
      vector[N2] f2_mu = mu + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 = gp_exp_quad_cov(x2, alpha1, rho1)
                              + gp_exp_quad_cov(x2, alpha2, rho2)
                              - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
  
  // Function to sample from a truncated normal (upper bound ub)
  real normal_ub_rng(real mu, real sigma, real ub) {
    real p_ub = normal_cdf(ub | mu, sigma);
    if (p_ub == 0)
      return ub;
    real u = uniform_rng(0, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
  
  // HS basis functions, from https://users.aalto.fi/~ave/casestudies/Motorcycle/motorcycle.html
  vector diagSPD_EQ(real alpha, real rho, real L, int M) {
    return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }

  matrix PHI(int N, int M, real L, vector x) {
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }
  
  // For parallelization with reduce_sum
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
                              vector alpha_stand,
                              vector beta_gdd,
                              vector beta_pre,
                              vector beta_vpd,
                              array[] vector delta_clim,
                              vector kappa_clim,
                              matrix f_tilde_ind,
                              matrix PHI_sp,
                              array[] vector sqrt_spd_ind,
                              real epsilon,
                              real sigma,
                              vector omega_conc_sck,
                              vector omega_shutdown,
                              vector logit_phi_sck,
                              real beta_phi_vpd,
                              real beta_phi_pre,
                              array[] vector thetas_idio,
                              real sigma_idio,
                              vector sigma_conc,
                              vector sigma_idio_conc){
    
    real lp = 0;
    profile("lkhd_in") {
      
      for (i in 1:(end-start+1)) {
    
        int st = stand_ids_slice[i];
        
        array[N_all_years] int stand_clim_idxs = linspaced_int_array(N_all_years,
          1+(st-1)*N_all_years, st*N_all_years);
        
        vector[N_stand_years[st]] lpd_nonconc = rep_vector(0, N_stand_years[st]);
        vector[N_stand_years[st]] lpd_conc = rep_vector(0, N_stand_years[st]);
        
        for (t in stand_tree_idxs[stand_trees_start_idxs[st]:stand_trees_end_idxs[st]]) {
          
          int sp = species_idxs[t];
          int stsp = stand_species_idxs[t];
          
          int tree_start = tree_start_idxs[t];
          int tree_end = tree_end_idxs[t];
          
          array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_start:tree_end];
          array[N_years[t]] int tree_clim_idxs = stand_clim_idxs[all_years_idxs_tree];
  
          vector[N_years[t]] f_ind = PHI_sp[all_years_idxs_tree, ] * (sqrt_spd_ind[sp] .* f_tilde_ind[, t]);
          
          vector[N_years[t]] mu = alpha[sp] + alpha_stand[st]
            + beta_gdd[sp] * (gdd_obs[tree_clim_idxs] - gdd0)
            + beta_pre[sp] * (pre_obs[tree_clim_idxs] - pre0)
            + beta_vpd[sp] * (vpd_obs[tree_clim_idxs] - vpd0)
            // + kappa_sh[sp] * f_sh[st, all_years_idxs_tree]
            + kappa_clim[sp] * delta_clim[st, all_years_idxs_tree];
          
          for (y in 1:N_years[t]) {
            
            int idx = tree_start + y - 1;
            int ys = all_years_idxs_tree[y] - stand_start_years_idxs[st] + 1;
            
            real mu_f = mu[y] + f_ind[y];
            
            if (rw_obs[idx] >= epsilon){
              
              real log_rw = log(rw_obs[idx]);
              
              vector[4] lpds = [
                normal_lpdf(log_rw | mu_f, sigma), // no concordant shock, no idio. shutdown
                normal_lpdf(log_rw | mu_f, sigma_idio), // no concordant shock, idio. shock
                
                normal_lpdf(log_rw | mu_f, sigma_conc[sp]), // concordant depressed growth, no idio. shutdown
                normal_lpdf(log_rw | mu_f, sigma_idio_conc[sp]) // concordant depressed growth, idio. shock
              ]';
              
              vector[4] lambdas = [
                (1-omega_conc_sck[stsp])*thetas_idio[t][1], 
                (1-omega_conc_sck[stsp])*thetas_idio[t][2], 
                
                omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][1],
                omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][2]
              ]'; // implicit:  shutdown cannot happen (log(0))
              
              lpd_nonconc[ys] += log_sum_exp(
                log([thetas_idio[t][1], thetas_idio[t][2]]') +
                [lpds[1], lpds[2]]'
              );
              
              lpd_conc[ys] += log_sum_exp(log(lambdas) + lpds);
  
            }else{
              
              vector[9] lpds = [
                normal_lcdf(log(epsilon) | mu_f, sigma), // no concordant shock, no idio. shutdown
                normal_lcdf(log(epsilon) | mu_f, sigma_idio), // no concordant shock, idio. shock
                log(1), // no concordant shock, idiosync. shutdown
                
                normal_lcdf(log(epsilon)| mu_f, sigma_conc[sp]), // concordant depressed growth, no idio. shutdown
                normal_lcdf(log(epsilon) | mu_f, sigma_idio_conc[sp]), // concordant depressed growth, idio. shock
                log(1), // concordant depressed growth, idiosync. shutdown
                
                log(1), // concordant shutdown, no idiosync. shock
                log(1), // concordant shutdown, idio. shock
                log(1) // concordant shutdown, idiosync. shutdown
              ]';
  
              vector[9] lambdas = [
                (1-omega_conc_sck[stsp])*thetas_idio[t][1], 
                (1-omega_conc_sck[stsp])*thetas_idio[t][2], 
                (1-omega_conc_sck[stsp])*thetas_idio[t][3], 
                
                omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][1],
                omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][2],
                omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][3],
                
                omega_conc_sck[stsp]*omega_shutdown[stsp]*thetas_idio[t][1],
                omega_conc_sck[stsp]*omega_shutdown[stsp]*thetas_idio[t][2],
                omega_conc_sck[stsp]*omega_shutdown[stsp]*thetas_idio[t][3]
              ]';
              
              lpd_nonconc[ys] += log_sum_exp(
                log([thetas_idio[t][1], thetas_idio[t][2], thetas_idio[t][3]]') +
                [lpds[1], lpds[2], lpds[3]]'
              );
              
              lpd_conc[ys] += log_sum_exp(log(lambdas) + lpds);
            }
          }
        }
        
        for (y in 1:N_stand_years[st]) {
          
          int ys = (st - 1) * N_all_years + stand_start_years_idxs[st] + y - 1;
          
          real phi_sck_y =  inv_logit(logit_phi_sck[st]
          + beta_phi_vpd * (vpd_obs[ys] - vpd0)
          + beta_phi_pre * (pre_obs[ys] - pre0));
          lp += log_mix(phi_sck_y, lpd_conc[y], lpd_nonconc[y]);
        }
      }
    }
    return lp;
  }
  
  
}

data {
  
  int<lower=1> N; // number of observations
  int<lower=1> N_trees; // number of trees
  int<lower=1> N_all_years; // max. number of years observed
  int<lower=1> N_stands; // number of stands
  
  array[N_stands] int<lower=1, upper=N> N_stand_trees; // number of trees within each stand
  array[N_stands] int<lower=1, upper=N> N_stand_years; // number of years observed, for each stand
  array[N_stands] int<lower=1, upper=N> stand_start_years_idxs; // first year observed, for each stand
  
  array[N_trees] int<lower=1, upper=N_stands> stand_idxs; // stand idx, for each tree
  array[N_trees] int<lower=1, upper=N> N_years; // number of years observed, for each tree
  array[N_all_years] real all_years; // 'true' years (e.g 1900-2000)
  array[N] real years; // 'true' year , for each observaton (not used)
  array[N] int<lower=1, upper=N_all_years> all_years_idxs; // year idx, for each observaton
  
  array[N_trees] int<lower=1, upper=N> tree_start_idxs; //  for each tree, start idx for ragged vector of observations
  array[N_trees] int<lower=1, upper=N> tree_end_idxs; // ... and end idx 
  
  array[N_stands] int<lower=1, upper=N_trees> stand_trees_start_idxs; // for each stand, start idx for ragged vector of trees
  array[N_stands] int<lower=1, upper=N_trees> stand_trees_end_idxs; //  ... and end idx
  array[N_trees] int<lower=1, upper=N_trees> stand_tree_idxs;
  
  vector[N] rw_obs; // the observations! ring widths in mm
  
  // climate covariates, at the stand-level
  vector[N_stands*N_all_years] gdd_obs; // GDD during entire year (x100 degC, kdegC?!)
  vector[N_stands*N_all_years] pre_obs; // Winter precipitation, NDJFMA (dm)
  vector[N_stands*N_all_years] vpd_obs; // VPD, JJA (hPa)
  
  int<lower=1> grainsize; // to tweak reduce_sum efficiency
  
  // species!
  int<lower=1> N_species;
  int<lower=1> N_stand_species;
  int<lower=1> N_clades;
  
  array[N_trees] int<lower=1, upper=N_species> species_idxs; // species idx, for each tree
  array[N_trees] int<lower=1, upper=N_stand_species> stand_species_idxs; // stand_species (species within a stand), for each tree
  
}

transformed data {
  
  // baselines
  real gdd0 = 10;
  real pre0 = 5;
  real vpd0 = 23;
  
  real epsilon = 1e-3; // measurement precision threshold (could be estimated...)
  
  // HS approximation (Riutort-Mayol et al. 2023)
  int M = 20;  // number of basis functions
  vector[N_all_years] years_c = to_vector(all_years)
                                - mean(to_vector(all_years)); 
  real S = (max(years_c) - min(years_c)) / 2.0; // half range
  real c = 2.0; // boundary factor: L = c * S
  real L = c * S; // boundary
  matrix[N_all_years, M] PHI_sp = PHI(N_all_years, M, L, years_c);  // doesn't depend on parameters!
  
  array[N_stands] int stand_ids;
  for (s in 1:N_stands) stand_ids[s] = s;
  
}

parameters {
  
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[N_species] alpha;
  
  // real mu_alpha_stand;
  real<lower=0> sigma_alpha_stand;
  vector[N_stands] alpha_stand;
  
  real mu_beta_gdd;
  real<lower=0> sigma_beta_gdd;
  vector[N_species] beta_gdd;
  
  real mu_beta_pre;
  real<lower=0> sigma_beta_pre;
  vector[N_species] beta_pre;
  
  real mu_beta_vpd;
  real<lower=0> sigma_beta_vpd;
  vector[N_species] beta_vpd;
  
  array[N_stands] vector[N_all_years] delta_clim;
  real<lower=0> tau_clim;
  // array[N_species-1] real<lower=0> kappa_clim_free;
  sum_to_zero_vector[N_species] log_kappa_clim; // new constraint!
  
  // array[N_stands] vector[N_all_years] f_tilde_sh; // short-term proportional growth functional behavior
  // real<lower=1> rho_sh; // length scale
  // real<lower=0> gamma_sh; // marginal variation
  // array[N_species-1] real<lower=0> kappa_sh_free;
  
  matrix[M, N_trees] f_tilde_ind; // tree-level func. behavior, short-term + long-term (HSGP)
  real mu_log_rho;
  real<lower=0> sigma_log_rho;
  vector<lower=1>[N_species] rho_merged; // length scale
  
  real mu_log_gamma;
  real<lower=0> sigma_log_gamma;
  vector[N_species] log_gamma_merged; // marginal variation
  
  // The (original) shocks!
  // vector<lower=0>[N_stands] mu_conc; // log variation location
  real mu_log_tau_conc;
  real<lower=0> sigma_log_tau_conc; 
  vector<lower=0>[N_species] tau_conc; // log variation scale
  
  real mu_phi;
  real<lower=0> sigma_phi;
  vector[N_stands] logit_phi_sck; // baselines!
  real beta_phi_vpd; // effect of summer VPD anomaly on shock-year log-odds
  real beta_phi_pre; // effect of winter precip anomaly on shock-year log-odds
  
  real mu_omega_conc;
  real<lower=0> sigma_omega_conc;
  vector[N_stand_species] logit_omega_conc_sck;
  
  real mu_omega_shutdown;
  real<lower=0> sigma_omega_shutdown;
  vector[N_stand_species] logit_omega_shutdown;
  
  // Idiosyncratic shocks!
  array[N_trees] simplex[3] thetas_idio;
  real<lower=0> tau_idio; // log variation scale
  
  real<lower=0> sigma; // proportional measurement error
  
}

transformed parameters {
  
  vector<lower=0, upper=1>[N_stands] phi_sck = inv_logit(logit_phi_sck);
  vector<lower=0, upper=1>[N_stand_species] omega_conc_sck = inv_logit(logit_omega_conc_sck);
  vector<lower=0, upper=1>[N_stand_species] omega_shutdown = inv_logit(logit_omega_shutdown);
  
  vector<lower=0>[N_species] gamma_merged = exp(log_gamma_merged);
  
  vector[N_species] kappa_clim = exp(log_kappa_clim);
  // array[N_species] real kappa_clim = append_array({1}, kappa_clim_free);
  // array[N_species] real kappa_sh = append_array({1}, kappa_sh_free);
  
  // Stand-level short-term GP
  // array[N_stands] vector[N_all_years] f_sh;
  // matrix[N_all_years, N_all_years] L_cov_sh;
  // profile("L_cov_fsh") {
  //   {
  //     matrix[N_all_years, N_all_years] cov
  //       = gp_exp_quad_cov(all_years, gamma_sh, rho_sh)
  //       + diag_matrix(rep_vector(1e-8, N_all_years));
  //     L_cov_sh = cholesky_decompose(cov);
  //     for (s in 1:N_stands)
  //       f_sh[s] = L_cov_sh * f_tilde_sh[s];
  //   }
  // }
  
  // Tree-level GP (HSGP)
  array[N_species] vector[M] sqrt_spd_ind;
  for (sp in 1:N_species)
    sqrt_spd_ind[sp] = diagSPD_EQ(gamma_merged[sp], rho_merged[sp], L, M);
    
  real sigma_idio = sqrt(tau_idio^2 + sigma^2);
  vector[N_species] sigma_conc = sqrt(tau_conc^2 + sigma^2);
  vector[N_species] sigma_idio_conc = sqrt(tau_idio^2 + tau_conc^2 + sigma^2);
  
}

model {
  
  mu_alpha ~ normal(0, log(10)/2.32);
  sigma_alpha ~ normal(0, log(10)/2.32);
  alpha ~ normal(mu_alpha, sigma_alpha);
  
  // mu_alpha_stand ~ normal(0, log(10)/2.32);
  sigma_alpha_stand ~ normal(0, log(10)/2.32);
  alpha_stand ~ normal(0, sigma_alpha_stand);
  
  mu_beta_gdd ~ normal(0, log(1.3) / 2.57);
  sigma_beta_gdd ~ normal(0, log(1.3) / 2.57);
  beta_gdd ~ normal(mu_beta_gdd, sigma_beta_gdd);
  
  mu_beta_pre ~ normal(0, log(1.3) / 2.57);
  sigma_beta_pre ~ normal(0, log(1.3) / 2.57);
  beta_pre ~ normal(mu_beta_pre, sigma_beta_pre);
  
  mu_beta_vpd ~ normal(0, log(1.3) / 2.57);
  sigma_beta_vpd ~ normal(0, log(1.3) / 2.57);
  beta_vpd ~ normal(mu_beta_vpd, sigma_beta_vpd);
  
  mu_log_rho ~ normal(log(15), 0.5);
  sigma_log_rho ~ normal(0, 0.5);
  rho_merged ~ lognormal(mu_log_rho, sigma_log_rho);
  
  // was gamma ~ normal(0, log(5)/2.57) before partial pooling
  mu_log_gamma ~ normal(log(0.42), 0.58); 
  sigma_log_gamma ~ normal(0, 0.5); 
  log_gamma_merged ~ normal(mu_log_gamma, sigma_log_gamma); 
  
  to_vector(f_tilde_ind) ~ normal(0, 1);
  
  mu_log_tau_conc ~ normal(log(0.8), 0.5);
  sigma_log_tau_conc ~ normal(0, 0.5);
  tau_conc ~ lognormal(mu_log_tau_conc, sigma_log_tau_conc);

  // for (s in 1:N_stands)
  //   f_tilde_sh[s] ~ normal(0, 1);
  // rho_sh ~ lognormal(3.0, 0.42); // 10 <~ rho_sh <~ 40
  // gamma_sh ~ normal(0, log(2) / 2.57); // 0 < gamma_sh < log(3)
  // kappa_sh_free ~ lognormal(0, 0.41 / 2.32); // 2/3 <~ kappa_sh <~ 3/2
  
  tau_clim ~ normal(0, log(2)/2.57); 
  for (s in 1:N_stands)
    delta_clim[s] ~ normal(0, tau_clim);
  // kappa_clim_free ~ lognormal(0, 0.41 / 2.32); // 2/3 <~ kappa_sh <~ 3/2
  log_kappa_clim ~ normal(0, 0.41 / 2.32);
  
  mu_phi ~ normal(-2.17, 0.40); // 5-20%
  mu_omega_conc ~ normal(-0.49, 0.46); // 0-60%
  mu_omega_shutdown ~ normal(-3.40, 0.61); // 1-10%
  sigma_phi ~ normal(0, 1.0); // up to 70% (already a lot)
  sigma_omega_conc ~ normal(0, 1.5); // up to ~100%
  sigma_omega_shutdown ~ normal(0, 0.95); // up to 70% (already a lot)
  
  logit_phi_sck ~ normal(mu_phi, sigma_phi);
  logit_omega_conc_sck ~ normal(mu_omega_conc, sigma_omega_conc);
  logit_omega_shutdown ~ normal(mu_omega_shutdown, sigma_omega_shutdown);
  
  beta_phi_vpd ~ normal(0, 0.3); // not well thought
  beta_phi_pre ~ normal(0, 0.3); // 
  
  // Idiosyncratic shocks
  vector[3] thetas_baseline = [100, 20, 1]';
  real omega_thetas = 8;
  vector[3] ones = rep_vector(1, 3);
  vector[3] alphas = thetas_baseline / omega_thetas + ones;
  
  for (t in 1:N_trees) {
    thetas_idio[t] ~ dirichlet(alphas);
  } 
  
  tau_idio ~ normal(0, log(5) / 2.57); 
  
  sigma ~ normal(0, log(1.1) / 2.57); // we expected 10% max.
  
  profile("lkhd_out") {
    target += reduce_sum(
      loglikelihood_partial_sum,
      stand_ids,
      grainsize,
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
      alpha_stand,
      beta_gdd,
      beta_pre,
      beta_vpd,
      delta_clim,
      kappa_clim,
      f_tilde_ind,
      PHI_sp,
      sqrt_spd_ind,
      epsilon,
      sigma,
      omega_conc_sck,
      omega_shutdown,
      logit_phi_sck,
      beta_phi_vpd,
      beta_phi_pre,
      thetas_idio,
      sigma_idio,
      sigma_conc,
      sigma_idio_conc);
  }

}

