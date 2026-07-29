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
  vector[N_species] log_kappa_clim; // new constraint! // was sum_to_zero_vector, modified for standalone GQ
  
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
  array[N_trees] vector[3] thetas_idio; // was simplex, modified for standalone GQ
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


generated quantities {

  array[N_stands, N_all_years] int sck_year = rep_array(-1, N_stands, N_all_years); // latent concordant states, -1 are unobserved years

  array[N] real log_rw_pred;
  vector[N] f_ind;
  vector[N] mu_f; 
  vector[N] delta_sck = rep_vector(0,N); // latent amplitude of shock
  vector[N] shutdown = rep_vector(0,N);  // shutdown state
  
  vector[N] tree_conc_state = rep_vector(0,N); 
  vector[N] tree_idio_state = rep_vector(0,N); 
  
  // reconstruct the (latent) concordant state at the stand level
  for (st in 1:N_stands) {
    
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
      
      f_ind[tree_start:tree_end] = PHI_sp[all_years_idxs_tree, ] * (sqrt_spd_ind[sp] .* f_tilde_ind[, t]);
      
      mu_f[tree_start:tree_end] = alpha[sp] + alpha_stand[st]
        + beta_gdd[sp] * (gdd_obs[tree_clim_idxs] - gdd0)
        + beta_pre[sp] * (pre_obs[tree_clim_idxs] - pre0)
        + beta_vpd[sp] * (vpd_obs[tree_clim_idxs] - vpd0)
        + kappa_clim[sp] * delta_clim[st, all_years_idxs_tree]
        + f_ind[tree_start:tree_end];
      
      for (y in 1:N_years[t]) {
        
        int idx = tree_start + y - 1;
        int ys = all_years_idxs_tree[y] - stand_start_years_idxs[st] + 1;
        
        if (rw_obs[idx] >= epsilon) {
          
          real log_rw = log(rw_obs[idx]);
          
          vector[4] lpds = [
            normal_lpdf(log_rw | mu_f[idx], sigma),
            normal_lpdf(log_rw | mu_f[idx], sigma_idio),
            normal_lpdf(log_rw | mu_f[idx], sigma_conc[sp]),
            normal_lpdf(log_rw | mu_f[idx], sigma_idio_conc[sp])
          ]';
          
          vector[4] lambdas = [
            (1-omega_conc_sck[stsp])*thetas_idio[t][1], 
            (1-omega_conc_sck[stsp])*thetas_idio[t][2], 
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][1],
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][2]
          ]';
          
          lpd_nonconc[ys] += log_sum_exp(
            log([thetas_idio[t][1], thetas_idio[t][2]]') +
            [lpds[1], lpds[2]]'
          );
          
          lpd_conc[ys] += log_sum_exp(log(lambdas) + lpds);
          
        } else {
          
          vector[9] lpds = [
            normal_lcdf(log(epsilon) | mu_f[idx], sigma),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_idio),
            log(1),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_conc[sp]),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_idio_conc[sp]),
            log(1),
            log(1),
            log(1),
            log(1)
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
      int ay = stand_start_years_idxs[st] + y - 1;
      
      real phi = inv_logit(logit_phi_sck[st]
        + beta_phi_vpd * (vpd_obs[ys] - vpd0)
        + beta_phi_pre * (pre_obs[ys] - pre0));
        
      real p_sck = inv_logit(log(phi) - log1m(phi)
        + lpd_conc[y] - lpd_nonconc[y]);
      
      sck_year[st, ay] = bernoulli_rng(p_sck);
    }
  }
  
  for (t in 1:N_trees) {
    
    int sp = species_idxs[t];
    int st = stand_idxs[t];
    int stsp = stand_species_idxs[t];
    
    int tree_start = tree_start_idxs[t];
    int tree_end = tree_end_idxs[t];
    
    array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_start:tree_end];
    
    for (y in 1:N_years[t]) {
      
      int idx = tree_start + y - 1;
      int sck_cat;
      
      // no concordant event
      if (sck_year[st, all_years_idxs_tree[y]] == 0) {
        
        if (rw_obs[idx] >= epsilon) {
          
          real log_rw = log(rw_obs[idx]);
          
          vector[2] lpds = [
            normal_lpdf(log_rw | mu_f[idx], sigma),
            normal_lpdf(log_rw | mu_f[idx], sigma_idio)
          ]';
          
          vector[2] lambdas = [thetas_idio[t][1], thetas_idio[t][2]]';
          
          vector[2] jlps = log(lambdas) + lpds;
          
          real jlp_all = log_sum_exp(jlps);
          
          vector[2] jp = exp(jlps - jlp_all);
          jp /= sum(jp); 
          
          sck_cat = categorical_rng(jp);
          
        } else {
          
          vector[3] lpds = [
            normal_lcdf(log(epsilon) | mu_f[idx], sigma),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_idio),
            log(1)
          ]';
          
          vector[3] lambdas = thetas_idio[t];
          
          vector[3] jlps = log(lambdas) + lpds;
          
          real jlp_all = log_sum_exp(jlps);
          
          vector[3] jp = exp(jlps - jlp_all);
          jp /= sum(jp); 
          
          sck_cat = categorical_rng(jp);
        }
        
        if (sck_cat == 1) { // no idiosync. shock
          
          if (rw_obs[idx] >= epsilon) {
            log_rw_pred[idx] = normal_rng(mu_f[idx], sigma);
          } else { // sub-threshold observation
            log_rw_pred[idx] = normal_ub_rng(mu_f[idx], sigma, log(epsilon));
          }
          
        } else if (sck_cat == 2) { // idiosync. depressed growth
          
          real log_rw;
          real residual;
          
          if (rw_obs[idx] >= epsilon) {
            
            log_rw = log(rw_obs[idx]);
            residual = log_rw - mu_f[idx];
            
          } else { // sub-threshold observation
            
            // sample from a truncated normal distribution? between -inf and log(epsilon)
            log_rw = normal_ub_rng(mu_f[idx], sqrt(tau_idio^2 + sigma^2), log(epsilon));
            residual = log_rw - mu_f[idx];
            
          }
          
          // we can reconstruct shock posterior using the normal-normal conjugancy
          real conjugate_mean = (tau_idio^2 / (tau_idio^2 + sigma^2)) * residual;
          real conjugate_sd   = sqrt((tau_idio^2 * sigma^2) / (tau_idio^2 + sigma^2));
          
          delta_sck[idx] = normal_rng(conjugate_mean, conjugate_sd);
          tree_idio_state[idx] = 1;
          
          log_rw_pred[idx] = normal_rng(mu_f[idx] + delta_sck[idx], sigma);
          
        } else { // sck_cat == 3, idiosync. shutdown
          
          log_rw_pred[idx] = not_a_number();
          tree_idio_state[idx] = 1;
          shutdown[idx] = 1;
          
        }
        
      } else {
        // concordant event!
        
        if (rw_obs[idx] >= epsilon) {
          
          real log_rw = log(rw_obs[idx]);
          
          vector[4] lpds = [
            normal_lpdf(log_rw | mu_f[idx], sigma),
            normal_lpdf(log_rw | mu_f[idx], sigma_idio),
            normal_lpdf(log_rw | mu_f[idx], sigma_conc[sp]),
            normal_lpdf(log_rw | mu_f[idx], sigma_idio_conc[sp])
          ]';
          
          vector[4] lambdas = [
            (1-omega_conc_sck[stsp])*thetas_idio[t][1], 
            (1-omega_conc_sck[stsp])*thetas_idio[t][2], 
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][1],
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][2]
          ]';
          
          vector[4] jlps = log(lambdas) + lpds;
          
          real jlp_all = log_sum_exp(jlps);
          
          vector[4] jp = exp(jlps - jlp_all);
          jp /= sum(jp); // floating stability
          
          sck_cat = categorical_rng(jp);
          
        } else {
          
          vector[9] lpds = [
            normal_lcdf(log(epsilon) | mu_f[idx], sigma),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_idio),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_conc[sp]),
            normal_lcdf(log(epsilon) | mu_f[idx], sigma_idio_conc[sp]),
            log(1),
            log(1),
            log(1),
            log(1),
            log(1)
          ]';
          
          vector[9] lambdas = [
            (1-omega_conc_sck[stsp])*thetas_idio[t][1], 
            (1-omega_conc_sck[stsp])*thetas_idio[t][2], 
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][1],
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][2],
            
            (1-omega_conc_sck[stsp])*thetas_idio[t][3], 
            omega_conc_sck[stsp]*(1-omega_shutdown[stsp])*thetas_idio[t][3],
            
            omega_conc_sck[stsp]*omega_shutdown[stsp]*thetas_idio[t][1],
            omega_conc_sck[stsp]*omega_shutdown[stsp]*thetas_idio[t][2],
            omega_conc_sck[stsp]*omega_shutdown[stsp]*thetas_idio[t][3]
          ]';
          
          vector[9] jlps = log(lambdas) + lpds;
          
          real jlp_all = log_sum_exp(jlps);
          
          vector[9] jp = exp(jlps - jlp_all);
          jp /= sum(jp); 
          
          sck_cat = categorical_rng(jp);
        }
        
        if (sck_cat == 1) { // no conc. effect, no idiosync. shock
          
          if (rw_obs[idx] >= epsilon) {
            log_rw_pred[idx] = normal_rng(mu_f[idx], sigma);
          } else { // sub-threshold observation
            log_rw_pred[idx] = normal_ub_rng(mu_f[idx], sigma, log(epsilon));
          }
          
        } else if (sck_cat == 2) { // no conc. effect, idiosync. depressed growth
          
          real log_rw;
          real residual;
          
          if (rw_obs[idx] >= epsilon) {
            
            log_rw = log(rw_obs[idx]);
            residual = log_rw - mu_f[idx];
            
          } else { // sub-threshold observation
            
            log_rw = normal_ub_rng(mu_f[idx], sqrt(tau_idio^2 + sigma^2), log(epsilon));
            residual = log_rw - mu_f[idx];
            
          }
          
          // we can reconstruct shock posterior using the normal-normal conjugancy
          real conjugate_mean = (tau_idio^2 / (tau_idio^2 + sigma^2)) * residual;
          real conjugate_sd   = sqrt((tau_idio^2 * sigma^2) / (tau_idio^2 + sigma^2));
          
          delta_sck[idx] = normal_rng(conjugate_mean, conjugate_sd);
          tree_idio_state[idx] = 1;
          
          log_rw_pred[idx] = normal_rng(mu_f[idx] + delta_sck[idx], sigma);
          
        } else if (sck_cat == 3) { // concordant depressed growth, no idiosync. shock
          
          real log_rw;
          real residual;
          
          if (rw_obs[idx] >= epsilon) {
            
            log_rw = log(rw_obs[idx]);
            residual = log_rw - mu_f[idx];
            
          } else { // sub-threshold observation
            
            log_rw = normal_ub_rng(mu_f[idx], sqrt(tau_conc[sp]^2 + sigma^2), log(epsilon));
            residual = log_rw - mu_f[idx];
            
          }
          
          real conjugate_mean = (tau_conc[sp]^2 / (tau_conc[sp]^2 + sigma^2)) * residual;
          real conjugate_sd   = sqrt((tau_conc[sp]^2 * sigma^2) / (tau_conc[sp]^2 + sigma^2));
          
          delta_sck[idx] = normal_rng(conjugate_mean, conjugate_sd);
          tree_conc_state[idx] = 1;
          
          log_rw_pred[idx] = normal_rng(mu_f[idx] + delta_sck[idx], sigma);
          
        } else if (sck_cat == 4) { // concordant depressed growth, idiosync. depressed growth
          
          real log_rw;
          real residual;
          
          if (rw_obs[idx] >= epsilon) {
            
            log_rw = log(rw_obs[idx]);
            residual = log_rw - mu_f[idx];
            
          } else { // sub-threshold observation
            
            log_rw = normal_ub_rng(mu_f[idx], sqrt(tau_idio^2 + tau_conc[sp]^2 + sigma^2), log(epsilon));
            residual = log_rw - mu_f[idx];
            
          }
          
          real conjugate_mean = ((tau_idio^2 + tau_conc[sp]^2) / (tau_idio^2 + tau_conc[sp]^2 + sigma^2)) * residual;
          real conjugate_sd   = sqrt(((tau_idio^2 + tau_conc[sp]^2) * sigma^2) / (tau_idio^2 + tau_conc[sp]^2 + sigma^2));
          
          delta_sck[idx] = normal_rng(conjugate_mean, conjugate_sd);
          tree_conc_state[idx] = 1;
          tree_idio_state[idx] = 1;
          
          log_rw_pred[idx] = normal_rng(mu_f[idx] + delta_sck[idx], sigma);
          
        } else if (sck_cat == 5) { // idiosync. shutdown (without any hidden concordant signal)
          
          log_rw_pred[idx] = not_a_number();
          tree_idio_state[idx] = 1;
          shutdown[idx] = 1;
          
        } else if (sck_cat == 6) { // idiosync. shutdown (with a hidden concordant depressed growth)
          
          log_rw_pred[idx] = not_a_number();
          tree_conc_state[idx] = 1;
          tree_idio_state[idx] = 1;
          shutdown[idx] = 1;
          
        } else if (sck_cat == 7) { // concordant shutdown (without any hidden idiosync. signal)
          
          log_rw_pred[idx] = not_a_number();
          tree_conc_state[idx] = 1;
          shutdown[idx] = 1;
          
        } else { // sck_cat 8 or 9, concordant shutdown (with a hidden idiosync. signal)
          
          log_rw_pred[idx] = not_a_number();
          tree_conc_state[idx] = 1;
          tree_idio_state[idx] = 1;
          shutdown[idx] = 1;
          
        }
      }
    }
  }
  
  // below: what fraction of each stand's total growth was lost due to shocks
    vector[N_stands] ratio_shock; 
    vector[N_stands] shock_decrease_perc; 
    vector[N_stands] growth_obs = rep_vector(0, N_stands); // as realized
    vector[N_stands] growth_cf = rep_vector(0, N_stands); // no concordant shocks

    for (t in 1:N_trees) {
      int st = stand_idxs[t];

      for (idx in tree_start_idxs[t]:tree_end_idxs[t]) {

        real base = exp(mu_f[idx]); // growth with no shock (counterfactual)
        growth_cf[st] += base;

        if (tree_conc_state[idx] == 1) {

          if (shutdown[idx] == 1) {
            growth_obs[st] += 0;
          } else {
            growth_obs[st] += exp(mu_f[idx] + delta_sck[idx]);
          }
        } else {
          growth_obs[st] += base; // idiosyncratic effects ignored
        }
      }
    }

    for (st in 1:N_stands) {
      ratio_shock[st] = growth_obs[st] / growth_cf[st];
      shock_decrease_perc[st] = 100 * (1 - ratio_shock[st]);
    }
    
}
