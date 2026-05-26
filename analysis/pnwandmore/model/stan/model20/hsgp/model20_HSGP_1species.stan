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
                              vector rw_obs,
                              vector gdd_obs,
                              vector pre_obs,
                              vector vpd_obs,
                              real gdd0,
                              real pre0,
                              real vpd0,
                              real alpha,
                              real beta_gdd,
                              real beta_pre,
                              real beta_vpd,
                              array[] vector f_sh,
                              matrix f_tilde,
                              vector f_ind_tilde,
                              matrix L_cov_ind,
                              matrix PHI_sp,
                              vector sqrt_spd_sp,
                              real epsilon,
                              real sigma,
                              real tau_small,
                              real mu_sck,
                              real tau_sck,
                              array[] vector thetas_conc,
                              vector omega_shutdown,
                              vector phi_conc,
                              array[] vector thetas_idio,
                              real tau_idio){

    real lp = 0;
    
    for (i in 1:(end-start+1)) {
  
      int s = stand_ids_slice[i];
      
      array[N_all_years] int stand_clim_idxs = linspaced_int_array(N_all_years,
        1+(s-1)*N_all_years, s*N_all_years);
      
      vector[N_stand_years[s]] lpd_nonconc = rep_vector(0, N_stand_years[s]);
      vector[N_stand_years[s]] lpd_conc = rep_vector(0, N_stand_years[s]);
      
      for (t in stand_tree_idxs[stand_trees_start_idxs[s]:stand_trees_end_idxs[s]]) {
        
        int tree_start = tree_start_idxs[t];
        int tree_end = tree_end_idxs[t];
        
        array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_start:tree_end];
        array[N_years[t]] int tree_clim_idxs = stand_clim_idxs[all_years_idxs_tree];

        vector[N_years[t]] f = PHI_sp[all_years_idxs_tree, ] * (sqrt_spd_sp .* f_tilde[, t]);
        vector[N_years[t]] f_ind = L_cov_ind[1:N_years[t], 1:N_years[t]] * f_ind_tilde[tree_start:tree_end];
        
        vector[N_years[t]] mu = alpha
          + beta_gdd * (gdd_obs[tree_clim_idxs] - gdd0)
          + beta_pre * (pre_obs[tree_clim_idxs] - pre0)
          + beta_vpd * (vpd_obs[tree_clim_idxs] - vpd0)
          + f_sh[s, all_years_idxs_tree];
        
        for (y in 1:N_years[t]) {
          
          int idx = tree_start + y - 1;
          int ys = all_years_idxs_tree[y] - stand_start_years_idxs[s] + 1;
          
          real mu_f = mu[y] + f[y] + f_ind[y];
          
          if (rw_obs[idx] >= epsilon){
            
            real log_rw = log(rw_obs[idx]);
            
            vector[6] lpds = [
              normal_lpdf(log_rw | mu_f, sigma), // no concordant shock, no idio. shutdown
              normal_lpdf(log_rw | mu_f, sqrt(tau_idio^2 + sigma^2)), // no concordant shock, idio. shock
              
              normal_lpdf(log_rw | mu_f, sqrt(tau_small^2 + sigma^2)), // concordant small resp., no idio. shutdown
              normal_lpdf(log_rw | mu_f, sqrt(tau_idio^2 + tau_small^2 + sigma^2)), // concordant small resp., idio. shock
              
              normal_lpdf(log_rw | mu_f - mu_sck, sqrt(tau_sck^2 + sigma^2)), // concordant extreme resp., no idio. shutdown
              normal_lpdf(log_rw | mu_f - mu_sck, sqrt(tau_idio^2 + tau_sck^2 + sigma^2)) // concordant extreme resp., idio. shock
            ]';
            
            vector[6] lambdas = [
              thetas_conc[s][1]*thetas_idio[t][1], 
              thetas_conc[s][1]*thetas_idio[t][2], 
              
              thetas_conc[s][2]*thetas_idio[t][1],
              thetas_conc[s][2]*thetas_idio[t][2],
              
              thetas_conc[s][3]*(1-omega_shutdown[s])*thetas_idio[t][1],
              thetas_conc[s][3]*(1-omega_shutdown[s])*thetas_idio[t][2]
            ]'; // implicit:  shutdown cannot happen (log(0))
            
            lpd_nonconc[ys] += log_sum_exp(
              log([thetas_idio[t][1], thetas_idio[t][2]]') +
              [lpds[1], lpds[2]]'
            );
            
            lpd_conc[ys] += log_sum_exp(log(lambdas) + lpds);

          }else{
            
            vector[12] lpds = [
              normal_lcdf(log(epsilon) | mu_f, sigma), // no concordant shock, no idio. shutdown
              normal_lcdf(log(epsilon) | mu_f, sqrt(tau_idio^2 + sigma^2)), // no concordant shock, idio. shock
              log(1), // no concordant shock, idiosync. shutdown
              
              normal_lcdf(log(epsilon)| mu_f, sqrt(tau_small^2 + sigma^2)), // concordant small resp., no idio. shutdown
              normal_lcdf(log(epsilon) | mu_f, sqrt(tau_idio^2 + tau_small^2 + sigma^2)), // concordant small resp., idio. shock
              log(1), // concordant small resp., idiosync. shutdown
              
              normal_lcdf(log(epsilon)| mu_f - mu_sck, sqrt(tau_sck^2 + sigma^2)), // concordant extreme resp., no idio. shutdown
              normal_lcdf(log(epsilon) | mu_f - mu_sck, sqrt(tau_idio^2 + tau_sck^2 + sigma^2)), // concordant extreme resp., idio. shock
              log(1), // concordant extreme resp., idiosync. shutdown
              
              log(1), // concordant shutdown, no idiosync. shock
              log(1), // concordant shutdown, idio. shock
              log(1) // concordant shutdown, idiosync. shutdown
            ]';

            vector[12] lambdas = [
              thetas_conc[s][1]*thetas_idio[t][1], 
              thetas_conc[s][1]*thetas_idio[t][2], 
              thetas_conc[s][1]*thetas_idio[t][3], 
              
              thetas_conc[s][2]*thetas_idio[t][1], 
              thetas_conc[s][2]*thetas_idio[t][2], 
              thetas_conc[s][2]*thetas_idio[t][3], 
              
              thetas_conc[s][3]*(1-omega_shutdown[s])*thetas_idio[t][1],
              thetas_conc[s][3]*(1-omega_shutdown[s])*thetas_idio[t][2],
              thetas_conc[s][3]*(1-omega_shutdown[s])*thetas_idio[t][3],
              
              thetas_conc[s][3]*omega_shutdown[s]*thetas_idio[t][1],
              thetas_conc[s][3]*omega_shutdown[s]*thetas_idio[t][2],
              thetas_conc[s][3]*omega_shutdown[s]*thetas_idio[t][3]
            ]';
            
            lpd_nonconc[ys] += log_sum_exp(
              log([thetas_idio[t][1], thetas_idio[t][2], thetas_idio[t][3]]') +
              [lpds[1], lpds[2], lpds[3]]'
            );
            
            lpd_conc[ys] += log_sum_exp(log(lambdas) + lpds);
          }
        }
      }
      
      for (y in 1:N_stand_years[s]) {
        lp += log_mix(phi_conc[s], lpd_conc[y], lpd_nonconc[y]);
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
}

transformed data {
  // baselines
  real gdd0 = 10;
  real pre0 = 5;
  real vpd0 = 23;
  
  real epsilon = 1e-3; // measurement precision threshold (could be estimated...)
  
  // HS approximation (see Riutort-Mayol et al. 2023)
  int M = 20; // number of inducing points
  real L_sp = 1.5 * N_all_years; // boundary conditions
  matrix[N_all_years, M] PHI_sp = PHI(N_all_years, M, L_sp, to_vector(all_years)); // doesn't depend on parameters!
  
  array[N_stands] int stand_ids;
    for (s in 1:N_stands) stand_ids[s] = s;
}

parameters {
  real alpha; // log(ring width) baseline
  
  real beta_gdd; // GDD slope (1/kdegC)
  real beta_pre; // Precipitation slope (1/dm)
  real beta_vpd; // VPD slope (1/hPa)
  
  array[N_stands] vector[N_all_years] f_tilde_sh; // short-term proportional growth functional behavior
  real<lower=1> rho_sh; // length scale
  real<lower=0> gamma_sh; // marginal variation
  
  vector[N] f_ind_tilde; // short-term tree-level func. behavior(canopy dynamics)
  real<lower=1> rho_ind; // length scale
  real<lower=0> gamma_ind; // marginal variation
  
  matrix[M, N_trees] f_tilde; // mid- to long-term tree-level func. behavior (allometry) -- approx. with HS
  real<lower=rho_ind> rho_sp; // length scale (should vary per species) + ADDED lower constraint
  real<lower=0> gamma_sp; // marginal variation (should vary per species)
  
  real<lower=0> mu_sck;
  real<lower=0> tau_sck; // log variation scale (the shocks!)
  real<lower=0> tau_small; // log variation scale (the shocks!)
  
  // Probabilities of stand-level 'concordant state' (partially pooled by stands)
  vector<lower=0, upper=1>[N_stands] phi_conc;
  
  // Probabilities of tree-level shock given stand in 'concordant state' (partially pooled by stands)
  array[N_stands] simplex[3] thetas_conc;
  
  // Probabilities of tree-level shutdown given stand in 'concordant state' and tree in a 'shock state' (partially pooled by stands)
  vector<lower=0, upper=1>[N_stands] omega_shutdown;
  
  // Idiosyncratic shocks!
  array[N_trees] simplex[3] thetas_idio;
  real<lower=0> tau_idio; // log variation scale
  
  real<lower=0> sigma; // proportional measurement error
}

transformed parameters {
  
  array[N_stands] vector[N_all_years] f_sh;
  matrix[N_all_years, N_all_years] L_cov_sh;
  profile("L_cov_fsh") {
    {
      matrix[N_all_years, N_all_years] cov
        = gp_exp_quad_cov(all_years, gamma_sh, rho_sh)
        + diag_matrix(rep_vector(1e-10, N_all_years));
      L_cov_sh = cholesky_decompose(cov);
      for (s in 1:N_stands)
        f_sh[s] = L_cov_sh * f_tilde_sh[s];
    }
  }
  
  matrix[N_all_years, N_all_years] L_cov_ind;
  {
    matrix[N_all_years, N_all_years] cov_ind
      = gp_exp_quad_cov(all_years, gamma_ind, rho_ind)
      + diag_matrix(rep_vector(1e-10, N_all_years));
    L_cov_ind = cholesky_decompose(cov_ind);
  }
  
  // Hilbert space approximatino for long-term GP
  vector[M] sqrt_spd_sp = diagSPD_EQ(gamma_sp, rho_sp, L_sp, M);

}

model {
  
  alpha ~ normal(0, log(5)/2.32); // 0.2 mm < exp(alpha) * 1 mm < 5 mm
  
  beta_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) < beta_gsl < log(1.8)
  beta_pre ~ normal(0, log(1.8) / 2.57); // -log(1.8) < beta_pre < log(1.8)
  beta_vpd ~ normal(0, log(1.8) / 2.57); // -log(1.8) < beta_vpd < log(1.8)
  
  rho_sp ~ lognormal(3.7, 0.35); // 20 < rho < 90
  gamma_sp ~ normal(0, log(10) / 2.57); // 0 < gamma < log(10)
  
  mu_sck ~ lognormal(0.154, 0.546); // log(1.5) to log(30)
  tau_sck ~ normal(0, 1.5 / 2.57); // 
  
  tau_small ~ normal(0, log(1.5) / 2.57); // 

  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26); // 3 <~ rho_sh <~ 10
  // rho_sh ~ lognormal(0.4, 0.3);
  gamma_sh ~ normal(0, log(3) / 2.57); // 0 < gamma_sh < log(3)
  
  rho_ind ~  lognormal(1.4, 0.25); // 2 < rho_ind < 7
  gamma_ind ~ normal(0, log(3) / 2.57); // 0 < gamma_ind < log(3)
  
  phi_conc ~ beta(2, 6); // TO COMPLETE
  
  // Idiosyncratic shocks
  vector[3] thetas_baseline_conc = [5, 100, 20]';  // TO COMPLETE
  real omega_thetas_conc = 3;
  vector[3] ones = rep_vector(1, 3);
  vector[3] alphas_conc = thetas_baseline_conc / omega_thetas_conc + ones;
  for (s in 1:N_stands) {
    thetas_conc[s] ~ dirichlet(alphas_conc);
  }
  
  omega_shutdown ~ beta(1.15, 4.3); // 0% < omega_shutdown0 < 60%
  
  sigma ~ normal(0, log(1.1) / 2.57); // we expected 10% max.
  
  to_vector(f_tilde) ~ normal(0, 1);
  f_ind_tilde ~ normal(0, 1);
  
  // Idiosyncratic shocks
  vector[3] thetas_baseline = [100, 5, 0.5]';
  real omega_thetas = 4;
  vector[3] alphas = thetas_baseline / omega_thetas + ones;
  for (t in 1:N_trees) {
    thetas_idio[t] ~ dirichlet(alphas);
  }
  
  tau_idio ~ normal(0, log(2) / 2.57);  // TO COMPLETE
  
  profile("likelihood") {
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
      f_sh,
      f_tilde,
      f_ind_tilde,
      L_cov_ind,
      PHI_sp,
      sqrt_spd_sp,
      epsilon,
      sigma,
      tau_small,
      mu_sck,
      tau_sck,
      thetas_conc,
      omega_shutdown,
      phi_conc,
      thetas_idio,
      tau_idio);
  }
}
