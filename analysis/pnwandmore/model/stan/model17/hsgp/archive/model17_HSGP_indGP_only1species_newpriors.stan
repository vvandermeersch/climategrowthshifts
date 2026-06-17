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
  
  real normal_ub_rng(real mu, real sigma, real ub) {
    real p_ub = normal_cdf(ub | mu, sigma);
    if (p_ub == 0)
      return ub;
    real u = uniform_rng(0, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
  
  
  // HGSP functions
  vector diagSPD_EQ(real alpha, real rho, real L, int M) {
    return sqrt((alpha^2) * sqrt(2*pi()) * rho *
      exp(-0.5*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2));
  }

  matrix PHI_EQ(int N, int M, real L, vector x) {
    return sin(diag_post_multiply(
      rep_matrix(pi()/(2*L) * (x+L), M),
      linspaced_vector(M, 1, M)))/sqrt(L);
  }
  
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
                              real tau_sck,
                              vector omega_conc_sck,
                              real pi_idsc_sck,
                              vector phi_sck){

    real lp = 0;
    
    for (i in 1:(end-start+1)) {
  
      int s = stand_ids_slice[i];
      
      array[N_all_years] int stand_clim_idxs = linspaced_int_array(N_all_years,
        1+(s-1)*N_all_years, s*N_all_years);
      
      vector[N_stand_years[s]] log_p0 = rep_vector(0, N_stand_years[s]);
      vector[N_stand_years[s]] log_p1 = rep_vector(0, N_stand_years[s]);
      
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
          
          real lp_shock;
          real lp_noshock;
          real lp_doubleshock;
          
          if (rw_obs[idx] >= epsilon){
            real log_rw = log(rw_obs[idx]);
            lp_shock = normal_lpdf(log_rw | mu_f, sqrt(tau_sck^2 + sigma^2));
            lp_noshock = normal_lpdf(log_rw | mu_f, sigma);
            lp_doubleshock = normal_lpdf(log_rw | mu_f, sqrt(tau_sck^2 + tau_sck^2 + sigma^2));
          }else{
            lp_shock = normal_lcdf(log(epsilon) | mu_f, sqrt(tau_sck^2 + sigma^2));
            lp_noshock = normal_lcdf(log(epsilon) | mu_f, sigma);
            lp_doubleshock = normal_lcdf(log(epsilon) | mu_f, sqrt(tau_sck^2 + tau_sck^2 + sigma^2));
          }
          
          real lp_conc = log_mix(pi_idsc_sck, lp_doubleshock, lp_shock);
          real lp_nonconc = log_mix(pi_idsc_sck, lp_shock, lp_noshock);
          
          log_p0[ys] += lp_nonconc;
          log_p1[ys] += log_mix(omega_conc_sck[s], lp_conc, lp_nonconc);
        }
      }
      
      for (y in 1:N_stand_years[s]) {
        lp += log_mix(phi_sck[s], log_p1[y], log_p0[y]);
      }
    }

    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> N_trees;
  int<lower=1> N_all_years;
  int<lower=1> N_stands;
  int<lower=1> N_clades;
  
  array[N_stands] int<lower=1, upper=N> N_stand_trees;
  array[N_stands] int<lower=1, upper=N> N_stand_years;
  array[N_stands] int<lower=1, upper=N> stand_start_years_idxs;
  
  array[N_trees] int<lower=1, upper=N_stands> stand_idxs;
  array[N_trees] int<lower=1, upper=N> N_years;
  array[N_all_years] real all_years;
  array[N] real years;
  array[N] int<lower=1, upper=N_all_years> all_years_idxs;
  
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;
  
  array[N_stands] int<lower=1, upper=N_trees> stand_trees_start_idxs;
  array[N_stands] int<lower=1, upper=N_trees> stand_trees_end_idxs;
  array[N_trees] int<lower=1, upper=N_trees> stand_tree_idxs;
  
  vector[N] rw_obs;
  
  vector[N_stands*N_all_years] gdd_obs;
  vector[N_stands*N_all_years] pre_obs;
  vector[N_stands*N_all_years] vpd_obs;
  
  int<lower=1> grainsize;
}

transformed data {
  real gdd0 = 10;
  real pre0 = 5;
  real vpd0 = 23;
  real epsilon = 1e-3;
  
  int M = 20;
  real L_sp = 1.5 * N_all_years;
  matrix[N_all_years, M] PHI_sp = PHI_EQ(N_all_years, M, L_sp, to_vector(all_years));
  
  array[N_stands] int stand_ids;
  for (s in 1:N_stands) stand_ids[s] = s;
}

parameters {
  real alpha;
  
  matrix[M, N_trees] f_tilde;
  vector[N] f_ind_tilde;
  
  real beta_gdd;
  real beta_pre;
  real beta_vpd;
  
  array[N_stands] vector[N_all_years] f_tilde_sh;
  real<lower=1> rho_sh;
  real<lower=0> gamma_sh;
  
  real<lower=1> rho_ind;
  real<lower=0> gamma_ind;
  
  real<lower=1> rho_sp;
  real<lower=0> gamma_sp;
  
  real<lower=0> tau_sck;
  
  real<lower=0, upper=1> phi_sck0;
  real<lower=0> tau_phi_sck;
  vector[N_stands] alpha_phi_sck_tilde;
  
  real<lower=0, upper=1> omega_conc_sck0;
  real<lower=0> tau_omega_conc_sck;
  vector[N_stands] alpha_omega_conc_sck;
  
  // real mu_logdelta_omega_nonconc_sck;
  // real<lower=0> tau_logdelta_omega_nonconc_sck;
  // vector[N_stands] logdelta_omega_nonconc_sck_tilde;
  real<lower=0, upper=1> pi_idsc_sck;
  
  real<lower=0> sigma;
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
  
  // Hilbert space approximate GP
  vector[M] sqrt_spd_sp = diagSPD_EQ(gamma_sp, rho_sp, L_sp, M);
  
  real mu_phi_sck = logit(phi_sck0);
  vector[N_stands] alpha_phi_sck = mu_phi_sck + tau_phi_sck * alpha_phi_sck_tilde;
  vector<lower=0, upper=1>[N_stands] phi_sck = inv_logit(alpha_phi_sck);
  
  real mu_omega_conc_sck = logit(omega_conc_sck0);
  vector<lower=0, upper=1>[N_stands] omega_conc_sck = inv_logit(alpha_omega_conc_sck);
  
  // vector[N_stands] logdelta_omega_nonconc_sck = mu_logdelta_omega_nonconc_sck +
  //   tau_logdelta_omega_nonconc_sck * logdelta_omega_nonconc_sck_tilde;
  // vector[N_stands] delta_omega_nonconc_sck = exp(logdelta_omega_nonconc_sck);
  
  // vector[N_stands] alpha_omega_nonconc_sck = alpha_omega_conc_sck - delta_omega_nonconc_sck;
  // vector<lower=0, upper=1>[N_stands] omega_nonconc_sck = inv_logit(alpha_omega_nonconc_sck);

}

model {
  
  alpha ~ normal(0, log(5)/2.32);
  
  beta_gdd ~ normal(0, log(1.8) / 2.57);
  beta_pre ~ normal(0, log(1.8) / 2.57);
  beta_vpd ~ normal(0, log(1.8) / 2.57);
  
  rho_sp ~ lognormal(3.7, 0.35); // before: lognormal(3.55, 0.24)
  gamma_sp ~ normal(0, log(10) / 2.57);
  
  tau_sck ~ normal(0, 5 / 2.57); // before: normal(0, log(20) / 2.57)

  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26);
  gamma_sh ~ normal(0, log(3) / 2.57);
  
  rho_ind ~ lognormal(1.4, 0.35); // before: lognormal(0.80, 0.40)
  gamma_ind ~ normal(0, log(3) / 2.57);
  
  phi_sck0 ~ beta(2, 20);
  tau_phi_sck ~ normal(0, 1);
  alpha_phi_sck_tilde ~ normal(0, 1);
  
  omega_conc_sck0 ~ beta(12.28, 34); // before:  beta(30, 20)
  tau_omega_conc_sck ~ normal(0, 1);
  alpha_omega_conc_sck ~ normal(mu_omega_conc_sck, tau_omega_conc_sck);
  
  // mu_logdelta_omega_nonconc_sck ~ normal(log(5), log(3)/2.57);
  // tau_logdelta_omega_nonconc_sck ~ normal(0, 1/2.57);
  // logdelta_omega_nonconc_sck_tilde ~ normal(0, 1);
  
  pi_idsc_sck ~ beta(2, 100); // between 0 an 5%...
  
  sigma ~ normal(0, log(1.1) / 2.57);
  
  to_vector(f_tilde) ~ normal(0, 1);
  f_ind_tilde ~ normal(0, 1);
  
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
      tau_sck,
      omega_conc_sck,
      pi_idsc_sck,
      phi_sck);
  }
}
