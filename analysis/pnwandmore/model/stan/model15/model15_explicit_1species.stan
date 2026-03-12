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
                              vector f_tilde,
                              matrix L_cov,
                              real epsilon,
                              real sigma,
                              real tau_sck,
                              vector omega_conc_sck,
                              vector phi_sck,
                              vector delta_ind,
                              real pi_prone,
                              vector omega_ind_sck) {

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
            
            int tree_start = tree_start_idxs[t];
            int tree_end  = tree_end_idxs[t];
          
            vector[N_years[t]] f;
            profile("compute_f") {
              // f = block(L_cov[sp], 1, 1, N_years[t], N_years[t]) * f_tilde[tree_start:tree_end];
              f = L_cov[1:N_years[t], 1:N_years[t]] * f_tilde[tree_start:tree_end]; // 7% faster
            }
            
            array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_start:tree_end];
            array[N_years[t]] int tree_clim_idxs = stand_clim_idxs[all_years_idxs_tree];
            
            vector[N_years[t]] mu;
            profile("compute_mu") {
              mu = alpha
              + beta_gdd * (gdd_obs[tree_clim_idxs] - gdd0)
              + beta_pre * (pre_obs[tree_clim_idxs] - pre0)
              + beta_vpd * (vpd_obs[tree_clim_idxs] - vpd0)
              + f_sh[s,all_years_idxs_tree];
            }
            
            real log_prone = 0;
            real log_nonprone = 0; 
             
            for(y in 1:N_years[t]) {
              
              int idx = tree_start + y - 1;
              
              log_prone += log_mix(omega_ind_sck[s],
                normal_lpdf(delta_ind[idx] | 0, tau_sck),
                normal_lpdf(delta_ind[idx]  | 0, 1e-10));
              
              log_nonprone += normal_lpdf(delta_ind[idx] | 0, 1e-10);
              
            }
            
            lp += log_mix(pi_prone, log_prone, log_nonprone);
            
            
            profile("store_trees") {
              for(y in 1:N_years[t]) {
                
                int idx = tree_start + y - 1;
                int ys =  all_years_idxs[idx]-stand_start_years_idxs[s]+1;
                
                // real log_rw;
                real mu_f = mu[y] + f[y] + delta_ind[idx];
                
                if(rw_obs[idx] >= epsilon){
                  real log_rw = log(rw_obs[idx]);
                  log_p0[ys] += normal_lpdf(log_rw | mu_f, sigma);
                  // log_p0[ys] += log_mix(omega_nonconc_sck[stsp],
                  //                 normal_lpdf(log_rw | mu_f, sqrt(tau_sck[sp]^2 + sigma^2)),
                  //                 normal_lpdf(log_rw | mu_f, sigma));
                  log_p1[ys] += log_mix(omega_conc_sck[s],
                                  normal_lpdf(log_rw | mu_f, sqrt(tau_sck^2 + sigma^2)),
                                  normal_lpdf(log_rw | mu_f, sigma));
                }else{
                  log_p0[ys] += normal_lcdf(log(epsilon) | mu_f, sigma);
                  // log_p0[ys] += log_mix(omega_nonconc_sck[stsp],
                  //                 normal_lcdf(log(epsilon) | mu_f, sqrt(tau_sck[sp]^2 + sigma^2)),
                  //                 normal_lcdf(log(epsilon) | mu_f, sigma));
                  log_p1[ys] += log_mix(omega_conc_sck[s],
                                  normal_lcdf(log(epsilon) | mu_f, sqrt(tau_sck^2 + sigma^2)),
                                  normal_lcdf(log(epsilon) | mu_f, sigma));
                }
              }
            }
          }
        }
        
        for(y in 1:N_stand_years[s]) {
          lp += log_mix(phi_sck[s], log_p1[y], log_p0[y]);
        }
        
      }
    }

    return lp;
  }
  
}


data {
  
  int<lower=1> N;           // Total number of observations
  int<lower=1> N_trees;     // Number of trees
  int<lower=1> N_all_years; // Max number of years
  int<lower=1> N_stands;    // Number of stands
  int<lower=1> N_clades;    // Number of clades - for now, Gymnosperms=1 or Angiosperms=2
  
  array[N_stands] int<lower=1, upper=N> N_stand_trees; // Number of trees per stand
  array[N_stands] int<lower=1, upper=N> N_stand_years; // Max. number of observed years per stand
  array[N_stands] int<lower=1, upper=N> stand_start_years_idxs; // Indice of first year observed per stand
  
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
  real alpha;
  
  vector[N] f_tilde;
  
  // GDD slope (1/100degC)
  real beta_gdd;
  
  // SM slope (1/%)
  real beta_pre;
  
  // VPD slope (1/hPa)
  real beta_vpd;
  
  // Short-term proportional growth functional behavior
  array[N_stands] vector[N_all_years] f_tilde_sh; // Non-centered functional behavior
  real<lower=1> rho_sh;   // Time scale
  real<lower=0> gamma_sh; // Marginal variation - now fixed to 1! (and scaled by kappa)
  
  // Lifetime proportional growth scale (here I implement the hard contraint on both clade and species parameters?)
  real<lower=rho_sh> rho_sp; 
  
  // Lifetime proportional growth variation
  real<lower=0> gamma_sp; 
  
  // Growth shocks
  real<lower=0> tau_sck; // Outer yearly log variation scale (the shocks!)
  
  // Probability of stand-level shock 
  real<lower=0, upper=1> phi_sck0; // probability
  real<lower=0> tau_phi_sck; // log-odds
  vector[N_stands] alpha_phi_sck; // log-odds
  
  // Probability of tree-level shock given stand in shock (concordant shock)
  real<lower=0, upper=1> omega_conc_sck0; // probability
  real<lower=0> tau_omega_conc_sck; // log-odds
  vector[N_stands] alpha_omega_conc_sck; // log-odds
  
  // Individual shocks
  vector[N] delta_ind;
  
  // Probability of being prone to individual shocks
  real<lower=0, upper=1> pi_prone; 
  
  // Probability of tree-level shock when the tree is prone
  real<lower=0, upper=1> omega_ind_sck0; // probability
  real<lower=0> tau_omega_ind_sck; // log-odds
  vector[N_stands] alpha_omega_ind_sck; // log-odds
  
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
      =   gp_exp_quad_cov(all_years, gamma_sh, rho_sh)
      + diag_matrix(rep_vector(1e-10, N_all_years));
      L_cov_sh = cholesky_decompose(cov);
      
      for (s in 1:N_stands) {
        f_sh[s] = L_cov_sh * f_tilde_sh[s];
      }
    }
  }
  
  matrix[N_all_years, N_all_years] L_cov;
  profile("f_lg") {
        matrix[N_all_years, N_all_years] cov
        =  gp_exp_quad_cov(all_years, gamma_sp, rho_sp)
        + diag_matrix(rep_vector(1e-10, N_all_years));
        
        L_cov = cholesky_decompose(cov);
    
  }
  
  real mu_phi_sck = logit(phi_sck0); // log-odds
  // vector[N_stands] alpha_phi_sck = mu_phi_sck + tau_phi_sck*alpha_tilde_phi_sck; // log-odds
  vector<lower=0, upper=1>[N_stands] phi_sck = inv_logit(alpha_phi_sck); // probabilities
  
  real mu_omega_conc_sck = logit(omega_conc_sck0); // log-odds
  // vector[N_stand_species] alpha_omega_conc_sck = mu_omega_conc_sck + tau_omega_conc_sck*alpha_tilde_omega_conc_sck; // log-odds
  vector<lower=0, upper=1>[N_stands] omega_conc_sck = inv_logit(alpha_omega_conc_sck); // probabilities

  real mu_omega_ind_sck = logit(omega_ind_sck0); // log-odds
  vector<lower=0, upper=1>[N_stands] omega_ind_sck = inv_logit(alpha_omega_ind_sck); // probabilities
}

model {
  
  alpha ~ normal(0, log(5)/2.32); // 0.2 mm <~ exp(alpha) * 1 mm <~ 5 mm
  
  beta_gdd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_gsl <~ log(1.8)
  beta_pre ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_pre <~ log(1.8)
  beta_vpd ~ normal(0, log(1.8) / 2.57); // -log(1.8) <~ beta_vpd <~ log(1.8)
  
  rho_sp ~ lognormal(2.65, 0.135); // 10 <~ rho <~ 20
  
  gamma_sp ~ normal(0, log(10) / 2.57); // 0 <~ gamma <~ log(10)
  
  tau_sck ~ normal(0, log(20) / 2.57); 

  for (s in 1:N_stands)
    f_tilde_sh[s] ~ normal(0, 1);
  rho_sh ~ lognormal(1.7, 0.26);       // 3 <~ rho_sh <~ 10
  gamma_sh ~ normal(0, log(3) / 2.57); // 0 <~ gamma_sh <~ log(3)
  
  phi_sck0 ~ beta(2, 20); 
  tau_phi_sck ~ normal(0, 0.5); // 95Q of phi_sck ~ 0.25 
  alpha_phi_sck ~ normal(mu_phi_sck, tau_phi_sck);
  
  omega_conc_sck0 ~ beta(230, 14); // 0.9 <~ omega_conc_sck <~ 0.97 (most trees, but not ALL trees)
  tau_omega_conc_sck ~ normal(0, 0.3/2.57); // TO MODIFY!
  alpha_omega_conc_sck ~ normal(mu_omega_conc_sck, tau_omega_conc_sck); 
  
  pi_prone ~ beta(1, 72); 
  
  omega_ind_sck0 ~ beta(1, 72); // 0 <~ omega_conc_sck <~ 0.05 (should be rare, but... who knows?)
  tau_omega_ind_sck ~ normal(0, 0.3/2.57);
  alpha_omega_ind_sck ~ normal(mu_omega_ind_sck, tau_omega_ind_sck);
  
  sigma ~ normal(0, log(1.1) / 2.57);   // 0 <~ sigma <~ +log(1.1)
  
  profile("prepare") {
    f_tilde ~ normal(0,1);
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
      L_cov,
      epsilon,
      sigma,
      tau_sck,
      omega_conc_sck,
      phi_sck,
      delta_ind,
      pi_prone,
      omega_ind_sck);
   }
  
}
