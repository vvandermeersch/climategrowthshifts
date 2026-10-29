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

  // For parallelization with reduce_sum (version en boucles, structure du modele d'origine)
  real loglikelihood_partial_sum(array[] int stand_ids_slice,
                              int start, int end,
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
                              array[] int species_idxs,
                              array[] int stand_species_idxs,
                              vector rw_obs,
                              vector gdd_obs,
                              vector pre_obs,
                              vector vpd_obs,
                              real gdd0,
                              real pre0,
                              real vpd0,
                              vector alpha_stand,
                              vector beta_gdd,
                              vector beta_pre,
                              vector beta_vpd,
                              array[] vector delta_clim,
                              vector kappa_clim,
                              vector f_tilde_ind,
                              array[] matrix L_cov_ind,
                              real epsilon,
                              real sigma,
                              real mu_omega_conc,
                              vector alpha_omega_species,
                              vector alpha_omega_stand,
                              real sigma_omega_event,
                              array[] vector alpha_tilde_omega_event,
                              vector omega_shutdown,
                              vector logit_phi_sck,
                              real beta_phi_vpd,
                              real beta_phi_pre,
                              array[] vector thetas_idio,
                              real sigma_idio,
                              vector sigma_conc,
                              vector shape_conc){

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

          vector[N_years[t]] f_ind = L_cov_ind[sp][1:N_years[t], 1:N_years[t]] * f_tilde_ind[tree_start:tree_end];

          vector[N_years[t]] mu = alpha_stand[st]
            + beta_gdd[sp] * (gdd_obs[tree_clim_idxs] - gdd0)
            + beta_pre[sp] * (pre_obs[tree_clim_idxs] - pre0)
            + beta_vpd[sp] * (vpd_obs[tree_clim_idxs] - vpd0)
            + kappa_clim[sp] * delta_clim[st, all_years_idxs_tree];

          for (y in 1:N_years[t]) {

            int idx = tree_start + y - 1;
            int ys = all_years_idxs_tree[y] - stand_start_years_idxs[st] + 1;

            real mu_f = mu[y] + f_ind[y];

            // omega now varies per year
            real omega_conc_sck = inv_logit(mu_omega_conc + alpha_omega_species[sp] + alpha_omega_stand[st]
              + sigma_omega_event * alpha_tilde_omega_event[stsp][all_years_idxs_tree[y]]);

            // epsilon censoring thresholding changed slightly (based on Ruben comment)
            if (rw_obs[idx] > epsilon * (1 + 1e-9)){

              real log_rw = log(rw_obs[idx]);

              // if you have a concordant shock, you can't have an idiosyncratic variation on top
              vector[3] lpds;
              lpds[1] = normal_lpdf(log_rw | mu_f, sigma); // no concordant shock, no idio. shock
              lpds[2] = normal_lpdf(log_rw | mu_f, sigma_idio); // no concordant shock, idio. shock
              lpds[3] = skew_normal_lpdf(log_rw | mu_f, sigma_conc[sp], shape_conc[sp]); // concordant depressed growth

              vector[3] lambdas;
              lambdas[1] = (1-omega_conc_sck)*thetas_idio[t][1];
              lambdas[2] = (1-omega_conc_sck)*thetas_idio[t][2];
              lambdas[3] = omega_conc_sck*(1-omega_shutdown[stsp]);
              // implicit: shutdown cannot happen (log(0))

              lpd_nonconc[ys] += log_sum_exp(
                log([thetas_idio[t][1], thetas_idio[t][2]]') +
                [lpds[1], lpds[2]]'
              );

              lpd_conc[ys] += log_sum_exp(log(lambdas) + lpds);

            }else{

              vector[5] lpds;
              lpds[1] = normal_lcdf(log(epsilon) | mu_f, sigma); // no concordant shock, no idio. shock
              lpds[2] = normal_lcdf(log(epsilon) | mu_f, sigma_idio); // no concordant shock, idio. shock
              lpds[3] = log(1); // no concordant shock, idiosync. shutdown
              lpds[4] = skew_normal_lcdf(log(epsilon) | mu_f, sigma_conc[sp], shape_conc[sp]); // concordant depressed growth
              lpds[5] = log(1); // concordant shutdown

              vector[5] lambdas;
              lambdas[1] = (1-omega_conc_sck)*thetas_idio[t][1];
              lambdas[2] = (1-omega_conc_sck)*thetas_idio[t][2];
              lambdas[3] = (1-omega_conc_sck)*thetas_idio[t][3];
              lambdas[4] = omega_conc_sck*(1-omega_shutdown[stsp]);
              lambdas[5] = omega_conc_sck*omega_shutdown[stsp];

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

          real phi_sck_y =  inv_logit(
          logit_phi_sck[st]
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
  real gdd0 = 12;
  real pre0 = 8;
  real vpd0 = 16;

  real epsilon = 1e-3; // measurement precision threshold (could be estimated...)

  array[N_stands] int stand_ids;
  for (s in 1:N_stands) stand_ids[s] = s;

}

parameters {

  real mu_alpha;

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

  array[N_stands] vector[N_all_years-1] delta_tilde_clim; // first delta_clim set to 0, delta_tilde_clim start at year=2
  real<lower=0,upper=1> rho_prim_clim;
  real<lower=0> tau_clim;
  vector<lower=0>[N_species - 1] kappa_clim_free; // sensibilite relative a l'espece 1 (reference)

  vector[N] f_tilde_ind; // tree-level func. behavior, short-term + long-term

  real mu_log_rho_short;
  real<lower=0> sigma_log_rho_short;
  vector<lower=0>[N_species] rho_short; // length scale

  real mu_log_gamma_short;
  real<lower=0> sigma_log_gamma_short;
  vector[N_species] log_gamma_short; // marginal variation

  real mu_log_delta_rho;
  real<lower=0> sigma_log_delta_rho;
  vector[N_species] log_delta_rho; // length scale

  real mu_log_gamma_long;
  real<lower=0> sigma_log_gamma_long;
  vector[N_species] log_gamma_long; // marginal variation

  // The (original) shocks!
  real mu_log_tau_conc;
  real<lower=0> sigma_log_tau_conc;
  vector<lower=0>[N_species] tau_conc; // log variation scale

  real mu_phi;
  real<lower=0> sigma_phi;
  vector[N_stands] logit_phi_sck; // baselines!
  real beta_phi_vpd; // effect of summer VPD on phi
  real beta_phi_pre; // effect of winter pre on phi

  real mu_omega_conc;
  real<lower=0> sigma_omega_species;
  real<lower=0> sigma_omega_stand;
  vector[N_species] alpha_omega_species;
  vector[N_stands] alpha_omega_stand;
  real<lower=0> sigma_omega_event; 
  array[N_stand_species] vector[N_all_years] alpha_tilde_omega_event; 

  real mu_omega_shutdown;
  real<lower=0> sigma_omega_shutdown;
  vector[N_stand_species] logit_omega_shutdown;

  // Idiosyncratic shocks!
  array[N_trees] simplex[3] thetas_idio;
  real<lower=0> tau_idio; // log variation scale

  real<lower=0> sigma; // proportional measurement error

}

transformed parameters {

  vector<lower=0, upper=1>[N_stand_species] omega_shutdown = inv_logit(logit_omega_shutdown);

  vector<lower=0>[N_species] rho_long = rho_short + exp(log_delta_rho);

  vector<lower=0>[N_species] gamma_short = exp(log_gamma_short);
  vector<lower=0>[N_species] gamma_long = exp(log_gamma_long);

  vector[N_species] kappa_clim = append_row(1, kappa_clim_free); // species 1 is the reference (kappa=1)

  real sigma_idio = sqrt(tau_idio^2 + sigma^2);
  vector[N_species] sigma_conc = sqrt(tau_conc^2 + sigma^2);

  vector[N_species] shape_conc = -tau_conc / sigma;

  array[N_species] matrix[N_all_years, N_all_years] L_cov_ind;
  {
    for(sp in 1:N_species){
      matrix[N_all_years, N_all_years] cov_ind
      = gp_exp_quad_cov(all_years, gamma_short[sp], rho_short[sp])
      + gp_exp_quad_cov(all_years, gamma_long[sp], rho_long[sp])
      + diag_matrix(rep_vector(1e-8, N_all_years));

      L_cov_ind[sp] = cholesky_decompose(cov_ind);
    }
  }

  // Random-walk prior
  real rho_clim = rho_prim_clim*2-1;
  array[N_stands] vector[N_all_years] delta_clim;
  for (s in 1:N_stands){
    delta_clim[s][1] = 0; // interpretation for n>1: changes in the mesoenvironment in stand s after the first year of observations
    for(i in 2:N_all_years){
      delta_clim[s][i] = rho_clim * delta_clim[s][i-1] + tau_clim * delta_tilde_clim[s][i-1]; // second indexing also i-1, because delta_tilde_clim is N_years-1
    }
  }

}

model {

  mu_alpha ~ normal(0, log(5)/2.32);

  sigma_alpha_stand ~ normal(0, log(5)/2.57);
  alpha_stand ~ normal(mu_alpha, sigma_alpha_stand);

  mu_beta_gdd ~ normal(0, log(1.3) / 2.57);
  sigma_beta_gdd ~ normal(0, log(1.3) / 2.57);
  beta_gdd ~ normal(mu_beta_gdd, sigma_beta_gdd);

  mu_beta_pre ~ normal(0, log(1.3) / 2.57);
  sigma_beta_pre ~ normal(0, log(1.3) / 2.57);
  beta_pre ~ normal(mu_beta_pre, sigma_beta_pre);

  mu_beta_vpd ~ normal(0, log(1.3) / 2.57);
  sigma_beta_vpd ~ normal(0, log(1.3) / 2.57);
  beta_vpd ~ normal(mu_beta_vpd, sigma_beta_vpd);

  mu_log_rho_short ~ normal(log(3.5), 0.25); // between 1.5 and 7 years
  sigma_log_rho_short ~ normal(0, 0.3); // up to 7 years
  rho_short ~ lognormal(mu_log_rho_short, sigma_log_rho_short);

  mu_log_gamma_short ~ normal(log(0.22), 0.58);
  sigma_log_gamma_short ~ normal(0, 0.5);
  log_gamma_short ~ normal(mu_log_gamma_short, sigma_log_gamma_short);

  mu_log_delta_rho ~ normal(log(20), 0.5);
  sigma_log_delta_rho ~ normal(0, 0.5);
  log_delta_rho ~ normal(mu_log_delta_rho, sigma_log_delta_rho); // will allow up to 60 years for rho_long

  mu_log_gamma_long ~ normal(log(0.42), 0.58);
  sigma_log_gamma_long ~ normal(0, 0.5);
  log_gamma_long ~ normal(mu_log_gamma_long, sigma_log_gamma_long);

  f_tilde_ind ~ std_normal();

  mu_log_tau_conc ~ normal(log(0.7), 0.25);
  sigma_log_tau_conc ~ normal(0, 0.5/2.57);
  tau_conc ~ lognormal(mu_log_tau_conc, sigma_log_tau_conc);

  tau_clim ~ normal(0, log(1.5) / 2.57); 
  rho_prim_clim ~ beta(4, 2);
  for (s in 1:N_stands){
    delta_tilde_clim[s] ~ std_normal(); // random-walk prior, ncp
  }
  kappa_clim_free ~ lognormal(-0.25, 0.6/2.32); // 2/3 <~ kappa <~ 3/2

  mu_phi ~ normal(-1.5, 1.5/2.57); // 0-50% (at VPD=16 and prec=8, not harsh conditions...)
  sigma_phi ~ normal(0, 0.5/2.57);
  logit_phi_sck ~ normal(mu_phi, sigma_phi);

  mu_omega_conc ~ normal(-2.5, 1.2/2.57); // 0-25%
  sigma_omega_species ~ normal(0, 2/2.57);
  sigma_omega_stand ~ normal(0, 3/2.57); 
  alpha_omega_species ~ normal(0, sigma_omega_species);
  alpha_omega_stand ~ normal(0, sigma_omega_stand);
  sigma_omega_event ~ normal(0, 3.5/2.57);
  for (ss in 1:N_stand_species) {
    alpha_tilde_omega_event[ss] ~ normal(0,1);
  }

  mu_omega_shutdown ~ normal(-5, 2); // 1-10%
  sigma_omega_shutdown ~ normal(0, 5/2.57); 
  logit_omega_shutdown ~ normal(mu_omega_shutdown, sigma_omega_shutdown);

  // not well thought
  beta_phi_vpd ~ normal(0, 0.4/2.57); //
  beta_phi_pre ~ normal(0, 0.4/2.57); //

  // Idiosyncratic shocks
  vector[3] thetas_baseline = [78, 27, 1]';
  real omega_thetas = 10;
  vector[3] alphas = thetas_baseline / omega_thetas;
  thetas_idio ~ dirichlet(alphas);

  tau_idio ~ normal(0, log(1.5) / 2.57);

  sigma ~ normal(0, log(1.1) / 2.57); // we expected 10% max

  profile("lkhd_out") {
    target += reduce_sum(
      loglikelihood_partial_sum,
      stand_ids,
      grainsize,
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
      species_idxs,
      stand_species_idxs,
      rw_obs,
      gdd_obs,
      pre_obs,
      vpd_obs,
      gdd0,
      pre0,
      vpd0,
      alpha_stand,
      beta_gdd,
      beta_pre,
      beta_vpd,
      delta_clim,
      kappa_clim,
      f_tilde_ind,
      L_cov_ind,
      epsilon,
      sigma,
      mu_omega_conc,
      alpha_omega_species,
      alpha_omega_stand,
      sigma_omega_event,
      alpha_tilde_omega_event,
      omega_shutdown,
      logit_phi_sck,
      beta_phi_vpd,
      beta_phi_pre,
      thetas_idio,
      sigma_idio,
      sigma_conc,
      shape_conc);
  }

}
