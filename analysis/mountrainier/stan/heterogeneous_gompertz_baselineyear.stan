functions {
  real gompertz(real t, real d_max, real beta, real delta_t_half, real baseline) {
    return d_max * exp(-log2() * exp( -   (2 / log2()) 
                                      * (beta / d_max) 
                                      * (t - (baseline + delta_t_half))));
  }
}

data {
  int<lower=1> N_trees; // number of trees
  int<lower=N_trees> N; // number of observations
  int<lower=0, upper=N> Sobs[N_trees+1]; // range of no. of observations per trees
  vector[N] years;
  vector[N] dbhs;
  
  int<lower=1> N_year_grid;
  vector[N_year_grid] year_grid;
}

transformed data {
  
  vector[N_trees] baseline; // different baseline for each tree
  for (i in 1:N_trees) {
    
    int start_idx = Sobs[i]+1;
    baseline[i] = years[start_idx]; // first year of observations
    
  }
  
}

parameters {
  real<lower=0> d0[N_trees];      // Minimum DBH (cm)
  real<lower=0> delta_d[N_trees]; // Difference between minimum and maximum DBH (cm)
  real<lower=0> beta[N_trees];    // Intermediate linear growth rate (cm / year)
  real delta_t_half[N_trees];     // Time to half maximum DBH (years relative to baseline)
  real<lower=0> sigma;            // Measurement variability (cm)
}

model {
  d0 ~ normal(0, 150 / 2.57);          // 99% prior mass between 0 and 150 cm
  delta_d ~ normal(0, 50 / 2.57);      // 99% prior mass between 0 and 50 cm
  beta ~ normal(0, 5 / 2.57);          // 99% prior mass between 0 and 5 cm / year
  delta_t_half ~ normal(0, 40 / 2.32); // 99% prior mass between +/- 40 years
  sigma ~ normal(0, 0.5 / 2.57);      // 99% prior mass between 0 and 0.5 cm 
  
  for (i in 1:N_trees) {
    
    int start_idx = Sobs[i]+1;
    int end_idx = Sobs[i+1];
    int N_obs = Sobs[i+1]-Sobs[i];
    vector[N_obs] years_i = years[start_idx:end_idx];
    vector[N_obs] dbhs_i = dbhs[start_idx:end_idx];
    
    for (n in 1:N_obs) {
      
      real mu =  gompertz(years_i[n], delta_d[i], beta[i], delta_t_half[i], baseline[i]) 
               + d0[i];
      dbhs_i[n] ~ normal(mu, sigma);
    }
  }
}

generated quantities {
  real dbh_grid[N_trees, N_year_grid];
  real dbh_grid_pred[N_trees, N_year_grid];
  for (t in 1:N_trees) {
    for (n in 1:N_year_grid) {
      dbh_grid[t, n] = gompertz(year_grid[n], delta_d[t], 
                       beta[t], delta_t_half[t], baseline[t]) + d0[t];
      dbh_grid_pred[t, n] = normal_rng(dbh_grid[t, n], sigma);
    }
  }
}

