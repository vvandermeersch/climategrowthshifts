
data {
  
  int<lower=1> N;           // Total number of observations
  int<lower=1> N_species;   // Number of species
  int<lower=1> N_trees;     // Number of trees
  int<lower=1> N_all_years; // Max number of years
  int<lower=1> N_stands;    // Number of stands
  // int<lower=1> N_clades;    // Number of clades - for now, Gymnosperms=1 or Angiosperms=2
  
  array[N_stands] int<lower=1, upper=N> N_stand_trees; // Number of trees per stand
  array[N_stands] int<lower=1, upper=N_all_years> N_stand_years; // Max. number of observed years per stand
  array[N_stands] int<lower=1, upper=N_all_years> stand_start_years_idxs; // Indice of first year observed per stand
  
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
  vector<lower=0, upper=1>[N] w_sck;
  
}

transformed data {
  
  for (s in 1:N_stands) {
    
    array[N_stand_trees[s]] int stand_trees_idxs = linspaced_int_array(N_stand_trees[s], stand_trees_start_idxs[s], stand_trees_end_idxs[s]);
    vector[N_stand_years[s]] log_p0 = rep_vector(0, N_stand_years[s]);
    vector[N_stand_years[s]] log_p1 = rep_vector(0, N_stand_years[s]);
    
    for(ts in 1:N_stand_trees[s]){
      int t = stand_trees_idxs[ts];
      array[N_years[t]] int tree_idxs = linspaced_int_array(N_years[t], tree_start_idxs[t], tree_end_idxs[t]);
      array[N_years[t]] int all_years_idxs_tree = all_years_idxs[tree_idxs];
      
      for(y in 1:N_years[t]) {
        int ys = all_years_idxs_tree[y]-stand_start_years_idxs[s]+1;
        print("tree obs id=", tree_idxs[y], " - stand obs id=", ys);
        
      }
    }
    
    for(y in 1:N_stand_years[s]) {
      // target += log_mix(phi_sck, log_p1[y], log_p0[y]);
    }
  }
  
  
}

