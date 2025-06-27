//This is the beginning of phylo analysis
//Created by Avery, June 25, 2025

data {
  int<lower=0> N_trees; //number of trees
  int S; //number of species
  array[N_trees] real y; //ring width observations
  int<lower=0> N_neighbors; //total number of neighbors
  
  array[N_trees] int<lower=1> tree_sp; //species of focal tree
  array[N_trees] int<lower=1> tree_N_neighbs;
  array[N_trees] int<lower=1, upper=N_neighbors> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N_neighbors> tree_end_idxs;
  
  vector[N_neighbors] neighbor_BA; 
  array[N_neighbors] int<lower=1, upper = S> neighbor_sp; 
  real BA0;
  //array[N_neighbors] real BA0; //baseline BA for neighbors
  
  matrix[S,S] Sigma; //species trait covariance
}

transformed data {
  matrix[S,S] L = cholesky_decompose(Sigma);
  
  vector[N_neighbors] neighbor_relBA = neighbor_BA/BA0; 
  
}

parameters {
  array[S] real mu; //base growth mean
  array[S] real<lower = 0> base_var; //base growth variance
  real<lower=0> alpha; //effect of neighbor BA
  real beta1; //slope of competition on growth
  real beta2; //factor of trait distance in phylo term
  
  vector[S] traits; //latent species traits
  //vector basegrowth[N_trees]; //base unobserved growth, unneeded?
  real kappa; //effect of hylogenetic competition at no trait distance
}


model {
  //Priors
  alpha ~ normal(0, 0.1);
  beta1 ~ normal(1, 1);
  beta2 ~ normal(1, 1);
  for(i in 1:S){
      mu[i] ~ normal(0.7, 0.3);
      base_var[i] ~ normal(0.5, 0.3);
  }

  traits ~ multi_normal_cholesky(zeros_vector(S), L);
  kappa ~ normal(0.1, 0.1);

  target += multi_normal_cholesky_lpdf(traits | zeros_vector(S), L);

  for(i in 1:N_trees) {
    //For each tree, get the product of BAnb^alphas, the product of phylogenetic competition,
    //and multiply those
      vector[tree_N_neighbs[i]] BAcomp = neighbor_relBA[tree_start_idxs[i]:tree_end_idxs[i]]^(-alpha);
      vector[tree_N_neighbs[i]] phylocomp = 1- kappa*exp(-beta2*(abs(
        traits[neighbor_sp[tree_start_idxs[i]:tree_end_idxs[i]]]-rep_vector(traits[tree_sp[i]], tree_N_neighbs[i])
        )));
        
      real totalcomp = beta1*prod(BAcomp)*prod(phylocomp);
      
      #If B~lognormal(u, sd), then a*B = C~lognormal(u + log(a), sd)
      target += lognormal_lpdf(y[i] | mu[tree_sp[i]] + log(totalcomp), base_var[tree_sp[i]]);
  }
}
