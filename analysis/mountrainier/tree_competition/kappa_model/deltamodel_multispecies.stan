data{
    int<lower=0> Nf; // number of focal tree
    int<lower=0> S; // *number of species
    array[Nf] int<lower=0> N; // number of species around each focal tree
    array[Nf] int<lower=0> tree_sp; // *focal tree species
    int <lower=0> N_total; // total number of neighbouring species(sum N)
    array[N_total] real<lower=0> b; // total basal area data of species
    array[N_total] real<lower=0, upper=1> focal_corr; // focal tree's correlation with neighbors (phylo distance, should we use trait distance instead?)
    array[Nf] int<lower=1,upper=N_total> start_idx; // start index of each focal tree
    array[Nf] int<lower=1,upper=N_total> end_idx; // end index of each focal tree
    array[Nf] real<lower=0> ra; // basal radius of the focal tree (*10 cm)
    array[Nf] real<lower=1e-10> y; // tree growth (mm)
}

transformed data {
   // baselines to have a meaningful interpretation of your parameters
   real r0 = 2.5;
   //real bf0 = 2; // for a focal tree with a basal area of 2000cm2
   real BA_compet0 = 16; // surrounded by 8 trees of basal area 2000cm2 each
}

parameters{
    real<lower=0> sigma; //*standard deviation of tree growth
    real<lower=0,upper=1> k; // *strength of the species in competition
    array[S] real<lower=0> beta; // Impact of inter-species competition
    array[S] real<lower=0,upper=1> r; // scaling factor of basel area
   array[S] real<lower=1e-10> y0; // *baseline growth (mm)
   //real<lower=0> y0_mu;
   //real<lower=0> y0_tau;
   //real<lower=0> beta_mu;
  // real<lower=0> beta_tau;
   //real<lower=0,upper=1> r_mu;
   //real<lower=0> r_phi;
}

transformed parameters{
    // real<lower=1e-10> y0 = exp(y0_raw); // baseline growth (mm)
    array[Nf] real baphy; // phelogeny distance
    //array[Nf] real<lower=0> competition; // inter-species competition
    array[Nf] real competition;
    //array[Nf] real<lower=0> avails;//available resources of focal tree
    array[Nf] real avails;
    array[Nf] real mu;

    for (i in 1:Nf){
        int len = end_idx[i]-start_idx[i]+1; // the length of bn and corrn for each focal tree 
        vector[len] bn = segment(to_vector(b),start_idx[i],len); // the total BA of each neighboring species of each focal tree
        vector[len] corrn = segment(to_vector(focal_corr),start_idx[i],len); // each focal tree's correlation with neighbours
        baphy[i]=0;
        for (t in 1:N[i]) {
            //baphy[i] += bn[t]*pow(corrn[t],k);
            baphy[i] += bn[t] * pow(fmax(corrn[t], 1e-2), k);
        }
        competition[i]=beta[tree_sp[i]]*(baphy[i]-BA_compet0);
        avails[i]=(log(ra[i])-log(r0))*r[tree_sp[i]];
        mu[i]=log(y0[tree_sp[i]]) + avails[i]-competition[i];
    }
}

model{
    sigma ~ normal(0, 0.095 / 2.57);
    y0 ~ normal(0,10/2.57);
    //y0_mu ~ normal(0,9/2.57);
    //y0_tau ~ normal(0,2/2.57);

    //beta_mu ~ normal(0,log(1.1)/2.57);
    //beta_tau ~normal(0,log(1.1)/(25.7)); // ** these priors make sense? 

    //r_mu ~ beta(8,4);
    //r_phi ~ normal(13,2) T[0,];
    
    //for (i in 1:S) {
        //y0[i]~normal(y0_mu,y0_tau) T[0,];
        //beta[i]~normal(beta_mu,beta_tau) T[0,];
        //r[i]~beta(r_mu*r_phi,(1-r_mu)*r_phi);
    //}

    k ~ beta(2,2);
    r ~ beta(4,2);
    // regard 1000 cm^2 as 1 unit of BA
    beta ~ normal(0,log(1.1)/2.57);
    for (i in 1:Nf) {
        log(y[i]) ~ normal(mu[i],sigma);
    }
}

generated quantities {
    array[Nf] real y_rep;
    for (i in 1:Nf) {
        y_rep[i] = lognormal_rng(mu[i],sigma);
    }
}

