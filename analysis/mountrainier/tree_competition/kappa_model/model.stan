data{
    int<lower=0> Nf; // number of focal tree
    array[Nf] int<lower=0> N; // number of species around each focal tree
    int <lower=0> N_total; // total number of neighbouring species(sum N)
    array[N_total] real<lower=0> b; // total basal area data of species
    array[N_total] real<lower=0, upper=1> focal_corr; // focal tree's correlation with neighbors
    array[Nf] int<lower=1,upper=N_total> start_idx; // start index of each focal tree
    array[Nf] int<lower=1,upper=N_total> end_idx; // end index of each focal tree
    array[Nf] real<lower=0> bf; // basal area of the focal tree
    array[Nf] real<lower=1e-10> y; // tree growth (mm)
}

transformed data {
   // baselines to have a meaningful interpretation of your parameters
   real bf0 = 2; // for a focal tree with a basal area of 2000cm2
   real BA_compet0 = 16; // surrounded by 8 trees of basal area 2000cm2 each
}

parameters{
    real<lower=0> sigma; //standard deviation of tree growth
    real<lower=0,upper=1> k; // strength of the species in competition
    real<lower=0> beta; // Impact of inter-species competition
    real<lower=0,upper=1> r; // scaling factor of basel area
    real<lower=1e-10> y0; // baseline growth (mm)
    // real y0_raw;
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
            baphy[i] += bn[t] * pow(fmax(corrn[t], 1e-10), k);
        }
        competition[i]=beta*(baphy[i]-BA_compet0);
        avails[i]=(bf[i]-bf0)*r;
        mu[i]=log(y0) + avails[i]-competition[i];
    }
}

model{
    sigma ~ normal(0, 0.095 / 2.57);
    //sigma ~ normal(0,1);
    y0 ~ normal(0,10/2.57);
    k ~ beta(2,2);
    r ~ beta(4,2);
    // regard 1000 cm^2 as 1 unit of BA
    beta ~ normal(0,log(1.1)/2.57);
    log(y) ~ normal(mu,sigma);
}




