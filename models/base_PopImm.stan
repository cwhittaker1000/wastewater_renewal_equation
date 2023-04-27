data {
  // Inputs Relating to the Infection Process
  int<lower=1> N0;                 // number of days to populate with initial_infections
  int<lower=1> N2;                 // days of observed data which we'll do proper inference on

  // Inputs Relating to the Deaths Observation Process 
  array[N2] int deaths;            // reported deaths
  vector[N2] death_delay;          // death delay distribution (time between infection and death occurring)
  int Death_Fitting_Start;         // time at which point you start fitting to deaths
  array[N2] real SI;               // fixed, pre-calculated Serial Interval
  
  // Inputs Relating to Time-Varying Rt
  int W;                           // number of unique periods for different Rts
  array[N2] int week_index;        // index specifying which Rt period a particular day's observation belongs to
  
  // Miscellaneous Inputs
  int SI_cutoff;                   // number of days after which SI is cutoff (for computational efficiency purposes)
  real<lower=0,upper=1> ifr_mean;  // mean IFR for prior
  real<lower=0> ifr_sd;            // sd IFR for prior    
  real pop;                        // population for immunity adjustment
}

transformed data {
  vector[SI_cutoff] SI_rev;              // SI in reverse order
  vector[SI_cutoff] death_delay_rev;     // death delay distribution in reversed order
  
  for(i in 1:SI_cutoff)                  // reversing the SI and death delay distribution - done for convenience (see code below)
    SI_rev[i] = SI[SI_cutoff - i + 1];       
  
  for(i in 1:SI_cutoff) {
    death_delay_rev[i] = death_delay[SI_cutoff - i + 1];
  }
}


parameters {
  real<lower=0> mu;                    // R0 - basic reproduction number
  real<lower=0> kappa;                 // Relating to the prior on the reproduction number
  real<lower=0> initial_infections;    // Initial infections for first N0 days
  real<lower=0> tau;                   // Prior on dispersion for initial_infections
  real<lower=0.001> phi;               // Overdispersion in the deaths data observational model
  real<lower=0,upper=1> ifr;           // Infection fatality ratio
  vector[W+1] weekly_effect;           // Time-varying Rt weekly effect
  real<lower=0, upper=1> weekly_rho;   // Relating to variance of weekly effects for time-varying Rt
  real<lower=0, upper=1> weekly_rho1;  // Relating to variance of weekly effects for time-varying Rt
  real<lower=0> weekly_sd;             // Relating to variance of weekly effects for time-varying Rt
}

transformed parameters {
  
  // Setting Up the Storage for Inferred Quantities
  vector<lower=0>[N2] infections = rep_vector(0, N2);   // Number of infections at each timepoint
  vector<lower=0>[N2] E_deaths  = rep_vector(0, N2);    // Number of deaths at each timepoint
  vector[N2] Rt = rep_vector(0, N2);                    // Reproduction number over time
  vector<lower=0>[N2] Rt_adj = Rt;
  {
    
    infections[1:N0] = rep_vector(initial_infections, N0); // learn the number of infections in the first N0 days and populate infections vector
    vector[N2] cumm_sum = rep_vector(0, N2);               // cumulative incidence of infections
    cumm_sum[2:N0] = cumulative_sum(infections[2:N0]);     // --> required for population immunity adjustment (susceptible depletion)

    // Calculating Time-Varying Rt
    for (i in 1:N2) {
      Rt[i] = mu * 2 * inv_logit(-weekly_effect[week_index[i]]);
    }

    // Calculating infections over time
    for (i in (N0 + 1):N2) {
      real convolution = 0; // UNCLEAR ABOUT THIS BEING RIGHT
      if ((i - SI_cutoff - 1) < 1) {
        convolution = dot_product(infections[1:(i - 1)], tail(SI_rev, i - 1)); 
      } else {
        convolution = dot_product(infections[(i - SI_cutoff):(i - 1)], tail(SI_rev, SI_cutoff)); 
      }
      cumm_sum[i] = cumm_sum[i - 1] + infections[i - 1];
      Rt_adj[i] = ((pop - cumm_sum[i]) / pop) * Rt[i];
      infections[i] = Rt_adj[i] * convolution;
    }
    
    // Calculating deaths over time
    E_deaths[1] = 1e-15 * infections[1]; // small amount of noise to avoid 0s at very beginning of time-series. 
    for (i in 2:N2){
      if ((i - SI_cutoff - 1) < 1) {
        E_deaths[i] = ifr * dot_product(infections[1:(i - 1)], tail(death_delay_rev, i - 1));
      } else {
        E_deaths[i] = ifr * dot_product(infections[(i - SI_cutoff):(i - 1)], tail(death_delay_rev, SI_cutoff));
      }
    }
  }
}

model {
  
  // Priors Relating to Inferred Initial Infections
  tau ~ exponential(0.1);
  initial_infections ~ exponential(1 / tau);
  
  // Priors Relating to the Basic Reproduction Number R0
  kappa ~ normal(0, 0.8);
  mu ~ normal(2, kappa); 
  
  // Priors Relating to Time-Varying Rt
  weekly_sd ~ normal(0,0.2);
  weekly_rho ~ normal(0.8, 0.05);
  weekly_rho1 ~ normal(0.1, 0.05);
  weekly_effect[1] ~ normal(0, 0.01);
  weekly_effect[2] ~ normal(0, weekly_sd * sqrt(1 - pow(weekly_rho, 2) - pow(weekly_rho1, 2) - 2 * pow(weekly_rho, 2) * weekly_rho1/(1 - weekly_rho1)));
  weekly_effect[3:(W+1)] ~ normal(weekly_effect[2:W]* weekly_rho + weekly_effect[1:(W-1)] * weekly_rho1, 
                                  weekly_sd *sqrt(1 - pow(weekly_rho, 2) - pow(weekly_rho1, 2) - 2 * pow(weekly_rho,2) * weekly_rho1/(1-weekly_rho1)));
  
  // Priors Relating to Deaths Observation Process
  ifr ~ normal(ifr_mean, ifr_sd);
  phi ~ normal(0, 5);
  deaths[Death_Fitting_Start:N2] ~ neg_binomial_2(E_deaths[Death_Fitting_Start:N2], phi);
}
