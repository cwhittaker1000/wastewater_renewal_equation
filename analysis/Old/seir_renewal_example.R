# Load required libraries
library(odin)

# Loading SIR model
seir <- odin::odin({
                   
    ## Core equations for transitions between compartments:
    deriv(S) <- - beta * S * I / N
    deriv(E) <- if(t < num_days) init_infections + beta * S * I / N - gamma * E else beta * S * I / N - gamma * E 
    deriv(I) <- gamma * E - sigma * I
    deriv(R) <- sigma * I
    
    ## Time-varying beta
    Rt <- interpolate(tt_Rt, Rt_set, "constant")
    tt_Rt[] <- user()
    Rt_set[] <- user()
    dim(tt_Rt) <- user()
    dim(Rt_set) <- length(tt_Rt)
    beta <- Rt #* sigma * N / S # from unadjusted Rt create beta that gives that Rt irrespective of population-level immunity
    
    ## Total population size (odin will recompute this at each timestep)
    N <- S + E + I + R
    
    ## Initial states:
    initial(S) <- S_ini 
    initial(E) <- E_ini 
    initial(I) <- I_ini
    initial(R) <- 0
    
    ## User defined parameters - default in parentheses:
    S_ini <- user(1000)
    E_ini <- user(1)
    I_ini <- user(1)
    gamma <- user(0.1)
    sigma <- user(0.1)
    init_infections <- user(10)
    num_days <- user(10)
    
    incidence <- beta * S * I / N
    output(incidence) <- TRUE
    
    incidence2 <- if(t < num_days) init_infections + beta * S * I / N else beta * S * I / N
    output(incidence2) <- TRUE
    
    Rt_no_adj <- beta / sigma
    output(Rt_no_adj) <- TRUE
    Rt_adj <- (beta / sigma) * (S / N)
    output(Rt_adj) <- TRUE
    output(beta) <- TRUE
    output(N) <- TRUE
                   
})

dur <- 200
pop <- 10^6
seir_mod <- seir$new(Rt_set = c(0.15, 0.12, 0.05),
                     tt_Rt = c(0, 60, 100),
                     init_infections = 5,
                     gamma = 0.1, 
                     sigma = 0.1, 
                     S_ini = pop, 
                     E_ini = 5, 
                     I_ini = 0)
output <- seir_mod$run(1:dur, rtol = 1e-12, atol = 1e-12)
sir_col <- c("#8c8cd9", "#EDBF85", "#cc0044", "#83C166", "black")
par(mfrow = c(1, 3))
matplot(output[, 1], output[, "Rt_no_adj"], xlab = "Time", ylab = "Time-Varying Rt (Unadjusted for Pop Immunity)", type = "l", col = sir_col, lty = 1, ylim = c(0, max(output[, "Rt_no_adj"])))
matplot(output[, 1], output[, "Rt_adj"], xlab = "Time", ylab = "Time-Varying Rt (Adjusted for Pop Immunity)", type = "l", col = sir_col, lty = 1, ylim = c(0, max(output[, "Rt_adj"])))
matplot(output[, 1], output[, "beta"], xlab = "Time", lty = 1, type = "l", xlim = c(0, 100), ylim = c(0, 1))

par(mfrow = c(1, 1))
matplot(output[, 1], output[, 2:5], xlab = "Time", ylab = "Number of individuals", type = "l", col = sir_col, lty = 1)
legend("topright", lwd = 1, col = sir_col, legend = c("S", "I", "R", "R", "incidence"), bty = "n")

# Generating data for model fitting
set.seed(11)
noisy_data <- rnbinom(n = dur, mu = c(0, diff(output[1:dur, "R"])), size = 10^5) # low size = more overdispersion
matplot(output[1:dur, 1], c(0, diff(output[1:dur, "R"])), xlab = "Time", ylab = "Number of individuals", type = "l", col = sir_col[3], lty = 1, ylim = c(0, max(noisy_data)))
matpoints(output[1:dur, 1], noisy_data, xlab = "Time", ylab = "Number of individuals", col = sir_col[3], pch = 20)

# Generating data list
week_index <- c(rep(1, seir_mod$contents()$tt_Rt[2] - seir_mod$contents()$tt_Rt[1]),
                rep(2, seir_mod$contents()$tt_Rt[3] - seir_mod$contents()$tt_Rt[2]),
                rep(3, dur - seir_mod$contents()$tt_Rt[3]))
data_list <- list(N0 = 10, 
                  N2 = dur,
                  deaths = noisy_data,
                  Death_Fitting_Start = 10,
                  week_index = week_index, 
                  W = max(week_index),
                  SI_cutoff = 75,
                  pop = max(output[, "N"]))
infection_SI <- dgamma(x = 1:data_list$N2, shape = 2, scale = 1 / seir_mod$contents()$gamma)
death_SI <- dgamma(x = 1:data_list$N2, shape = 2, scale = 1 / seir_mod$contents()$gamma)
data_list$SI <- infection_SI # [1:data_list$SI_cutoff]
data_list$death_delay <- death_SI #[1:data_list$SI_cutoff]

# try to work out what's going on with the weird infections values (negatives)
mod <- cmdstanr::cmdstan_model(stan_file = "models/base_renewal_equation2.stan")
fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 100,
  list(
    list(mu = 2,
         kappa = 0.5,
         initial_infections = 2,
         tau = 5,
         phi = 5))
)
draws_df <- fit$draws(format = "df")
fit$summary()
par(mfrow = c(3, 2))
hist(draws_df$mu)
hist(draws_df$kappa)
hist(draws_df$ifr)
hist(draws_df$phi)
hist(draws_df$initial_infections)
hist(draws_df$tau)

# Infections Over Time
par(mfrow = c(2, 2))
infections_output <- draws_df[, grep("infections0*", colnames(draws_df))[1:dur]]
infections <- unname(colMeans(infections_output))
infections_lower <- unname(apply(infections_output, 2, quantile, 0.05))
infections_upper <- unname(apply(infections_output, 2, quantile, 0.95))

matplot(output[1:dur, 1], output[1:dur, "incidence2"], xlab = "Time", ylab = "Number of infections", pch = 20, col = sir_col[3], 
        lty = 1, 
        ylim = c(0, max(c(output[1:dur, "incidence2"], infections_upper[1:dur]))))
matlines(output[1:dur, 1], infections[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], infections_lower[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], infections_upper[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)

# Recoveries Over Time
recoveries_output <- draws_df[, grep("deaths0*", colnames(draws_df))[1:dur]]
recoveries <- unname(colMeans(recoveries_output))
recoveries_lower <- unname(apply(recoveries_output, 2, quantile, 0.05))
recoveries_upper <- unname(apply(recoveries_output, 2, quantile, 0.95))

matplot(output[1:dur, 1], noisy_data[1:dur], xlab = "Time", ylab = "Number of recoveries", col = sir_col[3], pch = 20,
        ylim = c(0, max(c(noisy_data[1:dur], recoveries[1:dur]))))
matlines(output[1:dur, 1], recoveries[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], recoveries_lower[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)
matlines(output[1:dur, 1], recoveries_upper[1:dur], xlab = "Time", ylab = "Number of individuals", col = "black", pch = 20)

Rt <- draws_df[grep("Rt\\[", colnames(draws_df))]
plot(apply(Rt, 2, quantile, 0.95)[1:dur], type = "l", ylab = "Reff", ylim = c(0, max(apply(Rt, 2, quantile, 0.95))))
lines(apply(Rt, 2, quantile, 0.05)[1:dur], type = "l", ylab = "Reff")
lines(colMeans(Rt)[1:dur], type = "l", ylab = "Reff")

Rt_adj <- draws_df[grep("Rt_adj\\[", colnames(draws_df))]
plot(apply(Rt_adj, 2, quantile, 0.95)[1:dur], type = "l", ylab = "Reff", ylim = c(0, max(apply(Rt_adj, 2, quantile, 0.95))))
lines(apply(Rt_adj, 2, quantile, 0.05)[1:dur], type = "l", ylab = "Reff")
lines(colMeans(Rt_adj)[1:dur], type = "l", ylab = "Reff")

