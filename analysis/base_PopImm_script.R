## Tester Example for Renewal Equation Model WITH Population Immunity

### Inputs for Generating Synthetic Data
duration <- 350
initial_infections <- 10
SI_cutoff <- 200
N0 <- 10
N2 <- 350
ifr <- 0.01
infection_SI <- dgamma(x = 1:N2, shape = 2, scale = 1 / 0.1)
infection_SI_rev <- rev(infection_SI)
death_SI <- dgamma(x = 1:N2, shape = 2, scale = 1 / 0.1)
death_SI_rev <- rev(death_SI)
Death_Fitting_Start <- 10
Rt_synth <- c(rep(2.2, 150), rep(1.2, 100), rep(1.3, 50), rep(1.6, 50))
Rt_adj_synth <- vector(mode = "numeric", length = N2)
pop <- 10^6

### Generating Synthetic Data 
infections <- vector(mode = "numeric", length = N2)
infections[1:N0] <- initial_infections
deaths <- vector(mode = "numeric", length = N2)
cumm_sum <- vector(mode = "numeric", length = N2)
for (i in N0:N2) {
  if ((i - SI_cutoff - 1) < 1) {
    convolution <- sum(infections[1:(i - 1)] * tail(infection_SI_rev, i - 1)) 
  } else {
    convolution <- sum(infections[(i - SI_cutoff):(i - 1)] * tail(infection_SI_rev, SI_cutoff)) 
  }
  cumm_sum[i] <- cumm_sum[i - 1] + infections[i - 1]
  Rt_adj_synth[i] <- ((pop - cumm_sum[i]) / pop) * Rt_synth[i]
  infections[i] <- Rt_adj_synth[i] * convolution
}
for (i in 2:N2) {
  if ((i - SI_cutoff - 1) < 1) {
    deaths[i] <- ifr * sum(infections[1:(i - 1)] * tail(death_SI_rev, i - 1)) 
  } else {
    deaths[i] <- ifr * sum(infections[(i - SI_cutoff):(i - 1)] * tail(death_SI_rev, SI_cutoff)) 
  }
}

# Plotting Synthetic Data to Check Looks Reasonable
par(mfrow = c(1, 4))
plot(Rt_synth, type = "l", ylim = c(0, max(Rt_synth)), main = "Rt (No Adj)", ylab = "Rt")
plot(Rt_adj_synth, type = "l", ylim = c(0, max(Rt_adj_synth)), main = "Rt (Adj)", ylab = "Rt")
plot(infections, type = "l", main = "Infections")
set.seed(11)
noisy_data <- rnbinom(n = N2, mu = deaths, size = 50) # low size = more overdispersion
matplot(deaths, xlab = "Time", ylab = "Number of deaths", type = "l", 
        col = "blue", lty = 1, ylim = c(0, max(noisy_data)), main = "Noisy Deaths")
matpoints(noisy_data, col = "blue", pch = 20)

### Generating Data List for Input to Stan
week_index <- c(rep(1, sum(Rt_synth == unique(Rt_synth)[1])), rep(2, sum(Rt_synth == unique(Rt_synth)[2])),
                rep(3, sum(Rt_synth == unique(Rt_synth)[3])), rep(4, sum(Rt_synth == unique(Rt_synth)[3])))
data_list <- list(N0 = N0, 
                  N2 = N2,
                  ifr_mean = 0.02,
                  ifr_sd = 0.02,
                  deaths = noisy_data,
                  Death_Fitting_Start = 10,
                  week_index = week_index, 
                  W = max(week_index),
                  SI_cutoff = 100,
                  SI = infection_SI,
                  death_delay = death_SI,
                  pop = pop)

### Compiling and Fitting Model to Synthetic Data
mod <- cmdstanr::cmdstan_model(stan_file = "models/base_PopImm.stan")
fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 100,
  list(list(mu = 2, kappa = 0.5, initial_infections = 2, tau = 5, phi = 5)))

### Processing and Plotting Model Outputs
draws_df <- fit$draws(format = "df")
fit$summary()
layout(matrix(c(1,2,3,0,
                4,5,6,0,
                7,8,9,10), 3, 4, byrow = TRUE))

#### Histograms of Posterior Distribution for Model Parameters
hist(draws_df$mu, xlab = "", ylab = "", main = "Mu (R0)")
abline(v = Rt_synth[1], col = "red", lwd = 5)
hist(draws_df$kappa, xlab = "", ylab = "", main = "Kappa")
hist(draws_df$ifr, xlab = "", ylab = "", main = "IFR")
abline(v = ifr, col = "red", lwd = 5)
hist(draws_df$phi, xlab = "", ylab = "", main = "Phi")
hist(draws_df$initial_infections, xlab = "", ylab = "", main = "Initial Infections")
abline(v = initial_infections, col = "red", lwd = 5)
hist(draws_df$tau, xlab = "", ylab = "", main = "tau")

#### Rt Over Time and Comparison to Synthetic Data
Rt <- draws_df[grep("Rt\\[", colnames(draws_df))]
plot(apply(Rt, 2, quantile, 0.95)[1:N2], type = "l", ylab = "Rt", 
     ylim = c(0, max(apply(Rt, 2, quantile, 0.95))), col = "red", main = "Rt")
lines(apply(Rt, 2, quantile, 0.05)[1:N2], type = "l", ylab = "Reff", col = "red")
lines(colMeans(Rt)[1:N2], type = "l", ylab = "Reff", col = "red")
lines(Rt_synth[1:N2], col = "black")

#### Rt Over Time and Comparison to Synthetic Data
Rt_adj <- draws_df[grep("Rt_adj\\[", colnames(draws_df))]
plot(apply(Rt_adj, 2, quantile, 0.95)[1:N2], type = "l", ylab = "Rt_adj", 
     ylim = c(0, max(apply(Rt_adj, 2, quantile, 0.95))), col = "red", main = "Rt_adj")
lines(apply(Rt_adj, 2, quantile, 0.05)[1:N2], type = "l", ylab = "Reff", col = "red")
lines(colMeans(Rt_adj)[1:N2], type = "l", ylab = "Reff", col = "red")
lines(Rt_adj_synth[1:N2], col = "black")

#### Infections Over Time and Comparison to Synthetic Data
infections_output <- draws_df[, grep("infections0*", colnames(draws_df))[1:N2]]
infections <- unname(colMeans(infections_output))
infections_lower <- unname(apply(infections_output, 2, quantile, 0.05))
infections_upper <- unname(apply(infections_output, 2, quantile, 0.95))
matplot(infections, xlab = "Time", ylab = "Number of infections", pch = 20, col = "black", 
        lty = 1, ylim = c(0, max(infections_upper[1:N2])), main = "Infections Over Time")
matlines(infections[1:N2], col = "red", pch = 20)
matlines(infections_lower[1:N2], col = "red", pch = 20)
matlines(infections_upper[1:N2], col = "red", pch = 20)

#### Deaths Over Time and Comparison to Synthetic Data
deaths_output <- draws_df[, grep("deaths0*", colnames(draws_df))[1:N2]]
deaths <- unname(colMeans(deaths_output))
deaths_lower <- unname(apply(deaths_output, 2, quantile, 0.05))
deaths_upper <- unname(apply(deaths_output, 2, quantile, 0.95))
matplot(noisy_data[1:N2], xlab = "Time", ylab = "Number of deaths", col = "black", pch = 20, 
        ylim = c(0, max(c(noisy_data[1:N2], deaths[1:N2]))), main = "Deaths Over Time")
matlines(deaths[1:N2], col = "red", pch = 20)
matlines(deaths_lower[1:N2], col = "red", pch = 20)
matlines(deaths_upper[1:N2], col = "red", pch = 20)
