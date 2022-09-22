
# set working directory 
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# source helper_functions
source("helper_functions.r")

# load raw data (https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data)
cases <- read.csv("global_cases.csv", stringsAsFactors = FALSE)
deaths <- read.csv("global_deaths.csv", stringsAsFactors = FALSE)

## a. Process raw data into cumulative and incidence time series ----
cases <- process_data(cases, "Zambia")
deaths <- process_data(deaths, "Zambia")

# Visualise data
par(mfrow = c(2, 1))
plot_inc_cum(cases, "Cases")
plot_inc_cum(deaths, "Deaths")
dev.off()
#----

## b. Calculate R0 ----

# visualise data for first wave and identify the period of exponential growth
x <- cases[cases$date <= as.Date("2020-12-15"), ]
index <- min(which(x$cumulative > 0))
x <- x[index:nrow(x), ]
x$day <- seq(0, nrow(x) - 1, by = 1)

plot_inc_cum(x, "Cases")

xx <- x[x$date <= as.Date("2020-06-15"), ]
xx$day <- seq(0, nrow(xx) - 1, by = 1)
plot_inc_cum(xx, "Cases")

# The mean serial interval (V) of Covid-19 ranges from 4.2 to 7.5 days
# with RE meta-analysis from 23 studies 5.2 (95%CI 4.9 - 5.5) and
# incubation period from 14 studies 6.5 (95%CI 4.9 5.5)
# (https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-021-05950-x).
# High CT values are on average lowest from 5 days pre- to 10 days post-symptoms
# (https://pubmed.ncbi.nlm.nih.gov/33160066/)
V <- c(4.9, 5.5)
mean_si <- 5.2
std_si <- sqrt(23) * (V[2] - V[1]) / 3.92

# Regress cumulative incidence over time
# split incidence data by 5 day intervals (mean SI)
time_step <- rep(1:(nrow(xx) / 5), each = 5)
y <- sapply(split(xx$incidence, time_step), sum)

y <- data.frame(
  time = as.numeric(names(y)),
  inc = unname(y),
  cum = cumsum(y)
)

log <- lm(data = y, log(cum) ~ time)

# Rate of exponential growth
r <- log$coefficients["time"]
# doubling time
log(2) / r

# R_0
R0 <- V * r + 1
R0

#----


## Compare results to EpiEstim with different methods ----

# parametric SI using mean and std from literature
# default config is 7-day sliding window for R estimates
# let's change this since to fortnightly estimates as data is quite noisy
t_start <- seq(2, nrow(xx) - 13)
t_end <- t_start + 13

estim <- EpiEstim::estimate_R(xx$incidence,
                                  method = "parametric_si",
                                  config = EpiEstim::make_config(list(
                                    mean_si = mean_si,
                                    std_si = std_si,
                                    t_start = t_start,
                                    t_end = t_end
                                    )))
plot(estim, "R")

# now let's explore the range for the serial interval
lb_estim <- EpiEstim::estimate_R(xx$incidence,
                                      method = "parametric_si",
                                      config = EpiEstim::make_config(list(
                                        mean_si = V[1],
                                        std_si = std_si,
                                        t_start = t_start,
                                        t_end = t_end
                                      )))

ub_estim <- EpiEstim::estimate_R(xx$incidence,
                                      method = "parametric_si",
                                      config = EpiEstim::make_config(list(
                                        mean_si = V[2],
                                        std_si = std_si,
                                        t_start = t_start,
                                        t_end = t_end
                                      )))

plot(estim$R$`Mean(R)`, type = "l", bty = "n",
     ylab = expression(R[t]), ylim = c(0, 5))
lines(lb_estim$R$`Quantile.0.025(R)`, col = grey(0.6))
lines(ub_estim$R$`Quantile.0.975(R)`, col = grey(0.6))
abline(h = 1, lty = 2, col = "red")

# now let's define discrete time intervals to calculate Rt which we can use for
# our model
t_start <- seq(2, nrow(xx) - 21, by = 21)
t_end <- t_start + 21

estim_discrete <- EpiEstim::estimate_R(xx$incidence,
                                       method = "parametric_si",
                                       config = EpiEstim::make_config(list(
                                         mean_si = mean_si,
                                         std_si = std_si,
                                         t_start = t_start,
                                         t_end = t_end
                                       )))
plot(estim_discrete, "R")

betas <- estim_discrete$R$`Mean(R)`

round((betas * R0[2]) / max(betas), 1)

#----



## c. Estimate the force of infection (beta) ----

# We assume that: a) cases over time follow a Poisson distribution; b) infections
# at time (t+1) are the result of a binomial draw; c) contacts happen at random;
# and d) the serial interval is the basic unit of time.

# give candidate values for beta and estimate the MLE
S0_candidate <- 18.4e6
beta_candidate <- seq(0, 5, by = 0.1)

ll_beta <- rep(NA, length(beta_candidate))

for (i in 1:length(beta_candidate)) {
  
  ll_beta[i] <- ll_infected(S0 = S0_candidate, beta = beta_candidate[i], I = y$inc)
    
}

mle_beta <- beta_candidate[which.min(ll_beta)]

plot(ll_beta ~ beta_candidate,
     ylab = "Negative log-likelihood", xlab = expression(beta),
     bty = "n")
abline(v = mle_beta, col = "red", lty = 2)
legend("right", bty = "n", legend = paste("MLE =", mle_beta),
       lty = 2, col = "red")


# give candidate values for S0 and estimate the MLE
S0_candidate <- seq(1.4e6, 16e6, length = 20)
beta_candidate <- mle_beta

ll_s0 <- rep(NA, length(S0_candidate))

for (i in 1:length(S0_candidate)) {
  
  ll_s0[i] <- ll_infected(S0 = S0_candidate[i],
                          beta = beta_candidate,
                          I = y$inc)
}

mle_S0 <- S0_candidate[which.min(ll_s0)]

plot(ll_s0 ~ S0_candidate,
     ylab = "Negative log-likelihood", xlab = expression(S[0]),
     bty = "n")
abline(v = mle_S0, col = "red", lty = 2)
legend("right", bty = "n", legend = paste("MLE =", mle_S0),
       lty = 2, col = "red")

#----


## MLE for parameters in a chain-binomial model ----

fit <- bbmle::mle2(ll_infected,
                   start = list(S0 = mle_S0, beta = mle_beta),
                   data = list(I = y$inc),
                   method = "Nelder-Mead")
bbmle::summary(fit)
cov2cor(bbmle::vcov(fit))

#----

## Stochastic simulation ----
sim_covid <- function(S0, beta, I0) {
  
  I <- I0
  S <- S0
  i = 1
  
  while (!any(I == 0)) {
    i <- i + 1
    I[i] <- rbinom(1, size = S[i - 1], prob = 1 - exp(-beta * I[i - 1] / S0))
    S[i] <- S[i - 1] - I[i]
  }
  
  data.frame(S = S, I = I)
  
}

plot(y$inc, type = "n", xlim = c(1, 18), bty = "n")
legend("topleft", legend = c("Data", "Model"), lty = c(1, 1), bty = "n",
       pch = c(1, NA), col = c("red", "grey"))
for (i in 1:100) {
  sim <- sim_covid(S0 = floor(bbmle::coef(fit)["S0"]),
                   beta = bbmle::coef(fit)["beta"],
                   I0 = 1)
  
  lines(sim$I, col = grey(0.7))
}
points(y, type = "b", col = "red")

#----
