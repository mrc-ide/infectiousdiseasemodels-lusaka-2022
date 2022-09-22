
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
