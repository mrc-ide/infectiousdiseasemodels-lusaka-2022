
## 0. Prelude, load packages and data ----
# make sure you have these necessary packages 
packages <- c("EpiEstim", "epitrix", "gridExtra", "ggplot2")
lapply(packages, require, character.only = TRUE)

# set working directory 
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# source helper_functions
source("helper_functions.R")

# load raw data (https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data)
cases <- read.csv("global_cases.csv", stringsAsFactors = FALSE)
deaths <- read.csv("global_deaths.csv", stringsAsFactors = FALSE)

#----

## 1. Process raw data and visualise as time series ----
cases <- process_data(cases, "Zambia")
deaths <- process_data(deaths, "Zambia")

# Visualise data
par(mfrow = c(2, 1))
plot_inc_cum(cases, "Cases")
plot_inc_cum(deaths, "Deaths")
par(mfrow = c(1, 1))
#----

## 2. Visualise exponential growth in cases data to use for R_0 calculations ----

# visualise data for first wave and identify the period of exponential growth
x <- cases[cases$date <= as.Date("2020-12-15"), ]
index <- min(which(x$cumulative > 0))
x <- x[index:nrow(x), ]
x$day <- seq(0, nrow(x) - 1, by = 1)

plot_inc_cum(x, "Cases", log = TRUE)

xx <- x[x$date <= as.Date("2020-06-15"), ]
xx$day <- seq(0, nrow(xx) - 1, by = 1)
plot_inc_cum(xx, "Cases", log = TRUE)

### As you can see from the plots above, daily cases data is very noisy! It is
### difficult to pin down a period of exponential growth. Let's aggregate data
### to weekly cases incidence!

# aggregate weekly
dt <- 7
tmp_wkly <- diff(cases$cumulative, lag = dt)
weekly <- data.frame(date = cases$date[7 + seq(1, length(tmp_wkly), dt)],
                     incidence = tmp_wkly[seq(1, length(tmp_wkly), dt)])
weekly$cumulative <- cumsum(weekly$incidence)

# check the dates align; here we defined weekly$cases as the number of cases
# in the week finishing on the date given in weekly$date
all(weekly$cumulative == cases$cumulative[cases$date %in% weekly$date])

plot_inc_cum(weekly, "Cases", log = TRUE)

### Data is now easier to visualise! Let's do the same to the data subsets for
### the first wave in 2020

x_weekly <- weekly[weekly$date <= as.Date("2020-12-15"), ]
index <- min(which(x_weekly$cumulative > 0))
x_weekly <- x_weekly[index:nrow(x_weekly), ]
x_weekly$day <- seq(0, nrow(x_weekly) - 1, by = 1)

plot_inc_cum(x_weekly, "Cases", log = TRUE)

xx_weekly <- weekly[weekly$date <= as.Date("2020-06-15"), ]
xx_weekly$day <- seq(0, nrow(xx_weekly) - 1, by = 1)
plot_inc_cum(xx_weekly, "Cases", log = TRUE)

### Hurray! During the first wave, particularly between mid-March and mid-June 
### 2020 we can clearly see a period where the epidemic grew exponentially.
### We can use these data to get estimates of R_0

#----

## 3. Estimate R0 from the data above ----

### As you may remember, to estimate R0 we need to get parameters or make 
### assumptions of the serial interval (SI)

# The mean serial interval (V) of Covid-19 ranges from 4.2 to 7.5 days
# with RE meta-analysis from 23 studies 5.2 (95%CI 4.9 - 5.5) and
# incubation period from 14 studies 6.5 (95%CI 4.9 5.5)
# (https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-021-05950-x).
# High CT values are on average lowest from 5 days pre- to 10 days post-symptoms
# (https://pubmed.ncbi.nlm.nih.gov/33160066/)
# Another study of infector/infectee pairs suggested a mean and STD of the 
# SI of 4.7 and 2.9 days; here we will use the mid-point, 3.8 days.
# (https://www.ijidonline.com/article/S1201-9712(20)30119-3/pdf)
# Additionally we will calculate the STD of the mean estimate from the first
# study (meta-analysis) of 23 studies, which we will use later for EpiEstim.
V <- c(4.9, 5.5)
mean_si <- 5.2
std_si <- 3.8
std_mean_si <- sqrt(23) * (V[2] - V[1]) / 3.92


# Regress cumulative incidence over time (remember our data is weekly, so each
# row represents a week)
# split incidence data by 5 day intervals (mean SI)
dt <- 7
time_step <- rep(1:(nrow(xx) / dt), each = dt)
time_step <- time_step[1:nrow(xx)]
y <- sapply(split(xx$incidence, time_step), sum)

y <- data.frame(
  time = as.numeric(names(y)),
  inc = unname(y),
  cum = cumsum(y)
)

# N.B. here we are doing log(inc + 1), given sometimes you have zeroes, which 
# will return NAs if calculating their log (very useful trick with noisy data!)
log <- lm(data = y, log(inc + 1) ~ time)

# Rate of exponential growth
r_weekly <- unname(log$coefficients["time"])
# remember we are using weekly time steps, so we need to divide
# by 7 here to transform the growth rate estimate back to days
r <- r_weekly / dt
r

# doubling time
log(2) / r

# visualise the fit:
plot(y$time, log$coefficients[["(Intercept)"]] + r_weekly * y$time,
     col = "blue", lty = 2)

# A 'ball-park' figure for R_0 using only the mean of the SI
R0_central <- mean_si * r + 1
R0_central

# A range for the R0_central, assuming a constant SI
R0_si_const <- V * r + 1
R0_si_const

# An improved R0 calculation, providing the whole distribution of the SI (w)
w <- EpiEstim::discr_si(0:20, mean_si, std_si)
R0_si_dist <- epitrix::r2R0(r, w)
R0_si_dist

# So we have learnt a lot here!
# 1. We have used Zambian cases data from the first wave to identify a period of
# expoenential growth and calculated that rate of exponential growth, r.
# 2. We investigated some values of the mean and std of the serial interval
# from the literature. 
# 2. We used our knowledge of the growth rate, r, and the SI to generate some
# estimates of R0 for the Zambian epidemic!

#----

## 4. Compare results to EpiEstim with different methods ----

## A. Estimating parametric SI using mean and std from literature

# The default config is 7-day sliding window for R estimates
# let's change this since to fortnightly estimates as data is quite noisy
window_width <- 14
t_start <- seq(2, nrow(xx) - (window_width - 1))
t_end <- t_start + (window_width - 1)

estim <- EpiEstim::estimate_R(xx$incidence,
                              method = "parametric_si",
                              config = EpiEstim::make_config(list(
                                mean_si = mean_si,
                                std_si = std_si, 
                                t_start = t_start,
                                t_end = t_end)))
plot(estim, "R")

# now let's explore the range for the serial interval
lb_estim <- EpiEstim::estimate_R(xx$incidence,
                                 method = "parametric_si",
                                 config = EpiEstim::make_config(list(
                                   mean_si = V[1],
                                   std_si = std_si,
                                   t_start = t_start,
                                   t_end = t_end)))

ub_estim <- EpiEstim::estimate_R(xx$incidence,
                                 method = "parametric_si",
                                 config = EpiEstim::make_config(list(
                                   mean_si = V[2],
                                   std_si = std_si,
                                   t_start = t_start,
                                   t_end = t_end)))

plot(estim$R$`Mean(R)`, type = "l", bty = "n",
     ylab = expression(R[t]), ylim = c(0, 5))
lines(lb_estim$R$`Quantile.0.025(R)`, col = grey(0.6))
lines(ub_estim$R$`Quantile.0.975(R)`, col = grey(0.6))
abline(h = 1, lty = 2, col = "red")

### So here we have it! We have estimated Rt given a parametric SI, which we have
### inferred from SI mean and std from the literature
### Given EpiEstim, this is the estimated value of Rt (~R0) at the start of the epidemic:
c("mean" = estim$R$`Mean(R)`[1], "lb" = lb_estim$R$`Quantile.0.025(R)`[1],
  "ub" = ub_estim$R$`Quantile.0.975(R)`[1])



## B. Estimating parametric SI when there is uncertainty (i.e. no literature to
##     inform our thinkin)

# You will not always have good literature to infor your thinking on the mean
# and std of the SI. 
# In such a case, you can actually include uncertainty in the SI directly
# in EpiEstim
# EpiEstim vignette that Rebecca wrote? 
uncertain_estim <- EpiEstim::estimate_R(xx$incidence,
                                        method = "uncertain_si",
                                        config = EpiEstim::make_config(list(
                                          mean_si = mean_si,
                                          std_mean_si = std_mean_si,
                                          min_mean_si = V[1],
                                          max_mean_si = V[2],
                                          std_si = std_si,
                                          std_std_si = 0.001, ### just making this up for now to be very small in case all you want is to explore uncertainty in the mean
                                          min_std_si = std_si - 0.001, ### just making this up for now to be very small in case all you want is to explore uncertainty in the mean
                                          max_std_si = std_si + 0.001, ### just making this up for now to be very small in case all you want is to explore uncertainty in the mean
                                          t_start = t_start,
                                          t_end = t_end,
                                          n1 = 100,
                                          n2 = 100)))

# You can compare with the estimate with no uncertainty - there is very little 
# difference in the mean but there should be more uncertainty in the latter
# this is not very visible here because once discretised the distributions
# of the SI are actually quite similar, you can see this here:
gridExtra::grid.arrange(plot(estim, "R") + ggplot2::labs(title = "No uncertainty"),
                        plot(uncertain_estim, "R") + ggplot2::labs(title = "With uncertainty"))
# and here you can see the uncertainty in the SI distribution
matplot(t(uncertain_estim$si_distr), type = "l", lty = 1, xlim = c(0, 15), 
        col = scales::alpha("black", 0.1), 
        xlab = "Days", ylab = "Positive matrix factorisation")

#----


###### Well done!!! ######
## You've learnt a lot about the likely value of R0 and Rt for Covid-19 in Zambia!
## We will need these knowledge for building and parameterising our model.
## Keep track of the values for R0 you have estimated. Compare this initial value
## to how Rt has changed over the course of the epidemic until 2020-06-15. 
## This knowledge will allow you to fit your model to the epidemic data.
