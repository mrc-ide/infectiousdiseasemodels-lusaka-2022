# variables
deriv(S) <-
deriv(I) <-
deriv(R) <-

# initial conditions
initial(S) <- pop - I0
initial(I) <- I0
initial(R) <- 0

# User parameters
R0 <- user(2.5)            # basic reproduction number
recov_t <- user(14)        # mean time from onset to recovery (DAYS)
mean_birth <- user(1)      # mean number of births (life-time)

# Fixed parameters
pop <- 1e6            # total population size
life_exp <- 64        # mean life expectancy (YEARS)
I0 <- 1               # Initial # infectious cases (seed)

# Calculated parameters
beta <- R0 * (gamma + mu)      # transmission rate
gamma <- 1 / (recov_t)         # recovery rate
mu <- 1 / (life_exp * 365)     # death date (match model time scale)
b <- mu * mean_birth           # birth rate

N <- S + I + R               # total population per unit of time

# Additional model outputs
output(Rt) <- (R0 * S) / N    # reproduction number per unit of time
output(N) <- TRUE
output(R0) <- TRUE
