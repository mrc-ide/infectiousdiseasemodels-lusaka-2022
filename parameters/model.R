# initial conditions
initial(I) <- 1e6
initial(R) <- 0
initial(D) <- 0

# equations
deriv(I) <- -gamma * I - mu * I
deriv(R) <- gamma * I
deriv(D) <- mu * I

# parameter values
recovery_time <- user()
death_time <- user()

# convert parameters
gamma <- 1 / recovery_time   # recovery rate
mu <- 1 / death_time         # death rate
