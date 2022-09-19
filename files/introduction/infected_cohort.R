
# state variables
deriv(I) <- -gamma * I
deriv(R) <- gamma * I

# initial conditions of the variables
initial(I) <- 1000
initial(R) <- 0

# input parameters
recovery_time <- 10     # mean number of days to recovery

# calculated parameters
gamma <- 1 / recovery_time    # recovery rate
