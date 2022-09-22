# variables
deriv(S) <- -beta * S * I / N
deriv(I) <- beta * S * I / N - gamma * I
deriv(R) <- gamma * I
deriv(V) <- 0

# initial conditions of the variables
initial(S) <- pop * (1 - p_vacc) - I0
initial(I) <- I0
initial(R) <- 0
initial(V) <- p_vacc * pop

# input parameter values
I0 <- user(1)
days_to_onset <- user(4)
days_to_recovery <- user(7)
p_vacc <- user(0)

# calculated parameter values
N <- S + I + V + R
beta <- 1 / days_to_onset
gamma <- 1 / days_to_recovery

# constant parameters
pop <- 1e4
