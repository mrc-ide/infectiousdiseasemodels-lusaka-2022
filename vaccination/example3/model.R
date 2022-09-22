# variables
deriv(S) <- -beta * S * I / N
deriv(I) <- beta * S * I / N - gamma * I
deriv(R) <- gamma * I
deriv(V) <- 0

# initial conditions of the variables
initial(S) <- pop * (1 - p_vacc * ve) - I0
initial(I) <- I0
initial(R) <- 0
initial(V) <- pop * p_vacc * ve

# input parameter values
p_vacc <- user(0)
ve <- user(0.7)

# calculated parameter values
N <- S + I + V + R
beta <- 1 / days_to_onset
gamma <- 1 / days_to_recovery

# constant parameters
pop <- 1e4
I0 <- 1
days_to_onset <- 4
days_to_recovery <- 7
