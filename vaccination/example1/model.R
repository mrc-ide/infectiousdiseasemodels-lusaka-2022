# code in the model equations
deriv(S) <- 
deriv(I) <- 
deriv(R) <- 

# initial conditions of the variables
initial(S) <- N - I0 - N * p
initial(I) <- I0
initial(R) <- N * p

# input parameter values
R0 <- user(4)
p <- user(0.25)

# calculated parameter values
mu <- b                     # assume constant population
beta <- R0 * (gamma + mu)   # transmission rate

# fixed parameters
I0 <- 1
b <- 1 / 52
gamma <- 1 / 10
N <- 1e6

# additional model outputs
output(Reff) <- (R0 * S) / (S + I + R)
