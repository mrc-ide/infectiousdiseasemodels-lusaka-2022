# initial conditions
initial(S) <- N - I_0
initial(I) <- I_0
initial(R) <- 0

# equations
deriv(S) <- -beta * S * (I / N)
deriv(I) <- beta * (I / N) * S - gamma * I
deriv(R) <- gamma * I

# parameter values
R_0 <- user(1.5)
D <- user(1)
I_0 <- 1 # default value
N <- 370

# convert parameters
beta <- R_0 / D
gamma <- 1 / D
