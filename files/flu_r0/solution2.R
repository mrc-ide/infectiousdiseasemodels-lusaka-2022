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
R_0_closure <- user(1.5)
D <- user(1)
I_0 <- 1 
N <- 370

# convert parameters
gamma <- 1 / D
beta <- if (t > 18 && t < 25) R_0_closure * gamma else R_0 * gamma

