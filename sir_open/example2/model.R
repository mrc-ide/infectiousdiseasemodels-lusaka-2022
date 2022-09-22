# age group 1
deriv(S_1) <- - lambda_1 * S_1 + sigma * R_1
deriv(I_1) <- lambda_1 * S_1 - gamma * I_1
deriv(R_1) <- gamma * I_1 - sigma * R_1

# age group 2
deriv(S_2) <- - lambda_2 * S_2 + sigma * R_2
deriv(I_2) <- lambda_2 * S_2 - gamma * I_2
deriv(R_2) <- gamma * I_2 - sigma * R_2

# age group 3
deriv(S_3) <- - lambda_3 * S_3 + sigma * R_3
deriv(I_3) <- lambda_3 * S_3 - gamma * I_3
deriv(R_3) <- gamma * I_3 - sigma * R_3


### Initial conditions ###

# age group 1
initial(S_1) <- N_1
initial(I_1) <- 0
initial(R_1) <- 0

# age group 2
initial(S_2) <- N_2 - I0
initial(I_2) <- I0
initial(R_2) <- 0

# age group 3
initial(S_3) <- N_3
initial(I_3) <- 0
initial(R_3) <- 0

# input parameters
N <- user(1e6)             # total population size
prop_age_1 <- user(0.2)    # population of age group 1
prop_age_2 <- user(0.65)   # population of age group 2
prop_age_3 <- user(0.15)   # population of age group 3
I0 <- user(1)              # num infectious cases at start of epidemic

b <- user(0.02)            # probability of transmission per contact
gamma <- user(0.07)      # recovery rate
sigma <- user(0.003)     # waning rate

c_1_1 <- user(7)
c_1_2 <- user(5)
c_1_3 <- user(1)

c_2_1 <- user(2)
c_2_2 <- user(9)
c_2_3 <- user(1)

c_3_1 <- user(1)
c_3_2 <- user(3)
c_3_3 <- user(2)


### Calculated parameters ###

# population by age
N_1 <- N * prop_age_1
N_2 <- N * prop_age_2
N_3 <- N * prop_age_3

# force of infection by age
lambda_1 <- (b * c_1_1 * I_1 / N_1) +
  (b * c_1_2 * I_2 / N_2) +
  (b * c_1_3 * I_3 / N_3)

lambda_2 <- (b * c_2_1 * I_1 / N_1) +
  (b * c_2_2 * I_2 / N_2) +
  (b * c_2_3 * I_3 / N_3)

lambda_3 <- (b * c_3_1 * I_1 / N_1) +
  (b * c_3_2 * I_2 / N_2) +
  (b * c_3_3 * I_3 / N_3)

# population by age over time
pop_1 <- S_1 + I_1 + R_1
pop_2 <- S_2 + I_2 + R_2
pop_3 <- S_3 + I_3 + R_3
pop <- pop_1 + pop_2 + pop_3

# additional model outputs
output(pop_1) <- TRUE
output(pop_2) <- TRUE
output(pop_3) <- TRUE

output(lambda_1) <- TRUE
output(lambda_2) <- TRUE
output(lambda_3) <- TRUE