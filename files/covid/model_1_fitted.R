###############################
###     Model dynamics      ###
###############################

deriv(S) <- - S/N * (beta_r * I_r + beta_d * I_d)

deriv(E) <- S/N * (beta_r * I_r + beta_d * I_d) - gamma * E

deriv(I_d) <- ifr * gamma * E - sigma_d * I_d

deriv(I_r) <- (1 - ifr) * gamma * E - sigma_r * I_r

deriv(R) <- sigma_r * I_r

deriv(Dead) <- sigma_d * I_d

deriv(cumul_onset) <- gamma * E

cumul_death <- Dead

### create delayed variables in order to compute non cumulative (weekly) 
### incidence and deaths variables
cumul_onset_delayed <- delay(cumul_onset, 7) # 7 for weekly delay
cumul_death_delayed <- delay(cumul_death, 7) # 7 for weekly delay

### useful variables to output
Exposed <- E
Infectious <- I_d + I_r 

###############################
###   Initial conditions    ###
###############################

initial(S) <- N
initial(E) <- 0
initial(I_d) <- I0 / 2
initial(I_r) <- I0 / 2
initial(R) <- 0 
initial(Dead) <- 0
initial(cumul_onset) <- 0

###############################
### User defined parameters ###
###############################

N <- user(18.4e+6, min = 0)           # population size - Zambia
I0 <- user(5, min = 0)              # initial number of infected individuals 
L <- user(6.5, min = 0)             # mean incubation period
mu_d <- user(14.4, min = 0)           # time from onset to death in days
mu_r <- user(14, min = 0)           # mean time from onset to recovery in days
ifr <- user(0.0055, min = 0, max = 1) # infection fatality ratio
R0 <- user(2.7, min = 0)            # Rt (assume same for those who recover or die
# beta (Rt) values
beta_1 <- user(2.35, min = 0)
beta_2 <- user(1.2, min = 0)
beta_3 <- user(0.51, min = 0)
beta_4 <- user(3.8, min = 0)
beta_5 <- user(3.55, min = 0)
beta_6 <- user(2.5, min = 0)
beta_7 <- user(1.4, min = 0)
beta_8 <- user(0.5, min = 0)

t_beta_1 <- user(100, min = 0)  # time of beta change points
t_beta_2 <- user(130, min = 0)
t_beta_3 <- user(145, min = 0)
t_beta_4 <- user(267, min = 0)
t_beta_5 <- user(280, min = 0)
t_beta_6 <- user(300, min = 0)
t_beta_7 <- user(315, min = 0)
t_beta_8 <- user(335, min = 0)
end_time <- user()

###############################
###  Calculated parameters  ###
###############################

gamma <- 1 / L
sigma_d <- 1 / mu_d
sigma_r <- 1 / mu_r
beta_r <- Rt / mu_r
beta_d <- Rt / mu_d
Rt <- if (t >= t_beta_8)
  beta_8 else if (t >= t_beta_7)
    beta_7 else if (t >= t_beta_6)
      beta_6 else if (t >= t_beta_5)
        beta_5 else if (t >= t_beta_4)
          beta_4 else if (t >= t_beta_3)
            beta_3 else if (t >= t_beta_2)
              beta_2 else if (t >= t_beta_1)
                beta_1 else
                  R0

### additional things to output
output(Rt) <- TRUE
output(cumul_death) <- TRUE
output(weekly_onset) <- cumul_onset - cumul_onset_delayed
output(weekly_death_h) <- cumul_death - cumul_death_delayed
output(Exposed) <- TRUE
output(Infectious) <- TRUE
