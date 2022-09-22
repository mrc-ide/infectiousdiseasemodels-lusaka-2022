#### Imperial College London – Zambart  
#### Workshop on Analysing and modelling epidemic data 


#### Practical 2 Part 1: Estimating parameters from data.  
#### Ruth McCabe and Dr Pablo N Perez-Guzman 



## Example 1: Estimating the mean time to recovery 

## 1. Read the data into R using the following code: 

recovery_data <- c(1, 12, 6, 7, 10, 9, 21, 9, 8, 7, 3, 9, 2, 8, 25, 12, 30, 3, 3, 6, 9,
                   8, 4, 13, 7, 6, 4, 13, 37, 6, 4, 78, 6, 12, 10, 5, 21, 7, 5, 15, 7, 4,
                   23, 13, 7, 19, 8, 2, 5, 4, 1, 22, 3, 22, 3, 59, 3, 11, 20, 8, 4, 5, 16,
                   2, 23, 4, 2, 17, 3, 3, 16, 5, 2, 10, 4, 9, 2, 5, 9, 1, 2, 7, 12, 8, 8,
                   15, 8, 8, 5, 4, 7, 4, 4, 10, 16, 12, 4, 11, 11, 10)


## 2. Use a histogram to look at the distribution of the data. What continuous probability distribution best describes them? (Eg. Normal, Exponential, Gamma, …) 
hist(recovery_data, freq = FALSE,
     breaks = 20, bty = "n", xlab = "Time to recovery (days)", main="")

## 3. What is the mean recovery time? 
mean(recovery_data)



## 4. Consider how the mean recovery time changes as the sample size increases. Do you notice any convergence? 
n_sample <- 1:length(recovery_data)
estimator_over_n <- cumsum(recovery_data) / n_sample

plot(estimator_over_n,
     xlab = "Sample size",
     ylab = "Estimated mean time to recovery (days)",
     bty = "n")
abline(h = mean(recovery_data), lty = 2)

##############################################################################################################

## Example 2: Estimating the mean times to recovery or death 

## 1. Sample 1000 observations from each of these distributions. (Note: the set.seed() function ensures that we all get the same sample each time. Without this, there would be slight variability in the results.)  
set.seed(2904)

data <- data.frame(
  recovery = round(rexp(1000, 1 / 14)),
  death = round(rexp(1000, 1 / 7)))

## 2. Use the head() function in R to look at the first 20 simulated times to recovery and death. 
head(data, 20)

## 3. Using the formulae above, what are the MLEs of your two samples (time to recovery and time to death)?  
length(data$recovery)/sum(data$recovery)
length(data$death)/sum(data$death)

## 4. What does this translate to in terms of the mean number of days until the outcome? Verify your calculation by using the mean() function.  
1/(length(data$recovery)/sum(data$recovery))
1/(length(data$death)/sum(data$death))

mean(data$recovery)
mean(data$death)

## 5. Consider how the mean time to recovery and mean time to death change as the sample size increases. Do you notice any convergence? 
n_sample <- 1:nrow(data)
estimator_recovery <- cumsum(data$recovery) / n_sample
estimator_death <- cumsum(data$death) / n_sample

plot(estimator_recovery,
     xlab = "Sample size",
     ylab = "Estimated mean time to outcome (days)",
     bty = "n",
     col = "blue",
     ylim = c(0, max(estimator_recovery)))
points(estimator_death,
       col = "red")
abline(h = mean(data$recovery), col = "blue", lty = 2)
abline(h = mean(data$death), col = "red", lty = 2)

##############################################################################################################


## Example 3: Analysing truncated data 

## 1. Truncate the data generated in Example 2 to observations less than 14 days. 
truncated_recov <- data$recovery[which(data$recovery <= 14)]
truncated_death <- data$death[which(data$death <= 14)]

## 2a. How many observations are there, what is the MLE and what is the mean of: Time to recovery? 

length(truncated_recov)
length(truncated_recov)/sum(truncated_recov)
mean(truncated_recov)

## 2b. How many observations are there, what is the MLE and what is the mean of: Time to death? 

length(truncated_death)
length(truncated_death)/sum(truncated_death)
mean(truncated_death)

## 3. Consider how the mean time to recovery and mean time to death change as the sample size increases. Do you notice any convergence? 

n_recov <- 1:length(truncated_recov)
n_death <- 1:length(truncated_death)

biased_estimator_recov <- cumsum(truncated_recov) / n_recov
biased_estimator_death <- cumsum(truncated_death) / n_death

plot(biased_estimator_recov,
     xlab = "Sample size",
     ylab = "Estimated mean time to outcome (days)",
     bty = "n",
     col = "blue",
     ylim = c(0, max(estimator_recovery)))
points(biased_estimator_death,
       col = "red")
abline(h = mean(data$recovery), col = "blue", lty = 2)
abline(h = mean(data$death), col = "red", lty = 2)

##############################################################################################################


## Example 4: Adjusting for truncated data 

## 1. Run the function titled truncation_adjusted_MLE() using the code below. 

truncation_adjusted_MLE <- function(data, trunc_cut_off = 14){
  
  dL <- function(lambda, n = length(data), sum_obs = sum(data)){
    n / lambda - sum_obs - ((n*trunc_cut_off*exp(-lambda*trunc_cut_off))/(1 - exp(-lambda*trunc_cut_off)))
  }
  
  f_zero <- function(lambda){
    dL(lambda, n = length(data),sum_obs = sum(data))
  }
  
  ML_sol <- uniroot(f_zero, interval = c(1e-6, 1e6))
  
  return(list("MLE" = ML_sol$root,
              "average_outcome" = 1/ML_sol$root))
}


## 2a. Use the function to work out the MLE and average outcome time for recovery and death, and compare this to what we observed in Example 2: Time to recovery. 

truncation_adjusted_MLE(data = truncated_recov)

## 2b. Use the function to work out the MLE and average outcome time for recovery and death, and compare this to what we observed in Example 2: Time to death. 

truncation_adjusted_MLE(data = truncated_death)


##############################################################################################################


### Extension: Biased estimates under different truncation times 


truncation_adjusted_MLE_extension <- function(df = data, trunc_cut_off){
  
  # truncate original data
  data_recov_trunc <- df$recovery[which(df$recovery <= trunc_cut_off)]
  data_death_trunc <- df$death[which(df$death <= trunc_cut_off)]
  
  ## naive MLE and outcome time 
  
  ##recovery
  MLE_naive_recov <- length(data_recov_trunc)/sum(data_recov_trunc)
  average_outcome_naive_recov <- 1/MLE_naive_recov
  
  #death
  MLE_naive_death <- length(data_death_trunc)/sum(data_death_trunc)
  average_outcome_naive_death <- 1/MLE_naive_death
  
  ### adjusted for truncation
  
  ## recovery
  dL_recov <- function(lambda, n = length(data_recov_trunc), sum_obs = sum(data_recov_trunc)){
    n / lambda - sum_obs - ((n*trunc_cut_off*exp(-lambda*trunc_cut_off))/(1 - exp(-lambda*trunc_cut_off)))
  }
  
  f_zero_recov <- function(lambda){
    dL_recov(lambda, n = length(data_recov_trunc),sum_obs = sum(data_recov_trunc))
  }
  
  ML_sol_recov <- uniroot(f_zero_recov, interval = c(1e-6, 1e6))
  
  ## death
  dL_death <- function(lambda, n = length(data_death_trunc), sum_obs = sum(data_death_trunc)){
    n / lambda - sum_obs - ((n*trunc_cut_off*exp(-lambda*trunc_cut_off))/(1 - exp(-lambda*trunc_cut_off)))
  }
  
  f_zero_death <- function(lambda){
    dL_death(lambda, n = length(data_death_trunc),sum_obs = sum(data_death_trunc))
  }
  
  ML_sol_death <- uniroot(f_zero_death, interval = c(1e-6, 1e6))
  
  return(c("trunc_cut_off" = trunc_cut_off,
           "MLE_naive_recov" = MLE_naive_recov,
           "average_outcome_naive_recov" = average_outcome_naive_recov,
           "MLE_naive_death" = MLE_naive_death,
           "average_outcome_naive_death" = average_outcome_naive_death,
           "MLE_adjust_recov" = ML_sol_recov$root,
           "average_outcome_adjust_recov" = 1/ML_sol_recov$root,
           "MLE_adjust_death" = ML_sol_death$root,
           "average_outcome_adjust_death" = 1/ML_sol_death$root))
}


trunc_times <- seq(5,21,1)

sens_analysis <- c()

for(i in trunc_times){
  
  output <- truncation_adjusted_MLE_extension(df = data, trunc_cut_off = i)
  sens_analysis <- rbind(sens_analysis,output)
  
}

sens_analysis <- data.frame(sens_analysis)



ylims <- c(0, 40)
plot(sens_analysis$trunc_cut_off, sens_analysis$average_outcome_naive_recov,
     col = "blue", bty = "n", type = "l", ylim = ylims,
     ylab = "Average time to recovery (days)", xlab = "Truncation time (days)")
lines(sens_analysis$trunc_cut_off, sens_analysis$average_outcome_adjust_recov,
      col = "red")
legend("topright", legend = c("Naive", "Adjusted"), col = c("blue", "red"),
       lty = 1, bty = "n")



ylims <- c(0, max(sens_analysis$average_outcome_adjust_death))
plot(sens_analysis$trunc_cut_off, sens_analysis$average_outcome_naive_death,
     col = "blue", bty = "n", type = "l", ylim = ylims,
     ylab = "Average time to death (days)", xlab = "Truncation time (days)")
lines(sens_analysis$trunc_cut_off, sens_analysis$average_outcome_adjust_death,
      col = "red")
legend("topright", legend = c("Naive", "Adjusted"), col = c("blue", "red"),
       lty = 1, bty = "n")







