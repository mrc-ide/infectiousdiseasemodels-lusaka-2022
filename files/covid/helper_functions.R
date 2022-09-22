
# function to process the data
process_data <- function(df, geography, province = FALSE) {
  
  # Get column names and format names to manipulate easily
  columns <- tolower(colnames(df))
  columns <- gsub("x", "", columns)
  
  # Assign these column names to data frame
  colnames(df) <- columns
  
  # Get a vector of dates from corresponding columns
  dates <- columns[5:length(columns)]
  dates <- lubridate::mdy(dates)
  
  # Get data for desired geography
  if (province) {
    data <- as.numeric(df[df$province.state == geography, -c(1:4)])
  } else {
    data <- as.numeric(df[df$country.region == geography, -c(1:4)])
  }
  
  data.frame(
    date = dates,
    cumulative = data,
    incidence = c(0, diff(data)))
  
}

# function to plot time series incidence + cumulative
plot_inc_cum <- function(df, what) {
  
  # some graphical parameters to make plots look better
  par(mar = c(5, 5, 2, 5))
  
  # plot
  plot(df$date, df$incidence, main = paste(what),
       type = "b", bty = "n", xlab = "Date", ylab = "Incidence")
  par(new = TRUE)
  plot(df$date, df$cumulative, type = "l", col = "red",
       axes = FALSE, ylab = "", xlab = "")
  axis(side = 4)
  mtext(side = 4, line = 4, "Cumulative", col = "red")
  
}

# function to return negative log-likelihood of I given model parameters
ll_infected <- function(S0, beta, I) {
  
  n <- length(I)
  
  S <- floor(S0 - cumsum(I[-n]))
  
  p <- 1 - exp(-beta * I[-n] / S0)
  
  ll <- -sum(dbinom(I[-1], S, p , log = TRUE))
  
  ll
}
