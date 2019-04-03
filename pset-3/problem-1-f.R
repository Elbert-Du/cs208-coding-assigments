library(rmutil)
library(sfsmisc)

PUMSdata <- read.csv(file="data/MAPUMS5full.csv")
n = dim(PUMSdata)[1]
epsilon = 1
D = 1000000
x = matrix(nrow = n, ncol = 2)
for (i in 1:n) {
  puma = PUMSdata[i,2]
  income = PUMSdata[i,6]
  x[i,1] = puma
  x[i,2] = income
}

num_pums = length(unique(x[,1]))

my_data <- list()
pums_list = unique(x[,1])

for (i in 1:num_pums) {
  my_data[[toString(pums_list[i])]] = vector()
}

for (i in 1:n) {
  my_data[[toString(x[i,1])]] = c(my_data[[toString(x[i,1])]],x[i,2])
}

clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}

lapRelease <- function(x, lower = 0, upper, epsilon) {
  n = length(x)
  GS = upper - lower
  x.clipped = vector(mode = "numeric", length = n)
  for (i in 1:n) {
    if (x[i] < lower) {
      x.clipped[i] = lower
    } else if (x[i] > upper) {
      x.clipped[i] = upper
    } else {
      x.clipped[i] = x[i]
    }
  }
  noise = rlaplace(n=1, m=0, s=GS/(n*epsilon))
  release = mean(x.clipped) + noise
  return(release)
}

percentileRelease <- function(x, t, lower, upper, nbins=0, epsilon){
  n <- length(x)
  if(nbins==0){
    bins <- floor(lower):ceiling(upper)    # For integers, this is just lower:upper
    nbins <- length(bins)
  } else {
    bins <- seq(floor(lower),ceiling(upper), by = (upper-lower)/(nbins-1))
  }
  x.clipped <- clip(x, lower, upper)
  sensitiveValue <- quantile(x, t/n)
  
  quality <- rep(NA, nbins)
  for(i in 1:length(quality)){
    quality[i] <- (-abs(t-sum(x.clipped<=bins[i])))/2
  }
  likelihoods <- exp(epsilon * quality) / 2
  probabilities <- likelihoods/sum(likelihoods)
  
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities) # See also rmultinom()
  DPrelease <- min(bins[flag]) 
  
  return(list(release=DPrelease, true=sensitiveValue))
}

trimmedRelease <- function(x, D, epsilon){
  n = length(x)
  t_05 = percentileRelease(x = x, t = .05*n, lower = 0, upper = D, nbins = 1000, epsilon = epsilon)$release
  t_95 = percentileRelease(x = x, t = .95*n, lower = 0, upper = D, nbins = 1000, epsilon = epsilon)$release
  sum = 0
  for (i in 1:n) {
    if (x[i] >= t_05 && x[i] <= t_95) {
      sum = sum + (x[i])/(0.9*n)
    }
  }
  sum = sum + rlaplace(n=1, m=0, s=3*(t_95-t_05)/(0.9*n*epsilon))
  return(c(t_05, t_95, sum))
}

get_releases <- function(i) {
  my_release = matrix(nrow = num_pums, ncol = 2)
  
  for (i in 1:num_pums) {
    my_release[i,1] = lapRelease(x = my_data[[toString(pums_list[i])]], upper = D, epsilon = epsilon)
    my_release[i,2] = trimmedRelease(x = my_data[[toString(pums_list[i])]], D = D, epsilon = epsilon)[3]
    #print(my_release[i,])
    }
  return(my_release)
}

num_trials <- 100
releases = lapply(1:num_trials, get_releases)


puma_means = vector(mode = "numeric", length = num_pums)
for (i in 1:num_pums) {
  puma_means[i] = mean(my_data[[toString(pums_list[i])]])
}

sorted_indices = sort(puma_means, index.return = TRUE)
pums_quartiles_laplace = matrix(nrow = num_pums, ncol = 5)
for (i in 1:num_pums) {
  temp_vec = vector(mode = "numeric", length = num_trials)
  index = sorted_indices$ix[i]
  for (j in 1:num_trials) {
    temp_vec[j] = releases[[j]][index,1]
  }
  pums_quartiles_laplace[i,1] = min(temp_vec)
  pums_quartiles_laplace[i,2] = quantile(temp_vec, 0.25)
  pums_quartiles_laplace[i,3] = median(temp_vec)
  pums_quartiles_laplace[i,4] = quantile(temp_vec, 0.75)
  pums_quartiles_laplace[i,5] = max(temp_vec)
}

pums_quartiles_trim = matrix(nrow = num_pums, ncol = 5)
for (i in 1:num_pums) {
  temp_vec = vector(mode = "numeric", length = num_trials)
  index = sorted_indices$ix[i]
  for (j in 1:num_trials) {
    temp_vec[j] = releases[[j]][index,2]
  }
  pums_quartiles_trim[i,1] = min(temp_vec)
  pums_quartiles_trim[i,2] = quantile(temp_vec, 0.25)
  pums_quartiles_trim[i,3] = median(temp_vec)
  pums_quartiles_trim[i,4] = quantile(temp_vec, 0.75)
  pums_quartiles_trim[i,5] = max(temp_vec)
}
par(mfrow=c(1,1))
laplace_df = data.frame(t(pums_quartiles_laplace))
boxplot(laplace_df, xlab = "PUMA", ylab = "income", xlim = c(1,52), ylim = c(10000,70000), range = Inf)
title(main = "laplace boxplots by PUMA")
dev.copy2pdf(file="figs/laplace-boxplot.pdf")
trim_df = data.frame(t(pums_quartiles_trim))
boxplot(trim_df, xlab = "PUMA", ylab = "income", xlim = c(1,52), ylim = c(10000,70000), range = Inf)
title(main = "trimmed boxplots by PUMA")
dev.copy2pdf(file="figs/trim-boxplot.pdf")