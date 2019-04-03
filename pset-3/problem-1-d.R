library(rmutil)
D = 100
epsilon = 0.05
x = sample.int(D, 1000, replace = TRUE)

clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
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

release_value <- function(x, D, epsilon){
  n = length(x)
  t_05 = percentileRelease(x = x, t = .05*n, lower = 1, upper = D, epsilon = epsilon)$release
  t_95 = percentileRelease(x = x, t = .95*n, lower = 1, upper = D, epsilon = epsilon)$release
  sum = 0
  for (i in 1:n) {
    if (x[i] >= t_05 && x[i] <= t_95) {
      sum = sum + (x[i])/(0.9*n)
    }
  }
  sum = sum + rlaplace(n=1, m=0, s=3*(t_95-t_05)/(0.9*n*epsilon))
  return(c(t_05, t_95, sum))
}

temp = release_value(x, D, epsilon)
t_05 = temp[1]
t_95 = temp[2]
meanRelease = temp[3]
true_mean = mean(x)
