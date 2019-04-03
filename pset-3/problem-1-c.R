rm(list=ls())		# Remove any objects in memory
par(mfrow=c(1,1))   # Rebuild fresh plot window, if previously divided

#### Parameters ####

set.seed(1234)
percentile = 0.4
n.sims <- 2000

## Bound/Censor/Clip a variable to a range
clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}

## Sample with replacement from a vector
bootstrap <- function(x, y=NULL, n){
  index <- sample(x=1:length(x), size=n, replace=TRUE) 
  
  if(is.null(y)){
    return(x[index])
  }else{
    return(list(x=x[index], y=y[index]))
  }
}

## Load the data

library("foreign")
PUMSdata <- read.csv(file="data/FultonPUMS5full.csv")   
data <- PUMSdata$educ    		# variable for means
populationTrue <- quantile(data, percentile)




percentileRelease <- function(x, t, lower, upper, nbins=0, epsilon){
  n <- length(x)
  if(nbins==0){
    bins <- floor(lower):ceiling(upper)    # For integers, this is just lower:upper
    nbins <- length(bins)
  }
  x.clipped <- clip(x, lower, upper)
  sensitiveValue <- quantile(x, t/n)
  
  quality <- rep(NA, nbins)
  for(i in 1:length(quality)){
    quality[i] <- (n - abs(t-sum(x.clipped<=bins[i])))/2
  }
  likelihoods <- exp(epsilon * quality) / 2
  probabilities <- likelihoods/sum(likelihoods)
  
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities) # See also rmultinom()
  DPrelease <- min(bins[flag]) 
  
  return(list(release=DPrelease, true=sensitiveValue))
}
sample.index <- sample(1:length(data), size=100, replace=FALSE)
x <- data[sample.index]

history <- rep(NA, n.sims)
for(i in 1:n.sims){
  history[i] <- percentileRelease(x=x, t = percentile, lower=1, upper=16, epsilon=1)$release
}

## Simulate the exponential mechanism across bootstrapped datasets to see utility									

my.seq <- seq(from=log10(50), to=log10(1000), length=20)  	# make evenly spaced in logarithmic space
n.seq  <- round(10^my.seq)                                 	# round to integers

my.seq <- seq(from=log10(1), to=log10(0.01), length=5)     	# make evenly spaced in logarithmic space
ep.seq <- round(10^my.seq * 100) /100						# round to two decimal places

rawhistory <- matrix(NA, nrow=length(n.seq)*length(ep.seq)*n.sims, ncol=4)  # matrix to store results
agghistory <- matrix(NA, nrow=length(n.seq)*length(ep.seq), ncol=4)         # matrix to store results
rawcount <- 0												# counter															
aggcount <- 0                                               # counter

for(i in 1:length(n.seq)){
  for(j in 1:length(ep.seq)){
    error <- utility <- NULL
    aggcount <- aggcount + 1
    for(k in 1:n.sims){
      rawcount <- rawcount + 1
      
      ## release
      bootdata <- bootstrap(x=data, n=n.seq[i])
      DPmedian <- percentileRelease(x=bootdata, t=40, lower=1, upper=16, epsilon=ep.seq[j])
      release <- DPmedian$release
      sampleTrue <- DPmedian$true
      
      error <- c(error, sampleTrue - release)
      utility <- c(utility, min( sum(bootdata<=release), sum(bootdata>=release)) )
      
      rawhistory[rawcount, 1] <- n.seq[i]
      rawhistory[rawcount, 2] <- ep.seq[j]
      rawhistory[rawcount, 3] <- release
      rawhistory[rawcount, 4] <- sampleTrue
      
    }
    agghistory[aggcount, 1] <- n.seq[i]
    agghistory[aggcount, 2] <- ep.seq[j]
    agghistory[aggcount, 3] <- sqrt( mean( (error)^2 ) )  # RMSE
    agghistory[aggcount, 4] <- mean(utility/n.seq[i])
  }
}

plot(agghistory[,2], agghistory[,3], xlab = "epsilon", ylab = "RMSE")