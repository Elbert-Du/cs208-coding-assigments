library(rlist)
##
##  membershipAttack.r
##
##  demonstrate membership attacks of Homer et al. (2008) and of Dwork, Smith, Steinke, Ullman, Vadhan (2015).
##
##  jH 2019.2.11
##


rm(list=ls())		# Remove any objects in memory
par(mfrow=c(1,1))   # Rebuild fresh plot window, if previously divided


pums <- read.csv(file="../../cs208/data/FultonPUMS5full.csv")

#### Parameters ####
set.seed(123)
num_people <- 100
k.attributes <- num_people*0.1
n.sims<- 1000
n.samples <- 50
null.sims <- 1000
prime <- 137
round_to <-  33/num_people #33 gave us the lowest correctness

sampleIndex <- sample(x=1:nrow(pums), size=num_people, replace = FALSE )
pums = pums[sampleIndex,]
## Generate underlying population attributes
#population.prob <- runif(n=k.attributes, min=0, max=1)
#population.mean <- 2*(population.prob-0.5)              # Because we are recoding to {-1,1} in next function

#instead we read in a sample of PUMS data for our null
#var <- "uscitizen"
#my.pi <- mean(pums[,var])
#sampleIndex <- sample(x=1:nrow(pums), size=k.attributes, replace = FALSE )



#making pums a list for lapply
pums_list <- list()
for (i in 1:dim(pums)[1]) {
  pums_list = list.append(pums_list, pums[i,])
}
#function to compute predicate in lapply
compute_predicate_value <- function(vec, r_nums=NA) {
  return((sum(vec * r_nums) %% prime) %% 2)
}
#compute a vector of predicates for each member of pums
compute_predicate <- function() {
  r_nums <- sample(1:prime, dim(pums)[2])
  predicates <- lapply(pums_list, FUN = compute_predicate_value, r_nums)
  return(predicates)
}
population.prob <- vector(mode = "numeric", length = k.attributes)
for (i in 1:k.attributes) {
  population.prob[i] <- mean(as.single(compute_predicate()), mode="integer")
}
population.mean <- 2*(population.prob-0.5)


## A utility function to create data from the population
rmvbernoulli <- function(n=1, prob){
  history <- matrix(NA, nrow=n, ncol=length(prob))
  for(i in 1:n){
    x<- rbinom(n=length(prob), size=1, prob=prob)
    x[x==0] <- -1      								# Transform from {0,1} to {-1,1}
    history[i,] <- x
  }
  return(history)
}

## Some potential test statistics
test.Homer <- function(alice, sample.mean, population.mean, referent){
  test.statistic <- sum(abs(alice - population.mean) - abs(alice - sample.mean))
  return(test.statistic)
}

# This is the Dwork et al. test using the population means
test.Dwork <- function(alice, sample.mean, population.mean, referent  = NA){
  test.statistic <- sum(alice * sample.mean) - sum(population.mean * sample.mean)
  return(test.statistic)
}

### This is the Dwork et al. test using only a single individual from the population as a referent
# test.Dwork <- function(alice, sample.mean, population.mean, referent){
# 	test.statistic <- sum(alice * sample.mean) - sum(referent * sample.mean)
# 	return(test.statistic)
# } 


## A null distribution and critical value generator
nullDistribution <- function(null.sims=1000, alpha=0.05, fun, population.prob){
  population.mean <- 2*(population.prob-0.5)
  hold <- rep(NA,null.sims)
  for(i in 1:null.sims){
    sample <- rmvbernoulli(n=n.samples, prob=population.prob)
    nullAlice <- rmvbernoulli(n=1, prob=population.prob)
    referent <- rmvbernoulli(n=1, prob=population.prob)
    sample.mean <- apply(sample, MARGIN=2, FUN=mean)
    hold[i] <- eval(fun(alice=nullAlice, sample.mean=sample.mean, population.mean=population.mean, referent=referent))
  }
  nullDistribution <- sort(hold, decreasing=TRUE)
  criticalValue <- nullDistribution[round(alpha*null.sims)]
  return(list(nullDist=nullDistribution, criticalVal=criticalValue))
}

add_noise <- function(mean) {
  new_mean <- round(mean/round_to)*round_to
  return(new_mean)
}

## Visualize the null distribution

showdist <- function(x,criticalvalue, main="", bw="nrd0"){
  testdens <- density(x, bw=bw)
  plot(testdens, main=main, xlab="Test Statistic")
  semi.blue <- rgb(0,90,239,180,maxColorValue=255)          # Slightly transparent colors
  semi.red  <- rgb(239,90,0,180,maxColorValue=255)
  flag <- testdens$x < criticalvalue
  polygon( c(min(testdens$x), testdens$x[flag], criticalvalue), y=c(0, testdens$y[flag], 0), col=semi.blue)
  polygon( c(criticalvalue, testdens$x[!flag], max(testdens$x)), y=c(0, testdens$y[!flag], 0), col=semi.red)
  abline(v=criticalvalue, lwd=1.5)
  accept.frac <- round(100*mean(x>criticalvalue))/100
  text(x=criticalvalue, y=0.8*max(testdens$y), labels=paste("OUT "), pos=2)
  text(x=criticalvalue, y=0.8*max(testdens$y), labels=paste(" IN"), pos=4)
  text(x=criticalvalue, y=0.7*max(testdens$y), labels=paste(accept.frac), pos=4)
}
# 
# 
# ## Find the null distribution for test1
# output <- nullDistribution(fun=test.Homer, population.prob=population.prob)
# testdist <- output$nullDist
# criticalValue <- output$criticalVal
# showdist(testdist, criticalValue, main="Null Distribution with Critical Value")
# 
# #### Export graph to .pdf ####
# dev.copy2pdf(file="nullDistribution.pdf")




## Simulate

history <- matrix(NA, nrow=n.sims, ncol=9)														# Matrix to store results

myalpha <- sqrt(log(num_people)/(num_people*log(10)))
#this is what we get from n=num_people, d = n, and delta = 1/(10n)


#nullDist.Homer<-nullDistribution(fun=test.Homer, population.prob=population.prob, alpha=myalpha)	# Find null distributions
nullDist.Dwork<-nullDistribution(fun=test.Dwork, population.prob=population.prob, alpha=myalpha)


for(i in 1:n.sims){
  # Simulate data
  sample <- rmvbernoulli(n=n.samples, prob=population.prob)
  sample.mean <- apply(sample, MARGIN=2, FUN=mean)
  alice <- sample[1,]
  nullAlice <- rmvbernoulli(n=1, prob=population.prob)
  #referent <- rmvbernoulli(n=1, prob=population.prob)
  
  noisy_sample.mean = add_noise(sample.mean)

  
  
  # Conduct tests
  #test.alice.Homer <- test.Homer(alice=alice, sample.mean=sample.mean, population.mean=population.mean, referent=referent)
  test.alice.Dwork <- test.Dwork(alice=alice, sample.mean=noisy_sample.mean, population.mean=population.mean)
  #test.nullAlice.Homer <- test.Homer(alice=nullAlice, sample.mean=sample.mean, population.mean=population.mean, referent=referent)
  test.nullAlice.Dwork <- test.Dwork(alice=nullAlice, sample.mean=noisy_sample.mean, population.mean=population.mean)
  
  # Store simulated values
  history[i,1]<-i
  #history[i,2]<-test.alice.Homer
  #history[i,3]<-test.alice.Homer>nullDist.Homer$criticalVal
  #history[i,4]<-test.nullAlice.Homer
  #history[i,5]<-test.nullAlice.Homer>nullDist.Homer$criticalVal
  history[i,6]<-test.alice.Dwork
  history[i,7]<-test.alice.Dwork>nullDist.Dwork$criticalVal
  history[i,8]<-test.nullAlice.Dwork
  history[i,9]<-test.nullAlice.Dwork>nullDist.Dwork$criticalVal
  
}


#### Export graph to .pdf ####

par(mfrow=c(1,2))
#showdist(history[,2], criticalvalue=nullDist.Homer$criticalVal, main="Homer Alice", bw=1)
#showdist(history[,4], criticalvalue=nullDist.Homer$criticalVal, main="Homer Null", bw=1)
Tvalue <- sqrt(8*k.attributes * log(1/myalpha))
showdist(history[,6], criticalvalue=nullDist.Dwork$criticalVal, main="Dwork Alice", bw=1)
abline(v=Tvalue, col="red", lty=2)
showdist(history[,8], criticalvalue=nullDist.Dwork$criticalVal, main="Dwork Null", bw=1)
abline(v=Tvalue, col="red", lty=2)

dev.copy2pdf(file="membershipAttack-rounding-100.pdf")





