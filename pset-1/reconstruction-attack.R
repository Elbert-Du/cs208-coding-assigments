pums <- read.csv(file="../../cs208/data/FultonPUMS5sample100.csv")

##
##  reconstruction attack using template by regression on sums of subsets, with randomization
##  in the form of adding gaussian noise
##
##  ED 2/23/19
##

#### Parameters ####
#set.seed(123)
n <- 100        # Dataset size
k.trials <- 2*n  # Number of queries
#q.size <- 50    # subset size, we're using predicates to query here though
sigma <- 1
num_samples <- 10
num_sigmas <- 101
prime <- 137

#### Get Data ####
## If we wanted to simulate data, we might try:
#my.pi <- 0.2    
#sensitiveData <- rbinom(n, size=1, prob=my.pi)

## But here we read it from the repository:
var <- "uscitizen"
my.pi <- mean(pums[,var])
sampleIndex <- sample(x=1:nrow(pums), size=n, replace = FALSE )
sensitiveData <- pums[sampleIndex, var]

#making pums a list for lapply
pums_list <- list()
for (i in 1:dim(pums)[1]) {
  pums_list = list.append(pums_list, pums[i,])
}

compute_predicate <- function(vec) {
  return((sum(vec * r_nums) %% prime) %% 2)
}

#### Here is our seemingly innocuous aggregated query ####
query <- function(n, data){
  r_nums <- sample(1:prime, dim(pums)[2])
  predicates <- lapply(pums_list, compute_predicate, r_nums)
  index <- which(predicates == 1)
  #index <- sample(1:length(data), n)
  subset <- data[index]
  sum <- sum(subset) + rnorm(n=1, mean=0, sd=sigma)
  return(list(sum=sum, index=index))
}



#### Here we run our query repeatedly and record results ####
true.frac = matrix(nrow = num_sigmas, ncol = num_samples)
#true0.frac = matrix(nrow = num_sigmas, ncol = num_samples)
test_sigmas = vector('numeric', num_sigmas)
for (i in 1:num_sigmas) {
  test_sigmas[i] = i-1
}
for (j in 1:num_sigmas) {
  for (k in 1:num_samples) {
    sigma = test_sigmas[j]
    
    history <- matrix(NA, nrow=k.trials, ncol=n+1)            # a matrix to store results in
    
    for(i in 1:k.trials){
      res <- query(n = q.size, data=sensitiveData)
      indicator <- 1:n %in% res$index                         # convert indices into a series of boolean/dummy variables
      indicator <- as.numeric(indicator)
      history[i,] <- c(res$sum, indicator)                    # save into our results matrix
    }
    
    #### Add (data augmentation) prior
    #s <- max(100 - k.trials, 0)
    #x <- matrix(rbinom(s*n, size=1, prob=0.5), nrow=s, ncol=n)
    #y <- apply(x, 1, sum)*my.pi
    
    
    #### Convert matrix into data frame
    xnames <- paste("x", 1:n, sep="")
    varnames<- c("y", xnames)
    releaseData <- as.data.frame(history)                     # convert matrix into data frame
    names(releaseData) <- varnames                   
    
    
    #### Run a linear regression ####
    formula <- paste(xnames, collapse=" + ")                  # construct the formula, y ~ x1 ... xn -1
    formula <- paste("y ~ ", formula, "-1")
    formula <- as.formula(formula)
    print(formula)
    
    output <- lm(formula, data=releaseData)                   # run the regression
    estimates <- output$coef                                  # save the estimates
    
    
    #### Plot results ####
    #delta <- 0.05                                             # Slight disturbance to add
    #jitterx <- runif(n=n, min=-delta, max=delta)
    #jittery <- runif(n=n, min=-delta, max=delta)
    #semi.blue <- rgb(0,90,239,200,maxColorValue=255)          # Slightly transparent colors
    #semi.red  <- rgb(239,90,0,200,maxColorValue=255)
    #col.values <- c(semi.red, semi.blue)
    
    true.1 <- (estimates>0.5) & (sensitiveData==1)            # Correctly predicted values
    true.0 <- (estimates<0.5) & (sensitiveData==0)
    true.frac[j,k] <- (round(sum(true.1, na.rm=TRUE)) + round(sum(true.0, na.rm=TRUE)))/100
    #true0.frac[j,k] <- round(sum(true.0, na.rm=TRUE)/sum(1-sensitiveData)*100)/100
    #truth.col<-1 + true.1 + true.0
  
  }
}
#filling in the sparse areas of the graph with this number of points
num_filling = 9
avg.frac = vector('numeric', num_sigmas + num_filling)
count = 0
for (j in 1:num_sigmas) {
  avg.frac[j + count] = mean(true.frac[j])
  if (j > 1 && j < num_filling + 2) {
    avg.frac[j+count+1] = (avg.frac[j+count] + mean(true.frac[j+1]))/2
    count = count+1
  }
}
inserted_values = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)
test_sigmas = sort(c(test_sigmas, inserted_values))

plot(test_sigmas, avg.frac, xlab = "sigma values", ylab = "fraction correct", ylim = c(0,1), col = "red")
#lines(test_sigmas, avg0.frac, col = "green")
#legend(1, 0.4, legend = c("ones correct", "zeros correct"), fill = c("red", "green"))
dev.copy2pdf(file="sigmaResults.pdf")

# plot(x=estimates + jitterx, y=sensitiveData + jittery, xlab="estimate", ylab="sensitive value", main="Reconstruction of uscitizen Variable", col=col.values[truth.col])    # Plot reconstruction against actual sensitive data
# abline(v=0.5, lty=2)
# text(x=0.5, y=0.8, labels=paste("fraction ones correct: ", true1.frac), pos=4)
# text(x=0.5, y=0.2, labels=paste("fraction zeros correct: ", true0.frac), pos=2)
# 
# #### Export graph to .pdf
# dev.copy2pdf(file="regAttack.pdf")





