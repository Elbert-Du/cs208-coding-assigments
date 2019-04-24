library(rmutil)



localRelease <- function(x, epsilon){
  draw <- runif(n=1, min=0, max=1)
  cutoff <- 1/(1+exp(epsilon))
  if(draw<cutoff){
    return(1-x)		
  }else{
    return(x)
  }
}

localBinaryData <- function(x, num_attributes, epsilon){
  n <- dim(x)[1]
  DPrelease = matrix(nrow = n, ncol = num_attributes)
  for(i in 1:n){
    for (j in 1:num_attributes) {
      temp_value <- (x[i,j] == 0 && x[i,num_attributes+1] == 1)
      DPrelease[i,j] <- localRelease(temp_value, epsilon=epsilon/(num_attributes))
    }
    
  }
  return(DPrelease)
}

compute_probs <- function(j, n, x_data) {
  my_val = 0
  for (i in 1:n) {
    if (x_data[i,j] == 1) {
      my_val = my_val + 1/n
    }
  }
  return(my_val)
}

epsilon = 1
data <- read.csv(file="data/CaPUMS5full.csv")
#data <- read.csv(file="data/hw4testdata.csv")
d <- length(data)-1
t = qnorm(0.1^(1/d), mean = (n/(1+exp(epsilon/d))), sd = sqrt(n*exp(epsilon/d)/(1+exp(epsilon/d))^2))/n

data[[2*d+1]] = data[[d+1]] 
for (i in 1:d) {
  data[[d+i]] = 1-data[[i]]
}
d = 2*d


#assumes that the last column is the y attribute
SQ_algorithm <- function(data, d, epsilon, t) {
  n = dim(data)[1]
  DP_data = localBinaryData(data, num_attributes = d, epsilon = epsilon)
  probs = sapply(1:d, compute_probs, n, DP_data)
  # for (i in 1:n) {
  #   for (j in 1:d) {
  #     if (x_data[i,j] == 0 && y_data[i] == 1) {
  #       probs[j] = probs[j] + 1/n
  #     }
  #   }
  # }
  return(list(indices = which(probs < t), full_vector = probs))
}

ns_to_test = c(10, 100, 1000, 10000, 100000, 1000000)
num_tests = 6
num_trials = 10

bootstrap_results = list()
for(i in 1:num_tests) {
  bootstrap_results[[i]] = list()
  for (j in 1:num_trials) {
    n = ns_to_test[i]
    #we have to recalculate t each time
    t = qnorm(0.9^(1/d), mean = (n/(1+exp(epsilon/d))), sd = sqrt(n*exp(epsilon/d)/(1+exp(epsilon/d))^2))/n
    bootstrap_indices = sample(1:dim(data)[1], n, replace = TRUE)
    bootstrap_data = data[bootstrap_indices,]
    bootstrap_results[[i]][[j]] = SQ_algorithm(bootstrap_data, d, epsilon, t)$indices
  }
}

target_indices =  c(6,8,10, 12, 22)

bootstrap_false_negative = vector(mode = "numeric", length = num_tests)
bootstrap_false_positive = vector(mode = "numeric", length = num_tests)
for (i in 1:num_tests) {
  for (j in 1:num_trials) {
    bootstrap_false_negative[i] = bootstrap_false_negative[i] + 5 - sum(target_indices %in% bootstrap_results[[i]][[j]])
    bootstrap_false_positive[i] = bootstrap_false_positive[i] + length(bootstrap_results[[i]][[j]]) - sum(target_indices %in% bootstrap_results[[i]][[j]])
  }
}
#a = SQ_algorithm(data, d, epsilon, t)
