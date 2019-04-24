library(rmutil)
epsilon = 1
#data <- read.csv(file="data/hw4testdata.csv")
data <- read.csv(file="data/CaPUMS5full.csv")
n = dim(data)[1]
d <- length(data)-1

data[[2*d+1]] = data[[d+1]] 
for (i in 1:d) {
  data[[d+i]] = 1-data[[i]]
}
d = 2*d
t = -(d*log(2*(1-(0.9)^(1/d))))/(epsilon*n)

compute_probs <- function(j, n, x_data, y_data) {
  my_val = 0
  for (i in 1:n) {
    if (x_data[i,j] == 0 && y_data[i] == 1) {
      my_val = my_val + 1/n
    }
  }
  return(my_val)
}

#assumes that the last column is the y attribute
SQ_algorithm <- function(data, d, epsilon, t) {
  n = dim(data)[1]
  y_data = data[,d+1]
  x_data = data[,1:d]
  probs = sapply(1:d, compute_probs, n, x_data, y_data)
  for (i in 1:d) {
    probs[i] = probs[i] + rlaplace(n=1, m=0, s=d/(epsilon*n))
  }
  return(list(indices = which(probs < t), full_vector = probs))
}

ns_to_test = c(10, 100, 1000, 10000, 100000)
num_tests = 5
num_trials = 10

bootstrap_results = list()
for(i in 1:num_tests) {
  bootstrap_results[[i]] = list()
  for (j in 1:num_trials) {
    n = ns_to_test[i]
    t = -(d*log(2*(1-(0.9)^(1/d))))/(epsilon*n)
    bootstrap_indices = sample(1:dim(data)[1], n, replace = TRUE)
    bootstrap_data = data[bootstrap_indices,]
    bootstrap_results[[i]][[j]] = SQ_algorithm(bootstrap_data, d, epsilon, t)$indices
    #We have to adjust t to match the value of n we're testing and n only appears in the denominator
  }
}

target_indices = c(6,8,10, 12, 22)

bootstrap_false_negative = vector(mode = "numeric", length = num_tests)
bootstrap_false_positive = vector(mode = "numeric", length = num_tests)
for (i in 1:num_tests) {
  for (j in 1:num_trials) {
    bootstrap_false_negative[i] = bootstrap_false_negative[i] + 5 - sum(target_indices %in% bootstrap_results[[i]][[j]])
    bootstrap_false_positive[i] = bootstrap_false_positive[i] + length(bootstrap_results[[i]][[j]]) - sum(target_indices %in% bootstrap_results[[i]][[j]])
  }
}