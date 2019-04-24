library(rmutil)
epsilon = 1
#data <- read.csv(file="data/hw4testdata.csv")
data <- read.csv(file="data/CaPUMS5full.csv")
n = dim(data)[1]
d <- length(data)-1
t = -(d*log(2*(1-(0.1)^(1/d))))/(epsilon*n)

data[[2*d+1]] = data[[d+1]] 
for (i in 1:d) {
  data[[d+i]] = 1-data[[i]]
}
d = 2*d

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
y_data = data[,d+1]
x_data = data[,1:d]
probs = sapply(1:d, compute_probs, n, x_data, y_data)