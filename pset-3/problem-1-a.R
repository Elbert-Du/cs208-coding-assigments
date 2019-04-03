library(rmutil)
D = 100
epsilon = 0.05
x = sample.int(D, 1000, replace = TRUE)

release_value <- function(x, D, epsilon){
  n = length(x)
  t_05 = quantile(x, 0.05)
  t_95 = quantile(x, 0.95)
  sum = 0
  for (i in 1:n) {
    if (x[i] >= t_05 && x[i] <= t_95) {
      sum = sum + (x[i])/(0.9*n)
    }
  }
  sum = sum + rlaplace(n=1, m=0, s=D/(0.9*n*epsilon))
  return(sum)
}


print(release_value(x, D, epsilon))