library(rmutil)

##
##  problem-2.r
##
##  computing the optimal upper bound for clipping
##
##  ED 2019.3.12
##


num_pulls = 200
#GS = 20
eps = 0.5

private_mean <- function(data, GS, epsilon) {
  true_mean = mean(data)
  n = length(data)
  noise <- rlaplace(n=1, m=0, s=GS/(n*epsilon))
  release = min(GS, max(0, true_mean+noise)) #clipping the outputted value as in i.
  return(release)
}

clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}


#time_between <- rexp(n = num_pulls, rate = 10)
x = rpois(n = num_pulls, lambda = 10)
# for (i in 1:num_pulls) {
#   if (i == 1) {
#     x[i] = time_between[i]
#   } else {
#     x[i] = time_between[i] + x[i-1]
#   }
# }
#my_release = private_mean(x, GS, eps)
errors = vector()
for (i in 1:100) {
  #do 10 trials for each upper bound
  MSE = 0
  GS = i
  temp_x = clip(x, 0, GS)
  for (j in 1:10) {
    my_release = private_mean(temp_x, GS, eps)
    MSE = MSE + (my_release - mean(x))^2/10
  }
  errors = c(errors, sqrt(MSE))
}
opt = which.min(errors)
plot(errors)
dev.copy2pdf(file="error-plot.pdf")