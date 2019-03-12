library(rmutil)

##
##  problem-3.r
##
##  first part of the code tests many different epsilons, the second part runs it for the even split and produces a graph.
##
##  ED 2019.3.12
##



num_pulls = 1000
sigma = 1
alpha = 1
beta = 1
epsilon = 0.25 #We use this 3 times independently
num_samples = 100
my_epsilons = list()

compute_MSE <- function(a,b,x,y) {
  MSE = 0
  for (i in 1:num_pulls) {
    MSE = MSE + (y[i]-(a+b*x[i]))^2/num_pulls
  }
  return(MSE)
}

clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}


for (i in 20:30) {
  for (j in 20:30) {
    for (k in 20:30) {
      l = 100 - i - j - k
      if (20 <= l && l <= 30) {
        my_epsilons[[length(my_epsilons)+1]] <- list(c(i/100,j/100,k/100,l/100))
      }
    }

  }
}



mean_private_MSEs = vector()
mean_non_private_MSEs = vector()
for (epsilons in my_epsilons) {
  epsilon_1 = epsilons[[1]][1]
  epsilon_2 = epsilons[[1]][2]
  epsilon_3 = epsilons[[1]][3]
  epsilon_4 = epsilons[[1]][4]
  private_MSEs = vector()
  non_private_MSEs = vector()
  for (j in 1:num_samples) {

    x <- rpois(n = num_pulls, lambda = 10)
    
    x_small_sample <- sample(x, 10, replace = FALSE)
    GS <- qpois(1-1/num_pulls, lambda = mean(x_small_sample))
    x <- clip(x, 0, GS)
    
    Sxx = sum((x-mean(x))^2) + rlaplace(n=1, m=0, s=(1-1/num_pulls)/epsilon_1)
    
    noise <- rnorm(n = num_pulls, mean = 0, sd = sigma)
    y = beta*x+alpha+noise
    y_small_sample <- sample(y, 10, replace = FALSE)
    y_min_est <- qpois(1/num_pulls, lambda = mean(y_small_sample))
    y_max_est <- qpois(1-1/num_pulls, lambda = mean(y_small_sample))
    
    y <- clip(y, y_min_est, y_max_est)
    
    Sxy = sum((x-mean(x))*(y-mean(y))) + rlaplace(n=1, m=0, s=(1-1/num_pulls)/epsilon_2)
    b = Sxy/Sxx
    
    noisy_mean_x = min(GS, max(0, mean(x) + rlaplace(n=1, m=0, s=GS/(num_pulls*epsilon_3))))
    noisy_mean_y = min(y_max_est, max(y_min_est, mean(y) + rlaplace(n=1, m=0, s=(y_max_est-y_min_est)/(num_pulls*epsilon_4))))
    
    a = noisy_mean_y - b*noisy_mean_x
    
    private_MSEs = c(private_MSEs, compute_MSE(a,b,x,y))
    
    simple_model = lm(y ~ x)
    a1 = simple_model$coefficients[1]
    b1 = simple_model$coefficients[2]
    non_private_MSEs = c(non_private_MSEs, compute_MSE(a1,b1,x,y))
  }
  mean_private_MSEs = c(mean_private_MSEs, mean(private_MSEs))
  mean_non_private_MSEs = c(mean_non_private_MSEs, mean(non_private_MSEs))
}

o = which.min(mean_private_MSEs - mean_non_private_MSEs)
print((mean_private_MSEs - mean_non_private_MSEs)[[o]])

# private_MSEs = vector()
# non_private_MSEs = vector()
# for (j in 1:num_samples) {
#   #print(j)
#   x <- rpois(n = num_pulls, lambda = 10)
# 
#   x_small_sample <- sample(x, 10, replace = FALSE)
#   GS <- qpois(1-1/num_pulls, lambda = mean(x_small_sample))
#   x <- clip(x, 0, GS)
#   
#   Sxx = sum((x-mean(x))^2) + rlaplace(n=1, m=0, s=(1-1/num_pulls)/epsilon)
# 
#   noise <- rnorm(n = num_pulls, mean = 0, sd = sigma)
#   y = beta*x+alpha+noise
#   y_small_sample <- sample(y, 10, replace = FALSE)
#   y_min_est <- qpois(1/num_pulls, lambda = mean(y_small_sample))
#   y_max_est <- qpois(1-1/num_pulls, lambda = mean(y_small_sample))
#   
#   y <- clip(y, y_min_est, y_max_est)
#   
#   Sxy = sum((x-mean(x))*(y-mean(y))) + rlaplace(n=1, m=0, s=(1-1/num_pulls)/epsilon)
#   b = Sxy/Sxx
# 
#   noisy_mean_x = min(GS, max(0, mean(x) + rlaplace(n=1, m=0, s=GS/(num_pulls*epsilon))))
#   noisy_mean_y = min(y_max_est, max(y_min_est, mean(y) + rlaplace(n=1, m=0, s=(y_max_est-y_min_est)/(num_pulls*epsilon))))
# 
#   a = noisy_mean_y - b*noisy_mean_x
# 
#   private_MSEs = c(private_MSEs, compute_MSE(a,b,x,y))
#   
#   simple_model = lm(y ~ x)
#   a1 = simple_model$coefficients[1]
#   b1 = simple_model$coefficients[2]
#   non_private_MSEs = c(non_private_MSEs, compute_MSE(a1,b1,x,y))
# }
# 
# 
# par(mfrow=c(1,2))
# hist(private_MSEs, breaks = 10)
# abline(v=mean(private_MSEs), col="red", lty=2)
# hist(non_private_MSEs, breaks = 10)
# abline(v=mean(non_private_MSEs), col="red", lty=2)
# 
# print(mean(private_MSEs - non_private_MSEs))
# 
# dev.copy2pdf(file="linreg-plot.pdf")