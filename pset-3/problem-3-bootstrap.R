bootstrap <- function(x, y, z, n){
  index <- sample(x=1:length(x), size=n, replace=TRUE) 
  return(list(x=x[index], y=y[index], z = z[index]))
}





bootstrap_intercepts = vector(mode = "numeric", length = 100)
bootstrap_x_slopes = vector(mode = "numeric", length = 100)
bootstrap_z_slopes = vector(mode = "numeric", length = 100)

for (i in 1:100) {
  my_data = bootstrap(data.x, data.y, data.z, 10000)
  model = lm(my_data$y ~ my_data$x + my_data$z)
  bootstrap_intercepts[i] = model$coefficients[1]
  bootstrap_x_slopes[i] = model$coefficients[2]
  bootstrap_z_slopes[i] = model$coefficients[3]
}

MSE_intercept = mean(bootstrap_intercepts^2-true.intercept^2)
MSE_x_slope = mean(bootstrap_x_slopes^2-true.x.slope^2)
MSE_z_slope = mean(bootstrap_z_slopes^2-true.z.slope^2)