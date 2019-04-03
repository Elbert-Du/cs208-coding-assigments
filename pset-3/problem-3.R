library(rmutil)
library("foreign")
PUMSdata <- read.csv(file="data/MaPUMS5full.csv")   
D = 1000000
data <- PUMSdata$educ    		# variable for means
data.x <- PUMSdata$educ			# x-variable for regression
data.y <- PUMSdata$income		# y-variable for regression
data.y <- clip(data.y, lower = 1, upper = 1000000)
data.y <- log(data.y)
data.z <- PUMSdata$age


sgn <- function(x) {
  return(ifelse(x < 0, -1, 1))
}

normalize <- function(x){
  x[x<0] <- 0
  x <- x/sum(x)
  return(x)
}

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

showHist <- function(release, main="Histogram"){
  
  semi.blue <- rgb(0,90,239,150,maxColorValue=255)          # Slightly transparent colors
  semi.red  <- rgb(239,90,0,150,maxColorValue=255)
  
  DPrelease <- release$release
  codebook <- release$codebook
  true <- release$true
  
  allylim <- c(min(c(DPrelease,true), na.rm = TRUE), max(c(DPrelease, true), na.rm = TRUE))
  granularity <- (max(codebook) - min(codebook))/(length(codebook)-1)
  
  allxlim <- c(min(codebook) - 0.5*granularity, max(codebook + 0.5*granularity))
  
  # If stability threshold would be off the graph, extend range of graph
  if(!is.null(release$threshold)){
    if(release$threshold>allylim[2]){
      allylim[2]<-release$threshold
    }
  }
  
  # Build empty plot
  plot.new()
  plot.window( xlim=allxlim, ylim=allylim)
  title(main = main)
  axis( side=1 )
  axis( side=2 )
  
  tiny <- granularity*0.03 # slight spacing between bars
  overlap <- granularity*0.2 # some small overlap between sensitive and DP values
  
  for(i in 1:length(codebook)){
    rect(xleft=codebook[i]-overlap, ybottom=0, xright=codebook[i]+0.5*granularity-tiny, ytop=true[i], col=semi.red)
    rect(xleft=codebook[i]-0.5*granularity+tiny, ybottom=0, xright=codebook[i]+overlap, ytop=DPrelease[i], col=semi.blue)
  }
  
  # If present, show stability threshold
  if(!is.null(release$threshold)){
    abline(h=release$threshold, col="black", lty=2, lwd=1.5)
  }
}

## Differentially private histogram for continuous
xyzHistogramRelease <- regressionRelease <- function(z, y, x, zlower, zupper, ylower, yupper, xlower, xupper, znbins=0, xnbins=0, ynbins=0, epsilon){
  n <- length(x)
  if(xnbins==0){
    xlower <- floor(xlower)
    xupper <- ceiling(xupper)
    xbins <- xlower:(xupper+1)    
    xnbins <- length(xbins)-1
    xgranularity <- 1
    xcodebook <- xbins[1:xnbins]
  } else {
    xbins <- seq(from=xlower, to=xupper, length=xnbins+1)
    xgranularity <- (xupper-xlower)/xnbins
    xbins[xnbins+1] <-  xbins[xnbins+1] + xgranularity
    xcodebook <- xbins[1:xnbins] + 0.5*xgranularity
  }
  
  if(ynbins==0){
    ylower <- floor(ylower)
    yupper <- ceiling(yupper)
    ybins <- ylower:(yupper+1)    
    ynbins <- length(ybins)-1
    ygranularity <- 1
    ycodebook <- ybins[1:ynbins]
  } else {
    ybins <- seq(from=ylower, to=yupper, length=ynbins+1)
    ygranularity <- (yupper-ylower)/ynbins
    ybins[ynbins+1] <-  ybins[ynbins+1] + ygranularity
    ycodebook <- ybins[1:ynbins] + 0.5*ygranularity
  }
  
  if(znbins==0){
    zlower <- floor(zlower)
    zupper <- ceiling(zupper)
    zbins <- zlower:(zupper+1)    
    znbins <- length(zbins)-1
    zgranularity <- 1
    zcodebook <- zbins[1:znbins]
  } else {
    zbins <- seq(from=zlower, to=zupper, length=znbins+1)
    zgranularity <- (zupper-zlower)/znbins
    zbins[ynbins+1] <-  zbins[znbins+1] + zgranularity
    zcodebook <- zbins[1:znbins] + 0.5*zgranularity
  }
  
  x.clipped <- clip(x=x, lower=xlower, upper=xupper)
  y.clipped <- clip(x=y, lower=ylower, upper=yupper)
  z.clipped <- clip(x=z, lower=zlower, upper=zupper)
  
  sensitivity <- 2
  scale <- sensitivity / (epsilon)
  
  sensitiveValue <- DPrelease <- array(dim=c(znbins, ynbins, xnbins))
  
  for(i in 1:xnbins){
    for(j in 1:ynbins){
      for (k in 1:znbins) {
        sensitiveValue[k,j,i] <- sum(x.clipped >= xbins[i] & x.clipped < xbins[i+1] & y.clipped >= ybins[j] & y.clipped < ybins[j+1] & z.clipped >= zbins[k] & z.clipped < zbins[k+1])
        DPrelease[k,j,i] <- sensitiveValue[k,j,i] + rlaplace(n = 1, m=0, s=scale)
      }
    }
  }
  
  return(list(release=DPrelease, true=sensitiveValue, xcodebook=xcodebook, ycodebook=ycodebook))
}

#data1 <- bootstrap(data, n=200)
#out.1 <- xyzHistogramRelease(z = data.z, y = data.y, x=data.x, zlower = 1, zupper = 100, ylower = 0, yupper = 14, xlower = 1, xupper = 16, epsilon=0.5)

#private_counts <- out.1$release
#private_probabilities <- normalize(private_counts)
#true_counts <- out.1$true
#true_probabilities <- normalize(true_counts)



true.output <- lm(data.y~data.x + data.z)
nsims <- 20
myn<-1000


## Good Synthetic Data from example
goodhistory <-list()
for(i in 1:nsims){
  output <- xyzHistogramRelease(z = data.z, y = data.y, x=data.x, zlower = 1, zupper = 100, ylower = 0, yupper = 14, xlower = 1, xupper = 16, epsilon=0.5)
  
  syn.prob <- as.vector(normalize(output$release))
  syn.xyz <- rmultinom(n=myn, prob=syn.prob, size=1)
  xarray <- array(dim = c(dim(output$release)[3],dim(output$release)[2],dim(output$release)[1]))
  yarray <- array(dim = c(dim(output$release)[3],dim(output$release)[2],dim(output$release)[1]))
  zarray <- array(dim = c(dim(output$release)[3],dim(output$release)[2],dim(output$release)[1]))
  for (j in 1:dim(output$release)[3]) {
    for (k in 1:dim(output$release)[2]) {
      for (l in 1:dim(output$release)[1]) {
        xarray[j,k,l] = j
        yarray[j,k,l] = k
        zarray[j,k,l] = l
      }
    }
  }
  syn.x <- t(syn.xyz) %*% as.vector(xarray)
  syn.y <- t(syn.xyz) %*% as.vector(yarray)
  syn.z <- t(syn.xyz) %*% as.vector(zarray)
  
  lin_model <- lm(syn.y ~ syn.x + syn.z)
  goodhistory[i] <- lin_model
}

true.intercept <- true.output$coefficients[1]
true.x.slope <- true.output$coefficients[2]
true.z.slope <- true.output$coefficients[3]

intercepts = vector(mode = "numeric", length = nsims)
x.slopes = vector(mode = "numeric", length = nsims)
z.slopes = vector(mode = "numeric", length = nsims)
for (i in 1:nsims) {
  intercepts[i] = goodhistory[[i]][1]
  x.slopes[i] = goodhistory[[i]][2]
  z.slopes[i] = goodhistory[[i]][3]
}

intercept.bias = mean(intercepts) - true.intercept
x.slope.bias = mean(x.slopes) - true.x.slope
z.slope.bias = mean(z.slopes) - true.z.slope

intercept.variance = var(intercepts)
x.slope.variance = var(x.slopes)
z.slope.variance = var(z.slopes)

