# Install dependency packages if not current
rm(list=ls())		# Remove any objects in memory
james_update_packages <- function(packageList){
  availableRepos <- getCRANmirrors()
  flag <- availableRepos$Country=="USA" & grepl("https",availableRepos$URL,)
  useRepos <- sample(availableRepos$URL[flag],1)
  
  ## install missing packages, and update if newer version available
  for(i in 1:length(packageList)){
    if (!require(packageList[i],character.only = TRUE)){
      install.packages(packageList[i], repos=useRepos)
    }
  }
  
  update.packages(ask = FALSE, dependencies = c('Suggests'), oldPkgs=packageList, repos=useRepos)
}

packagelist <- c("devtools", "jsonlite", "openssl")
james_update_packages(packagelist)


# Install PSIlence from GitHub
#devtools::install_github("privacytoolsproject/PSI-Library", ref="develop") 
library("PSIlence")

num_ks = 100
epsilonGlobal <- 1
deltaGlobal <- 1e-9
best_eps = vector(mode = "numeric", length = num_ks)
for (k in 1:num_ks) {
  init <- rep(c(1/k, 0), k)
  params <- matrix(init, nrow=k, ncol=2, byrow=TRUE)
  best_eps[k] <- max(PSIlence:::update_parameters(params=params, hold=0, eps=epsilonGlobal, del=deltaGlobal)[,1]) #computing the max of the epsilons
}
noise_amounts_optimal = vector(mode = "numeric", length = num_ks)
noise_amounts_advanced = vector(mode = "numeric", length = num_ks)
noise_amounts_basic = vector(mode = "numeric", length = num_ks)
for (i in 1:num_ks) {
  noise_amounts_optimal[i] = 1/best_eps[i]
  noise_amounts_advanced[i] = sqrt(2*i*log(1/deltaGlobal))
  noise_amounts_basic[i] = i
}

plot(noise_amounts_optimal, xlab = "k", ylab = "1/epsilon (directly proportional to noise)", type="l")
lines(noise_amounts_advanced, col= "red", lty = 2)
lines(noise_amounts_basic, col = "blue", lty =3)
legend(80, 20, legend=c("optimal", "advanced", "basic"),
       col=c("black", "red", "blue"), lty=1:3, cex=0.8)
title(main = "noise comparison")
dev.copy2pdf(file="figs/psiPromises.pdf")
