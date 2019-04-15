# MxE.power (power calculations for DNA methylation association tests leveraging
# methylome-environment interaction) takes as input the following numeric values:
#
# - alpha.0: the intercept for simulated methylation model
# - alpha.e: square root of the proportion of variation in M (methylation) explained by E (exposure)
# - beta.0: the intercept for simulated trait model
# - beta.e: square root of the proportion of variation in Y (trait) explained by main effect of E
# - beta.m: square root of the proportion of variation in Y explained by main effect of M
# - beta.me: square root of the proportion of variation in Y explained by M-E interaction
# - sigma.m: standard deviation of error term for simulating M
# - sigma.y: standard deviation of error term for simulating Y
# - sigma.measError: standard deviation of measurement error in E
# - MCn: number of Monte Carlo draws; defaults to 1000000
# - n: sample size
# - prev.e: prevalence of exposure when binary; defaults to NA for continuous exposure
# - sens.e: sensitivity of exposure when binary (probability that a truly exposed subject is 
#           classified as exposed); defaults to NA for continuous exposure
# - spec.e: specificity of exposure when binary (probability that a truly unexposed subject
#           is classified as unexposed); defaults to NA for continuous exposure
#
# and returns a 1x3 matrix with elements:
#
# - power.M: power of the M test
# - power.ME: power of the ME test
# - power.MME: power of the joint M-ME test
#
# For a continuous exposure, prev.e, sens.e, and spec.e should all be set to NA, 
# and sigma.measError should be specified (even if it is set to 0). For a binary exposure, 
# prev.e, sens.e, and spec.e should be specified, and sigma.measError should be set to NA. 
# The distribution of a continuous exposure can be edited as denoted in the function.
#
# An example of implementation of this function is provided at the end. The example
# computes the power of the tests for a continuous exposure when there is no measurement
# error in the exposure and no methylation-environment dependence.


MxE.power = function(
  alpha.0,
  alpha.e,
  beta.0,
  beta.e,
  beta.m,
  beta.me,
  sigma.m,
  sigma.y,
  sigma.measError=NA, ## Will not be used if exposure is binary
  MCn=1000000,
  n,
  prev.e=NA, ## Set prev.e=NA to indicate exposure is continuous
  sens.e=NA, ## Will not be used if exposure is continuous
  spec.e=NA) ## Will not be used if exposure is continuous
{
  
  ## Generate exposure data
  
  if(is.na(prev.e)){
    # Continuous expsoure
    E <- rnorm(MCn, mean=0, sd=1) # Edit this line to change continuous distribution of exposure
    if(is.na(sigma.measError)){
      return("Error: Must specify standard deviation of measurement error of exposure.")
    }else{
      Estar <- E + rnorm(MCn, mean=0, sd=sigma.measError)
    } 
  }else{
    if(prev.e<=0 || prev.e>=1){
      return("Error: Prevalence for binary exposure must be between 0 and 1.")
    }else{
      # Binary exposure (uncentered)
      E.uncent <- rbinom(MCn, 1, prob=prev.e)
      Estar.unexp <- rbinom(sum(E.uncent==0), 1, prob=1-spec.e)
      Estar.exp <- rbinom(sum(E.uncent==1), 1, prob=sens.e)
      Estar.uncent = c(Estar.unexp,Estar.exp)
      
      # Standardize the binary exposure
      E <- sort((E.uncent - prev.e) / sqrt(prev.e*(1-prev.e)))
      Estar <- (Estar.uncent - prev.e) / sqrt(prev.e*(1-prev.e))
      if(sum(E.uncent==0)==MCn | sum(Estar.uncent==0)==MCn) return("Error: Prevlance for binary expsoure is too small.")
      if(sum(E.uncent==1)==MCn | sum(Estar.uncent==1)==MCn) return("Error: Prevalence for binary exposure is too large.")
    }
  }
  
  ## Generate methylation and trait data
  
  M <- alpha.0 + alpha.e*E + rnorm(MCn, mean=0, sd=sigma.m)
  Y <- beta.0 + beta.e*E + beta.m*M + beta.me*E*M + rnorm(MCn, mean=0, sd=sigma.y)
  
  ## Fit full model (Y ~ E + M + EM)
  
  X = cbind(1, Estar, M, Estar*M)
  betahat.full = solve(t(X)%*%X) %*% t(X) %*% Y
  muhat.full = X%*%betahat.full
  sigmasquaredhat.full = as.numeric((t(Y - muhat.full) %*% (Y - muhat.full)) / (MCn - ncol(X)))
  logL.full = sum(-log(sqrt(2*pi)) - log(sqrt(sigmasquaredhat.full)) - ((Y - muhat.full)^2/(2*sigmasquaredhat.full)))
  
  ## Fit main effects model (Y ~ E + M)
  
  X = cbind(1, Estar, M)
  betahat.maineffects = solve(t(X)%*%X) %*% t(X) %*% Y 
  muhat.maineffects = X%*%betahat.maineffects
  sigmasquaredhat.maineffects = as.numeric((t(Y - muhat.maineffects) %*% (Y - muhat.maineffects)) / (MCn - ncol(X)))
  logL.maineffects = sum(-log(sqrt(2*pi)) - log(sqrt(sigmasquaredhat.maineffects)) - ((Y - muhat.maineffects)^2/(2*sigmasquaredhat.maineffects)))
  
  ## Fit model with just exposure (Y ~ E)
  
  X = cbind(1, Estar)
  betahat.justE = solve(t(X)%*%X) %*% t(X) %*% Y 
  muhat.justE = X%*%betahat.justE
  sigmasquaredhat.justE = as.numeric((t(Y - muhat.justE) %*% (Y - muhat.justE)) / (MCn - ncol(X)))
  logL.justE = sum(-log(sqrt(2*pi)) - log(sqrt(sigmasquaredhat.justE)) - ((Y - muhat.justE)^2/(2*sigmasquaredhat.justE)))
  
  ## Power for M-ME test
  
  ncp.MME = 2*n*(logL.full - logL.justE)/MCn
  c = qchisq(0.95, df=2, ncp=0)
  power.MME = pchisq(c, df=2, ncp=max(ncp.MME,0), lower.tail=FALSE)
  
  ## Power for ME test
  
  ncp.ME = 2*n*(logL.full - logL.maineffects)/MCn
  c = qchisq(0.95, df=1, ncp=0)
  power.ME = pchisq(c, df=1, ncp=max(ncp.ME,0), lower.tail=FALSE)
  
  ## Power for M test
  
  ncp.M = 2*n*(logL.maineffects - logL.justE)/MCn
  c = qchisq(0.95, df=1, ncp=0)
  power.M = pchisq(c, df=1, ncp=max(ncp.M,0), lower.tail=FALSE)
  
  ## Return power of M, ME, and M-ME tests
  
  power = as.matrix(data.frame(power.M, power.ME, power.MME))
  return(power)
}



##
## Example for continuous exposure where there is no measurement
## error in the exposure and no methylation-environment dependence
##

my_alpha.0 <- 0
my_alpha.e <- 0
my_beta.0 <- 0
my_beta.e <- sqrt(0.025)
my_beta.m <- sqrt(0.01)
my_beta.me <- sqrt(0.01)
my_sigma.m <- sqrt(1-my_alpha.e^2) # this forces Var(M)=1
my_sigma.y <- sqrt(1-(my_beta.e^2+my_beta.m^2+(my_beta.me^2)*(1+my_alpha.e^2))) # this forces Var(Y)=1
my_sigma.measError <- 0
my_MCn <- 1000000
my_n <- 500
my_prev.e <- NA
my_sens.e <- NA
my_spec.e <- NA

MxE.power(
  alpha.0=my_alpha.0, 
  alpha.e=my_alpha.e, 
  beta.0=my_beta.0, 
  beta.e=my_beta.e, 
  beta.m=my_beta.m, 
  beta.me=my_beta.me, 
  sigma.m=my_sigma.m, 
  sigma.y=my_sigma.y, 
  sigma.measError=my_sigma.measError, 
  MCn=my_MCn, 
  n=my_n, 
  prev.e=my_prev.e, 
  sens.e=my_sens.e, 
  spec.e=my_spec.e)
