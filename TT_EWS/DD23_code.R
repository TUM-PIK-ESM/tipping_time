library(ggplot2)   ## For nice plotting
library(stats4)    ## For maximum likelihood facilities (mle())
# library(writexl)   ## For writing data frames into excel files

## Function returning 2 x the negative log-likelihood of a trace up to a constant
## The trace is assumed an Ornstein-Uhlenbeck with constant mean mu and rate
## parameter alpha. 
loglikOU.fun = function(pars, data, delta){
  ## pars = (alpha, mu, sigma) are the parameters to be estimated in the model:
  ## dXt = - alpha * ( Xt - mu ) dt + sigma * dWt
  ## data is a trace 
  ## delta is the time step between observations
  alpha0 = max(pars[1],0.001) #Ensuring alpha is positive
  mu0    = pars[2]            #Mean
  sigma2 = max(0,pars[3])     #Infinitisimal variance, should be positive
  n      = length(data)
  Xupp   = data[2:n]
  Xlow   = data[1:(n-1)]
  time   = delta*(1:(n-1))
  gamma2 = sigma2/(2*alpha0)  #Asymptotic variance
  rho0   = exp(-alpha0*delta) #Autocorrelation
  m.part = Xupp - Xlow*rho0 - mu0*(1-rho0) #Observation minus the mean of transition distribution
  v.part = gamma2*(1-rho0^2) #Variance of transition distribution
  loglik = - n*(log(v.part)) - sum(m.part^2/v.part)
  -loglik
}

## Function returning the estimated parameters
## Input is a data trace, the time step delta between observations
estimate.OU = function(data, delta, initial.values = NULL){
  if(is.null(initial.values)){
    n    = length(data)
    Xupp = data[2:n]
    Xlow = data[1:(n-1)]
    mu   = mean(data) ##Starting value for mu0 
    alpha = -log(sum((Xupp-mu)*(Xlow-mu))/sum((Xlow-mu)^2))/delta ##Starting value for alpha0 
    s2 = mean(diff(data)^2)/delta ##Quadratic variation, starting value for sigma^2
    par.init = c(alpha, mu, s2)   ## Collect starting values
  }
  if(!is.null(initial.values)) par.init = initial.values
  minuslogl = function(alpha,mu,sigma2){
    logLik = tryCatch(loglikOU.fun(pars = c(alpha,mu,sigma2), data = data, delta = delta))
    if(is.na(logLik))
      {
        return(10000)
      }else
      {
        return(logLik)
      }
    }
  temp = stats4::mle(minuslogl = minuslogl, 
             start = list(alpha = par.init[1], mu = par.init[2], sigma2 = par.init[3]))
  return(temp)
}

## Function returning 2 x the negative log-likelihood of a trace up to a constant from model
## dXt = - (a*(Xt - m)^2 + lambda) dt + sigma * dWt
## Assuming alpha0, mu0, sigma known
## Strang splitting estimator based on Taylor expansion around mu(lambda)
loglik.fun = function(pars, data, delta, alpha0, mu0, sigma20, pen = 0){
  ##pars are the parameters to be estimated, tau and a
  ##data is a trace under linear changing of lambda
  ##delta is the time step between observations
  ##alpha0, mu0, sigma20 are parameters already estimated in the OU process from stationary part
  tau     = pars[1]             #Ramping time
  a       = max(0.1,pars[2])    #Factor in front of (x-m)^2 in drift. Positive
  m       = mu0 - alpha0/(2*a)  #Constant mean shift
  lambda0 = -alpha0^2/(4*a)     #Stationary level of control parameter
  sigma2  = sigma20             #Infinitisimal variance
  n       = length(data)
  Xupp    = data[2:n]
  Xlow    = data[1:(n-1)]
  time    = delta*(1:(n-1))
  lam.seq    = lambda0*(1-time/tau)
  alpha.seq  = 2*sqrt(-a*lam.seq)
  gamma2.seq = sigma2/(2*alpha.seq)
  rho.seq    = exp(-alpha.seq*delta)
  mu.seq     = m + sqrt(-lam.seq/a)
  ## Calculating the Strang splitting scheme pseudo likelihood
    fh.half.tmp = a*delta*(Xlow - mu.seq)/2
    fh.half     = (mu.seq*fh.half.tmp+Xlow)/(fh.half.tmp+1)
    fh.half.inv = (mu.seq*fh.half.tmp-Xupp)/(fh.half.tmp-1)
    mu.h        = fh.half*rho.seq + mu.seq*(1-rho.seq)
    m.part      = fh.half.inv - mu.h
    var.part    = gamma2.seq*(1-rho.seq^2)
    det.Dfh.half.inv = 1/(a*delta*(Xupp-mu.seq)/2-1)^2
    loglik      = - sum(log(var.part)) - sum(m.part^2/var.part) + 
      2*sum(log(det.Dfh.half.inv)) - pen*n*(1/a - 1)*(a < 1)
  return(-loglik)
}

## Function returning the estimated parameters
## Input is a data trace, the time step delta between observations, and time passed
## at present since time t_0
estimate.tipping = function(data, delta, initial.values = c(100,1),  
                            alpha0, mu0, sigma20, pen = pen){
  par.init = initial.values
  minuslogl = function(pars){
    logLik = tryCatch(loglik.fun(pars = pars, data = data, delta = delta,
               alpha0 = alpha0, mu0 = mu0, sigma20 = sigma20, pen = pen))
    if(is.na(logLik))
      {
        return(10000)
      }else
      {
        return(logLik)
      }
    }
  temp = optim(par = c(tau = par.init[1], a = par.init[2]), 
               fn = minuslogl, method = "Nelder-Mead", hessian = FALSE)
  return(temp)
}

## Function returning a trajectory up to a crossing time of X(t)  
## over the barrier as a function of
## the initial condition X0, the size of the noise sigma, 
## the integration time step dt, and the length of the simulation N
## and with time varying lambda:
X.traj <- function(sigma = 0.1, lambda0 = -2, tau = 1000, m = 0, a = 1,  
                   T0 = 0, X0 = sqrt(2), dt = 0.1, Ymax = 1000000){
  ##T0: Time before ramping starts
  ##Ymax: Max number of simulated points, if tipping does not happen
  Y = 0  ## Counting integration steps up to tipping
  xbarrier = m - 2 ## One smaller than crossing point at start
  Xtraj = X0
  X = X0
  ## Simulation during stationary period, constant lambda
  while(X > xbarrier & Y < T0/dt){
    X = S.onestep(sigma = sigma, lambda = lambda0,  
                  m = m, a = a, X0 = X, dt = dt)
    Xtraj = c(Xtraj, X)
    Y = Y+1
  }
  ## Simulation after lambda starts increasing
  while(X > xbarrier & Y < Ymax){
    time = dt * (Y - T0/dt) ## Time after lambda starts increasing
    lambda = lambda0*(1-time/tau)
    X = S.onestep(sigma = sigma, lambda = lambda0*(1 - time/tau), 
                  m = m, a = a, X0 = X, dt = dt)
    Xtraj = c(Xtraj, X)
    Y = Y + 1
  }
  Y = Y*dt
  return(list(FPT = Y, X = Xtraj))
}

## Function returning a simulation of the simple model one time step
## using the Euler-Maruyama scheme
S.onestep <- function(sigma = 0.1, lambda = 0, m = 0, a = 1, X0 = 1, dt = 0.1){
  dWt = rnorm(1,0,1)    ## Draws a random increment
  Y   = X0 - a*(X0-m)^2*dt - lambda*dt + sqrt(dt)*sigma*dWt
  return(Y)
}


###############################
### Estimation on AMOC data ###
###############################

indices <- c('C18_2GMT', 'C18', 'dipole')
datasets <- c('HadISST', 'ERSSTv5', 'HadCRUT5')

ttimes <- matrix(data=NA,nrow=3,ncol=3)

for (i in 1:length(indices)) {
  for (j in 1:length(datasets)) {
    dataset <- datasets[[j]]
    index <- indices[[i]]
    file <- paste0('../TT_fingerprints/', dataset, '_amoc.txt')

    AMOC.data = read.table(file,header=TRUE)

    t0     = 1924
    data.0 = AMOC.data[AMOC.data$time <= t0, index]

    temp   = estimate.OU(data = data.0, delta = 1/12)
    alpha0 = unname(coef(temp)["alpha"])
    mu0    = unname(coef(temp)["mu"])
    sigma2 = unname(coef(temp)["sigma2"])

    data.2 = AMOC.data[AMOC.data$time > t0, index]

    temp = estimate.tipping(data = data.2, delta = 1/12, initial.values = c(100,1),
                            alpha0 = alpha0, mu0 = mu0, sigma20 = sigma2, 
                            pen = 0.004)
    tau     = unname(temp$par[1])
    a       = unname(temp$par[2])
    m       = mu0 - alpha0/(2 * a)
    lambda0 = -alpha0^2/(4 * a)
    tc      = tau + t0
    # print(dataset)
    # print(index)
    # print(tc)
    ttimes[i,j] <- tc
    # print(round(c(t0 = t0, alpha0 = alpha0, mu0 = mu0, sigma2 = sigma2, tau = tau, 
            # a = a, m = m, lambda0 = lambda0, tc = tc),2))
  }
}
rownames(ttimes) <- indices
colnames(ttimes) <- datasets
print(ttimes)

# ## Read data
# # AMOC.data = read.table("AMOCdata.txt", header = TRUE)
# AMOC.data = read.table('../TT_fingerprints/HadISST_amoc.txt',header=TRUE)
# ## time: calendar time in years
# ## AMOC0: SST (sea surface temperature) in subpolar gyre, subtracted the monthly mean
# ## AMOC1: AMOC0 subtracted the global mean SST
# ## AMOC2: AMOC0 subtracted two times the global mean SST (Arctic amplification)
# ## GM: global mean SST 

# ## Adding fingerprint with 3 times subtracted global warming for robustness analysis
# # AMOC.data$AMOC3 = AMOC.data$AMOC0 - 3 * AMOC.data$GM

# ## Estimating Ornstein-Uhlenbeck parameters from data up to year 1924
# ## before linear ramping of control parameter lambda starts

# ## Baseline data: Subset of data for the years 1870-1924
# t0     = 1924
# # data.0 = AMOC.data[AMOC.data$time <= t0, "AMOC2"]
# data.0 = AMOC.data[AMOC.data$time <= t0, "dipole"]

# ## Estimate parameters
# temp   = estimate.OU(data = data.0, delta = 1/12)
# alpha0 = unname(coef(temp)["alpha"])
# mu0    = unname(coef(temp)["mu"])
# sigma2 = unname(coef(temp)["sigma2"])

# ## Ramping data: Subset of data for the years 1925-2020
# ## Subtracted 2 times global mean
# # data.2 = AMOC.data[AMOC.data$time > t0, "AMOC2"]
# data.2 = AMOC.data[AMOC.data$time > t0, "dipole"]

# ## Estimate ramping parameters
# temp = estimate.tipping(data = data.2, delta = 1/12, initial.values = c(100,1),
#                         alpha0 = alpha0, mu0 = mu0, sigma20 = sigma2, 
#                         pen = 0.004)
# tau     = unname(temp$par[1])
# a       = unname(temp$par[2])
# m       = mu0 - alpha0/(2 * a)
# lambda0 = -alpha0^2/(4 * a)
# tc      = tau + t0

# ## Collect estimates
# round(c(t0 = t0, alpha0 = alpha0, mu0 = mu0, sigma2 = sigma2, tau = tau, 
#         a = a, m = m, lambda0 = lambda0, tc = tc),2)