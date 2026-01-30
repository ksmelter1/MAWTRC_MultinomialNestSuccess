#---
# title: Multinomial Logistic Exposure Nest Success Model Darrah et al. (2018)
# authors: Darrah et al.
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
#' **Purpose**: This script uses simulated data to fit the Bayesian logistic exposure nest success model in Darrah et al.

################################################################################
## Load Packages

library(jagsUI)

# define simulation function with the following arguments:
# n.marked = # nests entering study
# ni = number of interations for MCMC chains
# nb = number burnin for MCMC chains
# nt = thinning rate for MCMC chains
# nc = number of parallel MCMC chains to run

nest.fate.fn.fixedeffs <- function(n.marked = 200, ni = 3000, nb = 1000, nt = 1, nc = 3){
  
  # set nest initiation dates
  # most nests initiated around day 32 - drawn from binomial dist
  
  init <- c(1:n.marked)
  for(i in 1:length(init)){
    init[i] <- rnbinom(1, size = 84, mu = 32)
  }
  
  n.occasions <- 120  # length of breeding season in days
  
  # define intercepts and covariate effects 
  alpha.p <- -4.2  # predation intercept
  alpha.a <- -8.3  # abandonment intercept
  alpha.f <- -6.4  # flooding intercept
  beta.p.ex <- -2.1  # exclosure effect on predation
  beta.p.veg <- -1.0 # vegetation effect on predation
  beta.a.ex <- 1.8   # exclosure effect on abandonment
  beta.a.init <- -1.2 # nest init. date on abandonment
  
  # assign 16% of nests to the mod-high veg category 
  veg.vec <- rep(rbinom(n.marked, 1, 0.16))
  
  # randomly assign 1/2 nests to be exclosed
  exclose.nest <- rep(rbinom(n.marked,1,0.5))
  
  # initialize empty exclosure matrix 
  ex <- matrix(0,ncol=n.occasions, nrow=n.marked) 
  
  # fill in exclosure matrix with variable age at exclosure
  for(i in 1:n.marked){
    ex[i,(init[i] + 4 + sample(1:4,1)):n.occasions]<-exclose.nest[i]
  }
  
  # scale nest initiation date
  init.sc <- as.numeric(scale(init))
  
  # prepare fate matrices and linear predictors
  # empty fate probability matrices:
  pred <- matrix(0, ncol = n.occasions, nrow = n.marked)
  aban <- matrix(0, ncol = n.occasions, nrow = n.marked)
  flood <- matrix(0, ncol = n.occasions, nrow = n.marked)
  surv <- matrix(0, ncol = n.occasions, nrow = n.marked)
  
  # linear predictors (Equation 2)
  # predation
  ctp <- exp(alpha.p + beta.p.ex*ex + beta.p.veg*veg.vec)
  
  # abandonment
  cta <- exp(alpha.a + beta.a.ex*ex + beta.a.init*init.sc)
  
  # flooding
  ctf <- exp(alpha.f)
  
  # survival
  cts <- 1
  
  # denominator of Equation 3
  den <- ctp + cta + cts + ctf
  
  # daily survival probability
  survp<-cts/den
  
  # fill in fate matrices with daily fate probabilities
  pred <- ((ctp/den)/(1-survp)*(1-survp))  
  aban <- ((cta/den)/(1-survp)*(1-survp))
  flood <- ((ctf/den)/(1-survp)*(1-survp))
  surv <- survp  # daily survival probability
  
  # work probabilities into 3-dimensional array to be used
  # in drawing daily fates for each nest
  PHI <- array(NA, dim = c(n.marked,n.occasions,4))
  PHI[,,1] <- surv
  PHI[,,2] <- pred
  PHI[,,3] <- flood
  PHI[,,4] <- aban
  # nests hatch if they survive to day 34
  hatch <- 34
  
  # initialize empty matrix of daily nest fates
  Fate.mat <- matrix(0, ncol=n.occasions, nrow=n.marked)
  
  # create vector to keep track of final nest fates
  final.fate <- rep(0,n.marked)
  
  # Fill in fate matrix
  for(i in 1:n.marked){
    # all nests are alive (status == 1) on nest intiation date
    Fate.mat[i,init[i]] <- 1  
    if (init[i] == n.occasions) next
    for(t in (init[i]+1):n.occasions){
      # draw fate from multinomial distribution
      status <- which(rmultinom(1,1,PHI[i,t-1,])== 1) 
      Fate.mat[i,t]<-status 
      # terminate nest fate history and move to next nest
      # if nest fails
      if (status == 2) break # predation
      if (status == 3) break # flood
      if (status == 4) break # abandonment
      if (t-(init[i]-1) == hatch) break  #allow nests to hatch
    }# t
    
    final.fate[i]<-max(Fate.mat[i,])
    
  }# i
  
  ###observation process: probability of discovering a nest, and ##checking every other day on average###
  
  disc <- 0.97  # probability of discovering an active nest
  freq <- 0.5   # check every other day on average 
  p <- disc*freq  # average daily probability of discovering a nest
  
  # initialize empty encounter history 
  encounter.history <- matrix(0, ncol = n.occasions, nrow = n.marked)
  discovery.date <- rep(NA, n.marked)  
  last.check <- rep(NA, n.marked)
  
  # discover nests
  for(i in 1:n.marked){
    for(t in (init[i]):n.occasions){
      # only active nests can be discovered
      if (Fate.mat[i,t] == 1)
      {encounter.history[i,t] <- rbinom(1,1,p)}
      # known nests are not subject to discovery probability
      if (encounter.history[i,t] == 1) break
      
    }# t 
    # record discovery dates in vector
    # undiscovered nests assigned a "discovery date" as last day of #season	to avoid warnings
    if (sum(encounter.history[i,]) > 0) discovery.date[i] <-which(encounter.history[i,] == 1)
    else discovery.date[i] <- n.occasions 
  }# i
  
  
  # fill in remaining encounter history
  for(i in 1:n.marked){
    for(t in (discovery.date[i]):n.occasions){
      if (discovery.date[i] == n.occasions) break
      if (Fate.mat[i,t] == 1) encounter.history[i,t] <- rbinom(1,1,freq)
    }# t
    # deal with nests never encountered to avoid warnings
    if (sum(encounter.history[i,]) == 0) 
    {final.fate[i] <- 0
    last.check[i] <- 1} 
    else 
      # record final nest fate 0-2 days after it happens	
    {last.check[i]<-(max(which(encounter.history[i,] > 0)) + sample(0:2,1)) }  	
    if (last.check[i] > n.occasions) 
      last.check[i] <- n.occasions
    encounter.history[i,last.check[i]] <- final.fate[i]
    
  }# i
  
  # rearrange data into vectorized format with covariates
  
  fate.vec <- c(t(encounter.history)) # vector of observed fates
  exclose.vec <- c(t(ex)) # vector of exclosure status
  veg.vec <- rep(veg.vec, each = n.occasions) # vegetation
  init.vec <- rep(init.sc, each = n.occasions) # nest init. date
  # vector of season dates
  days <- rep(1:n.occasions, n.marked)
  nestID <- rep(c(1:n.marked), each = n.occasions)
  # site.vec <- c(t(site))  # site ID vector
  
  # combine fate histories and covariates into data frame
  nestdata <- as.data.frame(cbind(fate.vec, days, nestID, exclose.vec, veg.vec, init.vec))
  
  # remove dates for which nests were not checked
  # (recoreded as fate = 0)
  nestdata=nestdata[nestdata$fate.vec!=0,]  
  
  # calculate exposure days
  expose<-rep(NA, length(nestdata$nestID))
  for(i in 2:length(nestdata$nestID)){
    if (nestdata$nestID[i] == nestdata$nestID[i-1]) {expose[i] <- (nestdata$days[i] - nestdata$days[i-1])}
    else {expose[i] <- NA}
  }
  
  # bind exposure days to data and remove first encounters 
  nestdata$expose <- expose
  nestdata <- nestdata[!(is.na(nestdata$expose)),]
  
  # extract and define nest fate matrix
  n <- length(nestdata$expose)
  Surv <- rep(0,n)
  Aban <- rep(0,n)
  Pred <- rep(0,n)
  Flood <- rep(0,n)
  
  for (i in 1:n){
    Surv[i][nestdata$fate.vec[i] == 1] <- 1
    Aban[i][nestdata$fate.vec[i] == 4] <- 1
    Pred[i][nestdata$fate.vec[i] == 2] <- 1
    Flood[i][nestdata$fate.vec[i] == 3] <- 1
  }
  
  Fate<-cbind(Surv, Aban, Pred, Flood)
  
  # bundle data for analysis in JAGS
  win.data <- list(ex = nestdata$exclose.vec, n = n, interval = nestdata$expose, Fate = Fate, veg = nestdata$veg.vec, init = nestdata$init.vec)
  
  # define function to draw initial values for MCMC chains
  inits <- function() {list(
    alpha.p = rnorm(1, 0, 1), alpha.a = rnorm(1, 0, 1), alpha.f = rnorm(1,0,1), beta.a.ex = rnorm(1, 0, 1), beta.p.ex = rnorm(1,0,1,), beta.p.veg = rnorm(1,0,1), beta.a.init = rnorm(1,0,1))}
  # list of parameters to be monitored
  parameters <- c("alpha.p", "alpha.a", "alpha.f", "beta.a.ex",  "beta.a.init", "beta.p.ex", "beta.p.veg")
  
  # run JAGS
  out <- jags(win.data, inits, parameters, "basic_model.txt", n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
  
  # extract and store estimated values
  output <- list(out = out, aban = sum(Aban), pred = sum(Pred), flood = sum(Flood))
  return(output)
  
} # function end

# simulations
# initiate vectors to store output results
n.sims <- 1000 #number of simulations to run
alpha.a.vec2 <- rep(NA, n.sims)
alpha.p.vec2 <- rep(NA, n.sims)
alpha.f.vec2 <- rep(NA, n.sims)
beta.a.ex.vec2 <- rep(NA, n.sims)
beta.a.init.vec2 <- rep(NA, n.sims)
beta.p.ex.vec2 <- rep(NA, n.sims)
beta.p.veg.vec2 <- rep(NA, n.sims)
aban.count2 <- rep(NA, n.sims)
pred.count2 <- rep(NA, n.sims)
flood.count2 <- rep(NA, n.sims)

# initiate simulation (runtime ~ 8.5 min per simulation or 141
# hrs total runtime on an Intel Xeon E5-1650 v3 processor with
# 32 GB RAM)
for (i in 1: n.sims){
  sim2 <- nest.fate.fn.fixedeffs(n.marked = 400, ni = 5000, nb = 2000)  
	alpha.a.vec2[i] <- sim2$out$mean$alpha.a
	alpha.p.vec2[i] <- sim2$out$mean$alpha.p
	alpha.f.vec2[i] <- sim2$out$mean$alpha.f
	beta.p.ex.vec2[i] <-sim2$out$mean$beta.p.ex
	beta.p.veg.vec2[i] <- sim2$out$mean$beta.p.veg
	beta.a.ex.vec2[i] <- sim2$out$mean$beta.a.ex
	beta.a.init.vec2[i] <- sim2$out$mean$beta.a.init
	aban.count2[i] <- sim2$aban
	pred.count2[i] <- sim2$pred
	flood.count2[i] <- sim2$flood
	cat("Finished",i,"of", n.sims,"runs") 
	}

################################################################################
## Fit Bayesian Logistic Exposure Model

model{    
  
  # Priors for predation sources of mortality
  alpha.p ~ dnorm(0,0.001)  #intercept
  beta.p.ex ~ dnorm(0, 0.001)  #exclosure slope
  beta.p.veg ~ dnorm(0, 0.001) #vegetation slope
  
  #priors for abandonment
  alpha.a ~ dnorm(0,0.001) #intercept
  beta.a.ex ~ dnorm(0,0.001) #exclosure slope
  beta.a.init ~ dnorm(0, 0.001) #nest initiation date slope
  
  #priors for flooding
  alpha.f ~ dnorm(0,0.001) #intercept
  
  # Likelihood
  for (i in 1:n) {
    #linear predictors (Equation 2)
    #flooding  
    ctf[i] <- exp(alpha.f)
    
    #predation
    ctp[i] <- exp(alpha.p + beta.p.ex*ex[i] + beta.p.veg*veg[i])
    
    #abandonment
    cta[i] <- exp(alpha.a + beta.a.ex*ex[i] + beta.a.init*init[i])
    
    #survival
    cts[i] <- 1
    
    #Equation 3 denominator
    den[i] <- ctf[i] + ctp[i] + cta[i] + cts[i]
    
    #Equation 3 daily survival probability 
    survp[i] <- cts[i]/den[i]
    
    #interval probabilities of nest loss (Equation 4)
    #flooding
    p[i,4] <- ((ctf[i]/(den[i]))/(1 - survp[i]))*(1 - pow(survp[i], interval[i]))
    #abandonment
    p[i,2] <- ((cta[i]/(den[i]))/(1 - survp[i]))*(1 - pow(survp[i], interval[i])) 
    #flooding
    p[i,3] <- ((ctp[i]/(den[i]))/(1 - survp[i]))*(1-pow(survp[i], interval[i]))
    
    #interval survival probability
    p[i,1] <- pow(survp[i], interval[i]) 
    
    Fate[i,1:4] ~ dmulti(p[i,] , 1 )
    
  }
  
}

################################################################################
###############################################################################X

