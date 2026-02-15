#'---
#' title: Multinomial Logistic exposure Nest Success Model in Nimble
#' authors: K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'
#' **Purpose**: Multinomial competing risks nest survival model with site and year random effects
#' ** Last Updated**: 2/15/2026
#' **Note**: Random effects account for variability for each competing source of nest failure

################################################################################
## Load Packages

library(nimble)
library(coda)
library(MCMCvis)

################################################################################
## Simulate Data

set.seed(123)

# Sample sizes
# 400 nests
# 4 Study Units
# 4 Years
n     <- 400
nsite <- 4
nyear <- 4

# Random assignments
site <- sample(1:nsite, n, replace = TRUE)
year <- sample(1:nyear, n, replace = TRUE)

# Continuous predictors
vis       <- rnorm(n, 0, 1) # Visual Obstruction
veg      <- rnorm(n, 0, 1) # Woody Stem Density
init     <- rnorm(n, 0, 1) # Nest Initiation Date
interval <- sample(5:27, n, replace = TRUE) 

# Baseline intercepts
alpha.pm <- -2.0 # Mammalian Nest Predation
alpha.pa <- -3.5 # Avian Nest Predation
alpha.a  <- -2.5 # Abandonment
alpha.f  <- -3.0 # Hen Killed
alpha.m  <- -1.5 # Anthropogenic Destruction

# True SDs
sigma.pm <- 0.6; tau.pm <- 0.5
sigma.pa <- 0.6; tau.pa <- 0.5
sigma.a  <- 0.6; tau.a  <- 0.5
sigma.f  <- 0.6; tau.f  <- 0.5
sigma.m  <- 0.6; tau.m  <- 0.5

# Site RE
b.pm <- rnorm(nsite, 0, sigma.pm)
b.pa <- rnorm(nsite, 0, sigma.pa)
b.a  <- rnorm(nsite, 0, sigma.a)
b.f  <- rnorm(nsite, 0, sigma.f)
b.m  <- rnorm(nsite, 0, sigma.m)

# Year RE
u.pm <- rnorm(nyear, 0, tau.pm)
u.pa <- rnorm(nyear, 0, tau.pa)
u.a  <- rnorm(nyear, 0, tau.a)
u.f  <- rnorm(nyear, 0, tau.f)
u.m  <- rnorm(nyear, 0, tau.m)

# Covariate effects
beta <- c(0.7,0.5,0.6,
          0.8,0.9,0.8,
          0.9,0.1,0.9,
          0.7,1.0,0.7,
          0.8,0.3,0.8)

# Calculate log hazard for each nest and failure type
log_ctpm <- alpha.pm + b.pm[site] + u.pm[year] + beta[1]*vis  + beta[2]*veg  + beta[3]*init
log_ctpa <- alpha.pa + b.pa[site] + u.pa[year] + beta[4]*vis  + beta[5]*veg  + beta[6]*init
log_cta  <- alpha.a  + b.a[site]  + u.a[year]  + beta[7]*vis  + beta[8]*veg  + beta[9]*init
log_ctf  <- alpha.f  + b.f[site]  + u.f[year]  + beta[10]*vis + beta[11]*veg + beta[12]*init
log_ctm  <- alpha.m  + b.m[site]  + u.m[year]  + beta[13]*vis + beta[14]*veg + beta[15]*init

# Transform linear predictors to the count scale
# cts is the survival baseline
ctpm <- exp(log_ctpm)
ctpa <- exp(log_ctpa)
cta  <- exp(log_cta)
ctf  <- exp(log_ctf)
ctm  <- exp(log_ctm)
cts  <- rep(1, n)

# Denominator for daily probability calculation
den <- ctpm + ctpa + cta + ctf + ctm + cts

# Compute daily probabilities of each fate
# These sum to 1 for each nest on each day
p_surv_daily <- cts/den
p_pm_daily   <- ctpm/den
p_pa_daily   <- ctpa/den
p_a_daily    <- cta/den
p_f_daily    <- ctf/den
p_m_daily    <- ctm/den

# Convert daily survival to interval survival
# p_fail_total is the probability that a nest fails for any reason during interval
p_surv <- p_surv_daily^interval
p_fail <- 1 - p_surv
fail_denom <- p_pm_daily + p_pa_daily + p_a_daily + p_f_daily + p_m_daily

# Allocate competing risk probabilities conditional on each fate occurring
p_pm <- p_fail*(p_pm_daily/fail_denom)
p_pa <- p_fail*(p_pa_daily/fail_denom)
p_a  <- p_fail*(p_a_daily/fail_denom)
p_f  <- p_fail*(p_f_daily/fail_denom)
p_m  <- p_fail*(p_m_daily/fail_denom)

# Ensure numerical stability
probs <- cbind(p_surv,p_pm,p_pa,p_a,p_f,p_m)
probs <- probs/rowSums(probs)

# Simulate fate
Fate <- matrix(0,n,6)
for(i in 1:n) Fate[i,] <- rmultinom(1,1,probs[i,])
colnames(Fate) <- c("Survived","Pred_Mammal","Pred_Avian",
                    "Abandoned","Hen_Killed","Human")

################################################################################
## Nimble Model

nestCode <- nimbleCode({
  
  # Priors for each level of failure
  alpha.pm ~ dnorm(0,1)
  alpha.pa ~ dnorm(0,1)
  alpha.a  ~ dnorm(0,1)
  alpha.f  ~ dnorm(0,1)
  alpha.m  ~ dnorm(0,1)
  
  # Priors for Site RE
  sigma.pm ~ dunif(0,5)
  sigma.pa ~ dunif(0,5)
  sigma.a  ~ dunif(0,5)
  sigma.f  ~ dunif(0,5)
  sigma.m  ~ dunif(0,5)
  
  # Priors for Year RE
  tau.pm ~ dunif(0,5)
  tau.pa ~ dunif(0,5)
  tau.a  ~ dunif(0,5)
  tau.f  ~ dunif(0,5)
  tau.m  ~ dunif(0,5)
  
  # Site Effect
  for(s in 1:nsite){
    b.pm[s] ~ dnorm(0, sd=sigma.pm)
    b.pa[s] ~ dnorm(0, sd=sigma.pa)
    b.a[s]  ~ dnorm(0, sd=sigma.a)
    b.f[s]  ~ dnorm(0, sd=sigma.f)
    b.m[s]  ~ dnorm(0, sd=sigma.m)
  }
  
  # Year Effect
  for(y in 1:nyear){
    u.pm[y] ~ dnorm(0, sd=tau.pm)
    u.pa[y] ~ dnorm(0, sd=tau.pa)
    u.a[y]  ~ dnorm(0, sd=tau.a)
    u.f[y]  ~ dnorm(0, sd=tau.f)
    u.m[y]  ~ dnorm(0, sd=tau.m)
  }
  
  # Betas
  for(k in 1:15){
    beta[k] ~ dnorm(0,1)
  }
  
  for(i in 1:n){
    
    log_ctpm[i] <- alpha.pm + b.pm[site[i]] + u.pm[year[i]] +
      beta[1]*vis[i] + beta[2]*veg[i] + beta[3]*init[i]
    
    log_ctpa[i] <- alpha.pa + b.pa[site[i]] + u.pa[year[i]] +
      beta[4]*vis[i] + beta[5]*veg[i] + beta[6]*init[i]
    
    log_cta[i] <- alpha.a + b.a[site[i]] + u.a[year[i]] +
      beta[7]*vis[i] + beta[8]*veg[i] + beta[9]*init[i]
    
    log_ctf[i] <- alpha.f + b.f[site[i]] + u.f[year[i]] +
      beta[10]*vis[i] + beta[11]*veg[i] + beta[12]*init[i]
    
    log_ctm[i] <- alpha.m + b.m[site[i]] + u.m[year[i]] +
      beta[13]*vis[i] + beta[14]*veg[i] + beta[15]*init[i]
    
    ctpm[i] <- exp(log_ctpm[i])
    ctpa[i] <- exp(log_ctpa[i])
    cta[i]  <- exp(log_cta[i])
    ctf[i]  <- exp(log_ctf[i])
    ctm[i]  <- exp(log_ctm[i])
    cts[i]  <- 1
    
    den[i] <- ctpm[i]+ctpa[i]+cta[i]+ctf[i]+ctm[i]+cts[i]
    
    p_surv_daily[i] <- cts[i]/den[i]
    p_pm_daily[i] <- ctpm[i]/den[i]
    p_pa_daily[i] <- ctpa[i]/den[i]
    p_a_daily[i] <- cta[i]/den[i]
    p_f_daily[i] <- ctf[i]/den[i]
    p_m_daily[i] <- ctm[i]/den[i]
    
    p_surv[i] <- pow(p_surv_daily[i], interval[i])
    p_fail[i] <- 1 - p_surv[i]
    
    fail_denom[i] <- p_pm_daily[i]+p_pa_daily[i]+p_a_daily[i]+p_f_daily[i]+p_m_daily[i]
    
    probs[i,1] <- p_surv[i]
    probs[i,2] <- p_fail[i]*(p_pm_daily[i]/fail_denom[i])
    probs[i,3] <- p_fail[i]*(p_pa_daily[i]/fail_denom[i])
    probs[i,4] <- p_fail[i]*(p_a_daily[i]/fail_denom[i])
    probs[i,5] <- p_fail[i]*(p_f_daily[i]/fail_denom[i])
    probs[i,6] <- p_fail[i]*(p_m_daily[i]/fail_denom[i])
    
    Fate[i,1:6] ~ dmulti(probs[i,1:6], 1)
  }
})

################################################################################
## Fit Model

constants <- list(n=n, nsite=nsite, nyear=nyear)
data <- list(Fate=Fate, vis=vis, veg=veg, init=init,
             interval=interval, site=site, year=year)

inits <- list(
  alpha.pm=0, alpha.pa=0, alpha.a=0, alpha.f=0, alpha.m=0,
  sigma.pm=1, sigma.pa=1, sigma.a=1, sigma.f=1, sigma.m=1,
  tau.pm=1, tau.pa=1, tau.a=1, tau.f=1, tau.m=1,
  beta=rep(0,15),
  b.pm=rnorm(nsite), b.pa=rnorm(nsite),
  b.a=rnorm(nsite), b.f=rnorm(nsite), b.m=rnorm(nsite),
  u.pm=rnorm(nyear), u.pa=rnorm(nyear),
  u.a=rnorm(nyear), u.f=rnorm(nyear), u.m=rnorm(nyear)
)

model  <- nimbleModel(nestCode, constants, data, inits)
Cmodel <- compileNimble(model)

conf  <- configureMCMC(model)
mcmc  <- buildMCMC(conf)
Cmcmc <- compileNimble(mcmc, project=model)

samples <- runMCMC(Cmcmc, niter=10000, nburnin=3000,
                   nchains=3, thin=2)

################################################################################
## Diagnostics

mcmc.list.obj <- mcmc.list(
  mcmc(samples$chain1),
  mcmc(samples$chain2),
  mcmc(samples$chain3)
)

gelman.diag(mcmc.list.obj)
effectiveSize(mcmc.list.obj)
MCMCtrace(mcmc.list.obj, pdf = FALSE)

################################################################################
###############################################################################X
