#'---
#' title: Multinomial Logistic Exposure Nest Success Model in Nimble
#' authors: K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'
#' **Purpose**: This script creates tracks and extracts covariates using the amt package
#' ** Last Updated**: 2/14/2026


################################################################################
## Load Packages

library(nimble)
library(coda)

################################################################################
## Simulate Data

set.seed(123)

# Number of simulated nests
n <- 300

# Continous predictor simulations
ex    <- rnorm(n, 0, 1)          # predator exposure
veg   <- rnorm(n, 0, 1)          # vegetation
init  <- rnorm(n, 0, 1)          # initiation date
interval <- sample(5:20, n, replace = TRUE)  # exposure days needed to account for exposure period

# Baseline failure rates
alpha.p <- -2.0 # Predation
alpha.a <- -2.5 # Abandonment
alpha.f <- -3.0 # Flooding

# Covariates
beta.p.ex   <- 0.7 # Exclosure effect on predation
beta.p.veg  <- -0.5 # Vegetation effect on predation
beta.a.ex   <- 0.4 # Exclosure effect on abandonment
beta.a.init <- 0.6 # Initiation Date effect on abandonment

# Calculate log hazard for each nest and failure type
log_ctp <- alpha.p + beta.p.ex * ex + beta.p.veg * veg     # Predation
log_cta <- alpha.a + beta.a.ex * ex + beta.a.init * init   # Abandonment 
log_ctf <- rep(alpha.f, n)

# Transform linear predictors to the count scale
# cts is the survival baseline
ctp <- exp(log_ctp)
cta <- exp(log_cta)
ctf <- exp(log_ctf)
cts <- rep(1, n)

# Denominator for daily probability calculation
den <- ctp + cta + ctf + cts

# Compute daily probabilities of each fate: survival, predation, abandonment, flooding
# These sum to 1 for each nest on each day
p_surv_daily  <- cts / den
p_pred_daily  <- ctp / den
p_aban_daily  <- cta / den
p_flood_daily <- ctf / den

# Convert daily survival to interval survival
# p_fail_total is the probability that a nest fails for any reason during interval
p_surv <- p_surv_daily^interval
p_fail_total <- 1 - p_surv

# Allocate competing risk probabilities conditional on each fate occurring 
fail_denom <- p_pred_daily + p_aban_daily + p_flood_daily
p_pred  <- p_fail_total * (p_pred_daily / fail_denom)
p_aban  <- p_fail_total * (p_aban_daily / fail_denom)
p_flood <- p_fail_total * (p_flood_daily / fail_denom)
probs <- cbind(p_surv, p_aban, p_pred, p_flood)

# Ensure numerical stability
probs <- probs / rowSums(probs)

# Simulate the observed fate for each nest using a multinomial distribution.
# Each nest can only have one outcome
Fate <- matrix(0, n, 4)
for(i in 1:n){
  Fate[i,] <- rmultinom(1, size = 1, prob = probs[i,])
}

# Assign column names
colnames(Fate) <- c("Survived", "Abandoned", "Predated", "Flooded")

# Create list
nestdata <- list(
  Fate = Fate,
  ex = ex,
  veg = veg,
  init = init,
  interval = interval
)

# Create dataframe
nestdata.df <- as.data.frame(nestdata)


################################################################################
## Build Logistic Exposure Nest Success Model

nestCode <- nimbleCode({
  
  # Priors for each level of failure
  alpha.p ~ dnorm(0, sd = 1)
  alpha.a ~ dnorm(0, sd = 1)
  alpha.f ~ dnorm(0, sd = 1)
  
  # Priors for each covariate
  beta.p.ex ~ dnorm(0, sd = 1)
  beta.p.veg ~ dnorm(0, sd = 1)
  beta.a.ex ~ dnorm(0, sd = 1)
  beta.a.init ~ dnorm(0, sd = 1)
  
  for(i in 1:n){
    
    # Linear predictors
    log_ctp[i] <- alpha.p + beta.p.ex * ex[i] + beta.p.veg * veg[i]
    log_cta[i] <- alpha.a + beta.a.ex * ex[i] + beta.a.init * init[i]
    log_ctf[i] <- alpha.f
    
    ctp[i] <- exp(log_ctp[i])
    cta[i] <- exp(log_cta[i])
    ctf[i] <- exp(log_ctf[i])
    cts[i] <- 1
    
    den[i] <- ctp[i] + cta[i] + ctf[i] + cts[i]
    
    # Daily probabilities for each failure type
    p_surv_daily[i] <- cts[i] / den[i]
    p_pred_daily[i] <- ctp[i] / den[i]
    p_aban_daily[i] <- cta[i] / den[i]
    p_flood_daily[i] <- ctf[i] / den[i]
    
    # Interval survival (Exposure period survival)
    p_surv[i] <- pow(p_surv_daily[i], interval[i])
    p_fail_total[i] <- 1 - p_surv[i]  # Probability that a nest fails for any reason
    
    # Allocate total failure to specific failure types
    fail_denom[i] <- p_pred_daily[i] + p_aban_daily[i] + p_flood_daily[i]
    
    p_pred[i] <- p_fail_total[i] * (p_pred_daily[i] / fail_denom[i])
    p_aban[i] <- p_fail_total[i] * (p_aban_daily[i] / fail_denom[i])
    p_flood[i] <- p_fail_total[i] * (p_flood_daily[i] / fail_denom[i])
    
    # Combine probabilities into a matrix
    probs[i,1] <- p_surv[i]
    probs[i,2] <- p_aban[i]
    probs[i,3] <- p_pred[i]
    probs[i,4] <- p_flood[i]
    
    # Likelihood: Multinomial nest fate with 4 levels
    Fate[i,1:4] ~ dmulti(probs[i,1:4], 1)
  }
})

constants <- list(n = n)

# Data
data <- list(
  Fate = Fate,
  ex = ex,
  veg = veg,
  init = init,
  interval = interval
)

# Initial values set to zero
inits <- list(
  alpha.p = 0,
  alpha.a = 0,
  alpha.f = 0,
  beta.p.ex = 0,
  beta.p.veg = 0,
  beta.a.ex = 0,
  beta.a.init = 0
)

# Create nimble object
model <- nimbleModel(
  code = nestCode,
  constants = constants,
  data = data,
  inits = inits
)

# Build and compile MCMC algorithm for model
Cmodel <- compileNimble(model)

# Parameters to estimate
conf <- configureMCMC(model, monitors = c(
  "alpha.p", "alpha.a", "alpha.f",
  "beta.p.ex", "beta.p.veg",
  "beta.a.ex", "beta.a.init"
))

# Compile complete MCMC algorithm with parameters
mcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(mcmc, project = model)

# Fit Model
samples <- runMCMC(
  Cmcmc,
  niter = 10000,
  nburnin = 3000,
  nchains = 3,
  thin = 2,
  setSeed = TRUE
)

# Convert to coda to check model diagnostics
mcmc.list.obj <- mcmc.list(
  mcmc(samples$chain1),
  mcmc(samples$chain2),
  mcmc(samples$chain3)
)

# Model checks
gelman.diag(mcmc.list.obj)
effectiveSize(mcmc.list.obj)
traceplot(mcmc.list.obj) 

################################################################################
###############################################################################X

