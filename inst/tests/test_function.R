library(Rcpp)
library(microbenchmark)

set.seed(123)

############################################################
# Dimensions
############################################################

N <- 50      # individuals
T <- 8       # time periods
L <- 3       # number of registers parameter

############################################################
# Observation matrix
############################################################

y <- matrix(sample(1:10, N*T, replace = TRUE), N, T)

############################################################
# First capture time
############################################################

first <- rep(1, N)

############################################################
# Register combinations
############################################################

combins <- matrix(sample(1:4, N*T, replace = TRUE), N, T)

############################################################
# Covariates
############################################################

covariates <- matrix(rnorm(N*11), N, 11)

############################################################
# Age and time-in covariates
############################################################

age2 <- matrix(rbinom(N*T,1,.5), N, T)
age3 <- matrix(rbinom(N*T,1,.5), N, T)
age4 <- matrix(0, N, T)

tin2 <- matrix(rbinom(N*T,1,.5), N, T)
tin3 <- matrix(rbinom(N*T,1,.5), N, T)

############################################################
# Parameters (length must match original code)
############################################################

params <- rnorm(200)

############################################################
# Transition probabilities for fast / ultrafast versions
############################################################

inv_logit <- function(x) exp(x)/(1+exp(x))

s <- matrix(inv_logit(rnorm(N*T)), N, T)
e <- matrix(inv_logit(rnorm(N*T)), N, T)
i <- matrix(inv_logit(rnorm(N*T)), N, T)
lambda <- matrix(inv_logit(rnorm(N*T)), N, T)

############################################################
# Observation misclassification probabilities
############################################################

epsilon1 <- matrix(inv_logit(rnorm(N*T)), N, T)
epsilon2 <- matrix(inv_logit(rnorm(N*T)), N, T)
epsilon3 <- matrix(inv_logit(rnorm(N*T)), N, T)

############################################################
# Mixture proportions
############################################################

mixing_prop <- c(0.6,0.4)

############################################################
# BLB counts
############################################################

counts <- rep(1,N)

############################################################
# Initial distribution
############################################################

initial <- c(1,rep(0,7))

############################################################
# Dummy design matrix
############################################################

X <- matrix(0,1,1)

############################################################
# Dummy probabilities for emissions
############################################################

proportions <- rep(1,10)

reg_block <- 2^(L - 2)
obs_block <- 64 * reg_block
probs1 <- rep(1 / obs_block, obs_block)
probs2 <- rep(1 / obs_block, obs_block)

n <- N

############################################################
# Evaluate ORIGINAL likelihood
############################################################

ll_original <- loglikelihood_cpp(
  y = y,
  first = first,
  combins = combins,
  s = s,
  e = e,
  i = i,
  lambda = lambda,
  epsilon1 = epsilon1,
  epsilon2 = epsilon2,
  epsilon3 = epsilon3,
  mixing_prop = mixing_prop,
  probs1 = probs1,
  probs2 = probs2,
  counts = counts,
  L = L
)

############################################################
# Evaluate FAST likelihood
############################################################

ll_fast <- loglikelihood_fast_cpp(
  y = y,
  first = first,
  combins = combins,
  s = s,
  e = e,
  i = i,
  lambda = lambda,
  epsilon1 = epsilon1,
  epsilon2 = epsilon2,
  epsilon3 = epsilon3,
  mixing_prop = mixing_prop,
  probs1 = probs1,
  probs2 = probs2,
  counts = counts,
  L = L
)

############################################################
# Evaluate ULTRAFAST likelihood
############################################################

ll_ultra <- loglikelihood_ultrafast(
  y = y,
  first = first,
  combins = combins,
  s = s,
  e = e,
  i = i,
  lambda = lambda,
  epsilon1 = epsilon1,
  epsilon2 = epsilon2,
  epsilon3 = epsilon3,
  mixing_prop = mixing_prop,
  probs1 = probs1,
  probs2 = probs2,
  counts = counts,
  L = L
)

############################################################
# Print results
############################################################

print(ll_original)
print(ll_fast)
print(ll_ultra)

############################################################
# Equality checks
############################################################

cat("\nOriginal vs Fast:\n")
print(all.equal(ll_original, ll_fast, tolerance = 1e-8))

cat("\nFast vs Ultrafast:\n")
print(all.equal(ll_fast, ll_ultra, tolerance = 1e-8))

############################################################
# Benchmark
############################################################

microbenchmark(

  original =
    loglikelihood_cpp(
      y,first,combins,
      s,e,i,lambda,
      epsilon1,epsilon2,epsilon3,
      mixing_prop,probs1,probs2,counts,L
    ),

  fast =
    loglikelihood_fast_cpp(
      y,first,combins,
      s,e,i,lambda,
      epsilon1,epsilon2,epsilon3,
      mixing_prop,probs1,probs2,counts,L
    ),

  ultrafast =
    loglikelihood_ultrafast(
      y,first,combins,
      s,e,i,lambda,
      epsilon1,epsilon2,epsilon3,
      mixing_prop,probs1,probs2,counts,
      L
    ),

  times = 10

)
