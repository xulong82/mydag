## THE OCCAM PROJECT

## Objective: build causal Bayesian networks using observational data 

library(abn)
library(rstan)
library(MASS)
library(binaryLogic)
library(igraph)

set.seed(1)
y1 = rnorm(100, 1, 1)
y2 = y1 + rnorm(100, 0, 1) 
y3 = y1 + y2 + rnorm(100, 0, 1) 
y4 = y1 + y2 + y3 + rnorm(100, 0, 1) 

Y = as.data.frame(cbind(y1, y2, y3, y4))
y.stan = list(N = 100, K = 2, Y = Y$y4, X = Y[-c(1, 4)], Z = c(1, 1))

y.stan$Z = c(1, 1)
summary(lm(y4 ~ y2 + y3, data = Y))

y.stan$Z = c(0, 1)
summary(lm(y4 ~ y3, data = Y))

file = "~/bitbucket/mydag/test.stan"
mymodel = stan_model(file, model_name = "occam-1")

fit.mle = optimizing(mymodel, data = y.stan)

fit.mle$par["alpha"]
(est.beta = fit.mle$par[grep("^beta\\[", names(fit.mle$par))])

est.lp = fit.mle$par[grep("^lp\\[", names(fit.mle$par))]
sum(est.lp)
fit.mle$value

dnorm(y.stan$Y[1], mean = fit.mle$par["alpha"] + sum(est.beta * y.stan$X[1, ]), sd = fit.mle$par["sigma"], log = T)
fit.mle$par["lp[1]"]

# structure mcmc

# likelihood function: bayesian score

# prior: neighborhood proximity

# posterior 

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

# metropolis algorithm
# 1. starting at a random graph
# 2. choosing a new graph that is proximate to the old 
# 3. jumping to this new point with a probability p(new)/p(old), where p is the target function, and p>1 means jumping as well

metropolis <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1, 3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = rnorm(3, mean = chain[i,], sd = c(0.5,0.1,0.3))
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    
    if (runif(1) < probab) {
      chain[i+1,] = proposal
    } else {
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(0, 5, 10)
chain = run_metropolis_MCMC(startvalue, 1e4)

burnIn = 5e3
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
