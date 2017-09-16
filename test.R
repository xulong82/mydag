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
