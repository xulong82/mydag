## THE OCCAM PROJECT
## Is occam assumption sufficient to identify the true causal structure (DAG) in a 3-variable setup?

library(abn)
library(rstan)
library(MASS)
library(binaryLogic)
library(igraph)

## True model
## Z -> X & Z -> Y

set.seed(1)
z = rnorm(100, 1, 3)
x = 2 * z + rnorm(100, 0, 3) 
y = 3 * z + rnorm(100, 0, 3) 

summary(lm(y ~ x))
summary(lm(y ~ x + z))

## compare likelihood of different DAGs

lp = NA

(lp[1] = logLik(lm(x ~ 1)) + logLik(lm(y ~ 1)) + logLik(lm(z ~ 1))) # P(Z)P(X)P(Y)
(lp[2] = logLik(lm(z ~ 1)) + logLik(lm(x ~ z)) + logLik(lm(y ~ z))) # P(Z)P(X|Z)P(Y|Z)
(lp[3] = logLik(lm(x ~ 1)) + logLik(lm(z ~ x)) + logLik(lm(y ~ x + z))) # P(X)P(Z|X)P(Y|X,Z)
(lp[4] = logLik(lm(x ~ 1)) + logLik(lm(y ~ x)) + logLik(lm(z ~ x + y))) # P(X)P(Y|X)P(Z|X,Y)
(lp[5] = logLik(lm(z ~ 1)) + logLik(lm(x ~ z)) + logLik(lm(y ~ x + z))) # P(Z)P(X|Z)P(Y|X,Z)
(lp[6] = logLik(lm(x ~ 1)) + logLik(lm(z ~ x)) + logLik(lm(y ~ x))) # P(Z)P(X|Z)P(Y|Z)

## lp[4] and lp[5] are the same: annoying

logLik(lm(x ~ 1)); logLik(lm(y ~ x)); logLik(lm(z ~ x + y)) # P(X)P(Y|X)P(Z|X,Y)
logLik(lm(z ~ 1)); logLik(lm(x ~ z)); logLik(lm(y ~ x + z)) # P(Z)P(X|Z)P(Y|X,Z)

## compare BIC scores

n = 3
dag.k = c(0, 2, 3, 3, 3, 2)
(bic = log(n) * dag.k - 2 * lp)

plot(lp, type = "b")
plot(bic, type = "b", col = "red")

## it actually works
## it is possible to identify the true causal structures as bayesian networks from observational data
