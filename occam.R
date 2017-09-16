## THE OCCAM PROJECT
## Objective: build causal Bayesian networks using observational data 

library(abn)
library(rstan)
library(MASS)
library(binaryLogic)
library(igraph)

## Toy: BN of 4 variables
## How many DAG in the space w/o interaction?
## Enumerate the DAG by factorizing JPD of model variables: (2^3)4 w/o DAG constrain

## P(1|234) P(2|134) P(3|124) P(4|123)

design.vec = as.binary(x = (1:2^12 - 1), n = 12)

design.mat = lapply(design.vec, function(x) {
  y <- matrix(NA, ncol = 4, nrow = 4)
  diag(y) = 0
  y[lower.tri(y) | upper.tri(y)] = x
  y
})

## DAG constrains

design.dag.cons = sapply(design.mat, function(x) {
  io.only = which(rowSums(x) == 0 | colSums(x) == 0)
  igraph = graph_from_adjacency_matrix(x, mode = "directed")
  is.dag(igraph)
})

table(design.dag.cons)
design.mat.dag = design.mat[design.dag.cons]

## 1 DAG tryout
## P(1) P(2|1) P(3|21) P(4|321)
## 1 ~ 0
## 2 ~ 0 + 1
## 3 ~ 0 + 1 + 2
## 4 ~ 0 + 1 + 2 + 3

# V1: input only
# V4: output only

mydag <- matrix(c(0, 0, 0, 0, # input to V1
                  1, 0, 0, 0, # input to V2
                  1, 1, 0, 0, # input to V3
                  1, 1, 1, 0  # input to V4
               ), nrow = 4)

mydag

mydag.g = graph_from_adjacency_matrix(mydag, mode = "directed")
plot(mydag.g)
  
## Simulate a BN - 4 variables, all Gaussian process

set.seed(1)
y1 = rnorm(100, 1, 1)
y2 = y1 + rnorm(100, 0, 1) 
y3 = y1 + y2 + rnorm(100, 0, 1) 
y4 = y1 + y2 + y3 + rnorm(100, 0, 1) 

Y = as.data.frame(cbind(y1, y2, y3, y4))
y.stan = list(N = 100, K = 4, Y = Y, X = mydag)

file = "~/bitbucket/mydag/gaussian.stan"
mymodel = stan_model(file, model_name = "occam-1")

## MLE/MAP

fit.mle = optimizing(mymodel, data = y.stan)

(est.alpha = fit.mle$par[grep("^alpha\\[", names(fit.mle$par))])
(est.sigma = fit.mle$par[grep("^sigma\\[", names(fit.mle$par))])

est.beta = fit.mle$par[grep("^beta\\[", names(fit.mle$par))]
matrix(est.beta, nrow = 4)

## Check likelihood probability 

(est.lp = fit.mle$par[grep("^lp\\[", names(fit.mle$par))])
sum(est.lp)
fit.mle$value

dnorm(y.stan$Y[1, 1], mean = est.alpha[1], sd = est.sigma[1], log = T)
fit.mle$par["lp[1,1]"]

## check the MLE performance

fit.mle.dist <- sapply(1:1e2, function(i) {
  fit.mle = optimizing(mymodel, data = y.stan)
  fit.mle.par = fit.mle$par
  fit.mle.par.le = fit.mle.par[grep("^lp", names(fit.mle.par))]
  sum(fit.mle.par.le)
})

summary(fit.mle.dist)

## MCMC tryout

fit.mcmc = sampling(mymodel, data = y.stan)

## CI
stan_plot(fit.mcmc, pars = "alpha")
stan_plot(fit.mcmc, pars = "sigma")
stan_plot(fit.mcmc, pars = "beta")

## convergence
stan_trace(fit.mcmc, pars = "alpha")
stan_trace(fit.mcmc, pars = "sigma")
stan_trace(fit.mcmc, pars = "beta")

## All DAG posterior lp

init = list(alpha = rep(0, 4), sigma = rep(1, 4))
dag.lp = sapply(1:length(design.mat.dag), function(x) { cat(x, "\n")
  y.stan$X = design.mat.dag[[x]]
  fit.mle = optimizing(mymodel, init = init, data = y.stan)
  fit.mle.par = fit.mle$par
  fit.mle.par.lp = fit.mle.par[grep("^lp", names(fit.mle.par))]
  sum(fit.mle.par.lp)
})

## BIC approximation

n = 4
dag.k = sapply(design.mat.dag, sum)
dag.bic = log(n) * dag.k - 2 * dag.lp

hist(dag.lp)
hist(dag.bic)
summary(dag.bic)

## top DAG
dag.stan = design.mat.dag[[which.min(dag.bic)]]
dag.stan.g = graph_from_adjacency_matrix(dag.stan, mode = "directed")
plot(dag.stan.g) # looks good

## Comparison w the ABN package

y.abn <- Y
dist.abn <- as.list(rep("gaussian", 4))
max.par <- as.list(rep(4, 4))
names(max.par) = names(dist.abn) = names(y.abn)
ban <- matrix(rep(0, 4^2), nrow = 4)
rownames(ban) = colnames(ban) = names(y.abn)
mycache<-buildscorecache(data.df=y.abn, data.dists=dist.abn, dag.banned=ban, dag.retained=ban, max.parents=max.par)
mp.dag<-mostprobable(score.cache=mycache) 
myres<-fitabn(dag.m=mp.dag,data.df=y.abn,data.dists=dist.abn,create.graph=TRUE)
plot(myres$graph) # ?

## STAN results make more sense than ABN
