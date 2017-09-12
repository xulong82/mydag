## THE OCCAM PROJECT

## Objective: build causal Bayesian Networks using observational data 

library(abn)
library(rstan)
library(MASS)
library(binaryLogic)
library(igraph)

## Simulate a BN - 4 variables, all Gaussian process

y1 = rnorm(100, 1, 1)
y2 = 2 + 1 * y1 + rnorm(100, 0, 1) 
y3 = 3 + 2 * y1 + rnorm(100, 0, 1) 
y4 = 4 + 2 * y2 + 3 * y3 + rnorm(100, 0, 1) 

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
  dag = F
  if (sum(x) > 0 ) {
    ind = which(x == 1, arr.ind = T)
    ind.ord = apply(ind, 1, sort)
    ind.ord.pas = apply(ind.ord, 2, function(x) paste(x, collapse = "-"))
    if (! "TRUE" %in% names(table(duplicated(ind.ord.pas)))) dag = T
  }
  dag
})

table(design.dag.cons)
design.mat.dag = design.mat[design.dag.cons]

## 1 DAG tryout
## P(1) P(2|1) P(3|21) P(4|321)
## 1 ~ 0
## 2 ~ 0 + 1
## 3 ~ 0 + 1 + 2
## 4 ~ 0 + 1 + 2 + 3

dag1 <- matrix(c(0, 1, 1, 1,  # 1 to 234
                 0, 0, 1, 1,  # 2 to 34
                 0, 0, 0, 1,  # 3 to 4
                 0, 0, 0, 0   # 4 as final
              ), byrow = T, nrow = 4)

Y = rbind(y1, y2, y3, y4)
y.stan = list(N = 100, K = 4, Y = Y, X = dag1)

file = "~/bitbucket/mydag/test.stan"
mymodel = stan_model(file, model_name = "occam-1")

## Maximal likelihood estimation
xx = sapply(1:1e2, function(i) {
  fit.mle = optimizing(mymodel, data = y.stan)
  fit.mle.par = fit.mle$par
  fit.mle.par.le = fit.mle.par[grep("^lp", names(fit.mle.par))]
  sum(fit.mle.par.le)
})

## MCMC
fit.mcmc = sampling(mymodel, data = y.stan)

## All DAG posterior lp
init = list(alpha = rep(0, 4), sigma = rep(1, 4))
dag.lp = sapply(1:length(design.mat.dag), function(x) { cat(x, "\n")
  y.stan$X = design.mat.dag[[x]]
  fit.mle = optimizing(mymodel, init = init, data = y.stan)
  fit.mle.par = fit.mle$par
  fit.mle.par.le = fit.mle.par[grep("^lp", names(fit.mle.par))]
  sum(fit.mle.par.le)
})

## BIC approximation

n = 4
dag.k = sapply(design.mat.dag, sum)
dag.bic = log(n) * dag.k - 2 * dag.lp

hist(dag.lp)
hist(dag.bic)
summary(dag.bic)

## DAG
dag.stan = design.mat.dag[[which.min(dag.bic)]]

dag.stan.ind = as.data.frame(which(dag.stan == 1, arr.ind = T))
dag.igraph.dt = graph_from_data_frame(dag.stan.ind[c("col", "row")], directed = TRUE)
plot(dag.igraph.dt)

## Comparison w the ABN package

y.abn <- as.data.frame(cbind(y1, y2, y3, y4))
dist.abn <- as.list(rep("gaussian", 4))
max.par <- as.list(rep(4, 4))
names(max.par) = names(dist.abn) = names(y.abn)
ban <- matrix(rep(0, 4^2), byrow = T, nrow = 4)
rownames(ban) = colnames(ban) = names(y.abn)
## build cache 
mycache<-buildscorecache(data.df=y.abn, data.dists=dist.abn, dag.banned=ban, dag.retained=ban, max.parents=max.par); 
## find the globally best DAG 
mp.dag<-mostprobable(score.cache=mycache); 
## plot the best model - requires Rgraphviz 
myres<-fitabn(dag.m=mp.dag,data.df=y.abn,data.dists=dist.abn,create.graph=TRUE); 
plot(myres$graph)
