## THE OCCAM PROJECT

## Objective: build causal Bayesian networks using observational data 

# structure mcmc

# likelihood function: bayesian score

file = "~/bitbucket/mydag/gaussian.stan"
mymodel = stan_model(file, model_name = "occam-1")

## mle/map

fit.mle = optimizing(mymodel, data = y.stan)

init = list(alpha = rep(0, 4), sigma = rep(1, 4))
dag.lp = sapply(1:length(design.mat.dag), function(x) { cat(x, "\n")
  y.stan$X = design.mat.dag[[x]]
  fit.mle = optimizing(mymodel, init = init, data = y.stan)
  fit.mle.par = fit.mle$par
  fit.mle.par.lp = fit.mle.par[grep("^lp", names(fit.mle.par))]
  sum(fit.mle.par.lp)
})

## bic approximation

n = 4
dag.k = sapply(design.mat.dag, sum)
dag.bic = log(n) * dag.k - 2 * dag.lp

# prior: neighborhood proximity

design.vec = as.binary(x = (1:2^12 - 1), n = 12)

design.mat = lapply(design.vec, function(x) {
  y <- matrix(NA, ncol = 4, nrow = 4)
  diag(y) = 0
  y[lower.tri(y) | upper.tri(y)] = x
  y
})

## dag constrains

design.dag.cons = sapply(design.mat, function(x) {
  io.only = which(rowSums(x) == 0 | colSums(x) == 0)
  igraph = graph_from_adjacency_matrix(x, mode = "directed")
  is.dag(igraph)
})

table(design.dag.cons)
design.mat.dag = design.mat[design.dag.cons]

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
