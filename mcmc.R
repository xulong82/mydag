## THE OCCAM PROJECT

## Objective: build causal Bayesian networks using observational data 

library(dplyr)
library(rstan)
library(igraph)
library(binaryLogic)

# metropolis algorithm for structure inference

# 1. start from the null graph (no connection)
# 2. choosing a new graph that is proximate to the old (neighborhood proximity, randomly pick an element and switch) 
# 3. jumping to this new graph with a probability p(new)/p(old)

# simulate a bn - 4 nodes, all Gaussian process

set.seed(1)

y1 = rnorm(100, 1, 1) # input only
y2 = y1 + rnorm(100, 0, 1) # input-output
y3 = y1 + y2 + rnorm(100, 0, 1) # input-output 
y4 = y1 + y2 + y3 + rnorm(100, 0, 1) # output-only

mydag <- matrix(c(0, 0, 0, 0, # input to V1
                  1, 0, 0, 0, # input to V2
                  1, 1, 0, 0, # input to V3
                  1, 1, 1, 0  # input to V4
), nrow = 4)

mydag # true model
mydag.g = graph_from_adjacency_matrix(mydag, mode = "directed")
plot(mydag.g)

# stan model to calculate joint probability (map)

file = "~/GitHub/mydag/stan/gaussian.stan"
mymodel = stan_model(file, model_name = "occam-1")

mybic = function(dag) { # calculate bic score
  data.stan$X = dag
  fit.mle = optimizing(mymodel, init = init, data = data.stan)
  fit.mle.par = fit.mle$par
  fit.mle.par.lp = fit.mle.par[grep("^lp", names(fit.mle.par))]
  jp.map = sum(fit.mle.par.lp)
  n = nrow(dag)
  dag.k = sum(dag) # a bug?
  bic = log(n) * dag.k - 2 * jp.map
  bic
}
  
x0 <- matrix(0, ncol = 4, nrow = 4) # design matrix
y0 = as.data.frame(cbind(y1, y2, y3, y4)) # data matrix
data.stan = list(N = 100, K = 4, Y = y0, X = x0)

iter = 1e3 # iter number
init = list(alpha = rep(0, 4), sigma = rep(1, 4))

chain = list() # save dag and bic along the chain
bic = mybic(x0)
chain[[1]] = list(dag = x0, bic = bic)

for(i in 1:iter) { # if (i %% 10 == 0) cat(i, "\n")
  
  repeat{ # propose a new dag in nbh
    vec = c(dag[lower.tri(dag) | upper.tri(dag)])
    vec.idx = sample(1:length(vec), 1)
    vec[vec.idx] = 1 - vec[vec.idx]
    
    dag_new = x0 
    dag_new[lower.tri(dag_new) | upper.tri(dag_new)] = vec
    
    igraph = graph_from_adjacency_matrix(dag_new, mode = "directed")
#   print(is.dag(igraph))
    if(is.dag(igraph)) break
  } # repeat
  
  bic_new = mybic(dag_new)
  prob = exp(bic - bic_new)
    
  if (runif(1) < prob) {
    dag = dag_new
    bic = bic_new
  } else {
    dag = dag
    bic = bic
  }
  
  chain[[i+1]] = list(dag = dag, bic = bic)
}

allbic = sapply(chain, function(x) x$bic)

plot(allbic)
plot(allbic, ylim = c(1100, 1120))
plot(allbic, xlim = c(900, 1000), ylim = c(1100, 1120))

# mcmc quickly converged, found and trapped in the best dag
# the best found dag was similar and simpler than the true model
# more on bic approximation, maybe a bug
# per netfrag or per dag penalty (k)?
