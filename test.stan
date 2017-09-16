data {
  int N; # sample number 
  int K; # predictor number 
  vector[K] X[N]; # predictor 
  vector[K] Z; # design 
  vector[N] Y; # response
}

parameters {
  real alpha;
  real sigma;
  vector[K] beta0;
}

transformed parameters {
  vector[K] beta = Z .* beta0;
}

model {
  alpha ~ normal(0, 10);
  sigma ~ inv_gamma(.001, .001);
  beta0 ~ normal(1, 10); 
  
# print("lp before = ", target());

  for (n in 1:N) {
    target += normal_lpdf(Y[n] | alpha + sum(beta .* X[n]), sigma);
#   target += normal_lpdf(Y[n] | alpha + beta0[1] .* X[n, 1] + beta0[2] * X[n, 2], sigma);
  }

  print("lp after = ", target());
}

generated quantities {
  vector[N] lp;
  for (n in 1:N) {
#   lp[n] = normal_lpdf(Y[n] | alpha + beta0[1] .* X[n, 1] + beta0[2] * X[n, 2], sigma);
    lp[n] = normal_lpdf(Y[n] | alpha + sum(beta .* X[n]), sigma);
  }
}

