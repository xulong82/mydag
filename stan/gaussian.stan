data {
  int<lower=2> N; // sample number
  int<lower=2> K; // variable number
  vector[K] Y[N]; // data matrix
  vector[K] X[K]; // design matrix
}

parameters {
  vector[K] alpha;
  vector[K] sigma;
  vector[K] beta0[K];
}

transformed parameters {
  vector[K] beta[K];
  for (k in 1:K) beta[k] = beta0[k] .* X[k];
}

model {
  alpha ~ normal(0, 10);
  sigma ~ inv_gamma(.001, .001);

  for (k in 1:K)
    beta0[k] ~ normal(0, 10);

  for (n in 1:N) {
    for (k in 1:K) {
      target += normal_lpdf(Y[n, k] | alpha[k] + sum(beta[k] .* Y[n]), sigma[k]);
    }
  }
}

generated quantities {
  vector[K] lp[N];

  for (n in 1:N) {
    for (k in 1:K) {
        lp[n, k] = normal_lpdf(Y[n, k] | alpha[k] + sum(beta[k] .* Y[n]), sigma[k]);
    }
  }
}

