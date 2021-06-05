functions {
  real generalized_inverse_gaussian_lpdf(real x, int p,
                                        real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
 }
}

data {
  int<lower=1> N;
  real x[N];
  vector[N] y;

  int<lower = 1> K; // number of components
  vector[N] observed;// is 1 if observed and 0 if missing
  vector<lower = 0>[K] alpha;
  // real<lower = 0, upper = 1> pinit; //prior on initial state
  real mu0;
  real lambda;
  real shape;
  real rate;// normal inverse gamma

  real p;// for gig prior if used
  real a;
  real b;
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho; // gaussian process part
  real<lower=0> alpha;
  vector[N] eta;
  ordered[K] mu; // mean for two states
  real<lower = 0> sigma[K];// std for two states
  simplex[K] theta[K]; // transition matrix
}
model {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);

    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  // alternative
  // target += generalized_inverse_gaussian_lpdf(rho, p, a, b)
  sigma ~ std_normal();
  eta ~ std_normal();

  y ~ normal(f, sigma);
}
