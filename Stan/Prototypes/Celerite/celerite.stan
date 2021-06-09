functions {
   //real logLikSHO(vector t, vector y,real S0, real w0, real Q, real eps, vector diag);
   vector dotCholRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
   real logLikRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
}

data {
  int<lower=1> N;
  vector[N] t; // time
  vector[N] y; // light curve
  //log uniform priors
  vector[2] sigma_prior;
  vector[2] period_prior;
  vector[2] Q0_prior;
  vector[2] dQ_prior;
  vector[2] f_prior;
  vector[2] err_prior; // precision for the white noise
  // hyperpara
  vector[N] diag;
}

transformed data {
  real eps = 1e-9;
}

parameters{
   vector[N] eta;
   real<lower = sigma_prior[1], upper = sigma_prior[2]> lsigma;
   real<lower = period_prior[1], upper = period_prior[2]> lperiod;
   real<lower = Q0_prior[1], upper = Q0_prior[2]> lQ0;
   real<lower = dQ_prior[1], upper = dQ_prior[2]> ldQ;
   real<lower = f_prior[1], upper = f_prior[2]> f;
   real<lower = 0> err; //error variance
}

transformed parameters{
   real sigma;
   real period;
   real Q0;
   real dQ;
   sigma = exp(lsigma);
   period = exp(lperiod);
   Q0 = exp(lQ0);
   dQ = exp(ldQ);
}

model{
   vector[N] trend;
   trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
   err ~ inv_gamma(err_prior[1], err_prior[2]);
   y ~ normal(trend, err);
   eta ~ normal(0,1);
                  
   //target += logLikRotation(t, y, sigma, period, 
   //                     Q0, dQ, f, eps, diag + err);

}

generated quantities {
   vector[N] trend;
   trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
}
