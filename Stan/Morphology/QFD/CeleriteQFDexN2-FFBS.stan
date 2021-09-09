// Celerite with Linear Firing Exponential Decay
// Quiet-Firing-Decay (QFD model)
functions {
    real logLikRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
    vector dotCholRotation(vector t, vector y,real sigma, real period, real Q0, real dQ, real f,real eps, vector diag);
    real sigmoid(real x){
        return(1/(exp(-x)+1));
    }
    vector normalize(vector x){
        return(x / sum(x));
    }
}

data {
    int<lower=1> N;
    vector[N] t; // time
    vector[N] y; // light curve
    //log uniform priors for celerite
    vector[2] sigma_prior;
    vector[2] period_prior;
    vector[2] Q0_prior;
    vector[2] dQ_prior;
    vector[2] f_prior;
    // prior for random firing
    //  transitioning
    // quiet can only goto iself or firing 
    vector<lower = 0>[2] alpha_quiet;
    // firing can only goto itself or decay
    vector<lower = 0>[2] alpha_firing;
    // decay can go anywhere
    vector<lower = 0>[3] alpha_decay;
    //   prior for quite noise
    real mu0_quiet;
    real lambda_quiet;
    vector[2] gamma_noise; // shape_rate for noise
    // prior on linear increasing, a normal random walk, with positive offset and slope 1
    real mu0_rate_firing; // prior on how much increase
    real<lower=0> sigma_rate_firing; // prior on how much increase
    // prior on decreasing, a AR with no offset and slope from 0 to 1
    real mu0_rate_decay;// I will put a beta prior here
    real<lower = 0> sigma_rate_decay;
    
    vector[N] diag;// hyperpara, usually set to 0
}

transformed data {
    real eps = 1e-9;
}

parameters {
    // trend parameter
    vector[N] eta;
    real<lower = sigma_prior[1], upper = sigma_prior[2]> lsigma;
    real<lower = period_prior[1], upper = period_prior[2]> lperiod;
    real<lower = Q0_prior[1], upper = Q0_prior[2]> lQ0;
    real<lower = dQ_prior[1], upper = dQ_prior[2]> ldQ;
    real<lower = f_prior[1], upper = f_prior[2]> f;
    // quiet parameter
    simplex[2] theta_quiet; // transitioning probability, 1. to quiet, 2. to firing
    real mu_quiet; // work as grand mu
    real<lower = 0> sigma2_noise; // quiet state variance
    // firing parameter
    simplex[2] theta_firing;// 1. to firing, 2. to decay
    real lograte_firing; // must be increasing on average
    // decay parameter
    simplex[3] theta_decay;// 1. to rest, 2. to firing, 3. to decay
    real lograte_decay; 
}

transformed parameters{
    real sigma;
    real period;
    real Q0;
    real dQ;
    real rate_firing;
    real rate_decay;
    vector[N] trend;
    vector[N] yd;// detrended curve
    vector[3] gamma[N-1];
    sigma = exp(lsigma);
    period = exp(lperiod);
    Q0 = exp(lQ0);
    dQ = exp(ldQ);
    rate_firing = exp(lograte_firing);
    rate_decay = exp(lograte_decay);
    trend = dotCholRotation(t, eta, sigma, period, 
                          Q0, dQ, f, eps, diag);
    yd = y - trend - mu_quiet; // detrended light curve, also minus the usual mean
    // forward algo
    {
        real accu_quiet[2];// used for quiet
        real accu_firing[3];// used for firing, just to avoid -Inf firing can come from anywhere
        real accu_decay[2]; // used for decay
        // The forward algorithm, keep in mind we are going to start with time point 2 
        //   in stead of 1 since we have an AR(1) firing and decaying model

        // quiet state
        gamma[1,1] = normal_lpdf(yd[2]|0, sqrt(sigma2_noise));
        // firing state, a (positive tranded) Gaussian random walk 
        gamma[1,2] = exp_mod_normal_lpdf(yd[2] | yd[1] , sqrt(sigma2_noise) , rate_firing );
        // decay state, exponential decay
        gamma[1,3] = exp_mod_normal_lpdf(yd[2] | 0 , sqrt(sigma2_noise) , rate_decay );

        // then the forward algorithm will start from 3
        for(tt in 2:(N-1)){
            // at quiet state
            //  came from quiet state
            accu_quiet[1] = gamma[tt-1,1] + log(theta_quiet[1]) + 
                      normal_lpdf(yd[tt+1]|0, sqrt(sigma2_noise));
            //  came from decay state
            accu_quiet[2] = gamma[tt-1,3] + log(theta_decay[1]) + 
                      normal_lpdf(yd[tt+1]|0, sqrt(sigma2_noise));
        
            gamma[tt,1] = log_sum_exp(accu_quiet);

            // at firing state
            //  came from quiet state
            accu_firing[1] = gamma[tt-1,1] + log(theta_quiet[2]) + 
                      exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
            //  came from firing state
            accu_firing[2] = gamma[tt-1,2] + log(theta_firing[1]) + 
                      exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
            //  came from decay state, i.e. compound flaring
            accu_firing[3] = gamma[tt-1,3] + log(theta_decay[2]) + 
                      exp_mod_normal_lpdf(yd[tt+1] | yd[tt], sqrt(sigma2_noise) , rate_firing );
            gamma[tt,2] = log_sum_exp(accu_firing);  

            //  at decaying state
            //    came from firing state
            accu_decay[1] = gamma[tt-1, 2] + log(theta_firing[2]) + 
                      exp_mod_normal_lpdf(yd[tt+1] | 0, sqrt(sigma2_noise) , rate_decay );
            //    came from decay state
            accu_decay[2] = gamma[tt-1, 3] + log(theta_decay[3]) + 
                      exp_mod_normal_lpdf(yd[tt+1] | 0, sqrt(sigma2_noise) , rate_decay );
            gamma[tt,3] = log_sum_exp(accu_decay);

        }
    }
    
}

model{ 
    // get the trend, chol parameterization
    
    eta ~ normal(0,1);
    // prior settings 
      // No need of trend GP model parameters since coded in parameter section 
      // Firing HMM model parameters
        // AR of states
    sigma2_noise ~ inv_gamma(gamma_noise[1], gamma_noise[2]);
        // mean parameters
    mu_quiet ~ normal(mu0_quiet, sqrt(sigma2_noise/lambda_quiet));// serve as overall mean 
    //increm_firing ~ exponential(alpha_incre_firing); // the on average increase when firing
    //rate_decay ~ beta(alpha_rate_decay[1], alpha_rate_decay[2]);
    lograte_firing ~ normal(mu0_rate_firing, sigma_rate_firing);
    lograte_decay ~ normal(mu0_rate_decay, sigma_rate_decay);
        // transition 
    theta_quiet ~ dirichlet(alpha_quiet);
    theta_firing ~ dirichlet(alpha_firing);
    theta_decay ~ dirichlet(alpha_decay);
    
    // likelihood part
    
    // HMM firing with linear increase and exponential decrease AR structure
    
    target += log_sum_exp(gamma[N-1]);
}

generated quantities {
    int<lower = 1, upper = 3> FFBS[N-1]; 
    real Q1;
    real w1;
    real S1;
    real Q2;
    real w2;
    real S2;
    real sigma1;
    real rho1;
    real tau1;
    real sigma2;
    real rho2;
    real tau2;
    real amp;

    // parameters for celerite
    {
        amp =  sigma * sigma/ (1 + f);
        Q1 = 0.5 + Q0 + dQ;
        w1 = 4 * 3.1415926 * Q1 / (period * sqrt(4 * Q1 * Q1 - 1));
        S1 = amp / (w1 * Q1);
        Q2 = 0.5 + Q0;
        w2 = 8 * 3.1415926 * Q2 / (period * sqrt(4 * Q2 * Q2 - 1));
        S2 = f * amp / (w2 * Q2);
    
        rho1 = 2*3.1415926/w1;
        rho2 = 2*3.1415926/w2;
        tau1 = 2 * Q1/w1;
        tau2 = 2 * Q2/w2;
        sigma1 = sqrt(S1 * w1 * Q1);
        sigma2 = sqrt(S2 * w2 * Q2);
    }
     
    // FFBS state decoding:
    {
        {
            vector[2] p_temp2;
            vector[3] p_temp3;
            FFBS[N-1] = categorical_logit_rng(gamma[N-1]);
            for(tforward in 0:(N-1-2)){
                int tt = 0;
                tt = N-1-tforward; // backward...

                // next is quiet, previsou can only be quiet and decay
                if(FFBS[tt]==1){
                    // is quiet
                    p_temp2[1] = gamma[tt-1,1] + log(theta_quiet[1]);
                    // is decay
                    p_temp2[2] = gamma[tt-1,3] + log(theta_decay[1]);
                    FFBS[tt-1] = categorical_logit_rng(p_temp2) * 2 - 1; // this make 1 to 1 and 2 to 3
                }
                // next is firing, can be any one of it
                if(FFBS[tt]==2){
                    // is quiet:
                    p_temp3[1] = gamma[tt-1,1] + log(theta_quiet[2]);
                    // is firing
                    p_temp3[2] = gamma[tt-1,2] + log(theta_firing[1]);
                    // is decay
                    p_temp3[3] = gamma[tt-1,3] + log(theta_decay[2]);
                    FFBS[tt-1] = categorical_logit_rng(p_temp3);
                }
                // next is decay, can come from firing or decay
                if(FFBS[tt]==3){
                    // from firing
                    p_temp2[1] = gamma[tt-1,2] + log(theta_firing[2]);
                    // from decay
                    p_temp2[2] = gamma[tt-1,3] + log(theta_decay[3]);
                    FFBS[tt-1] = categorical_logit_rng(p_temp2) + 1;
                }
                
            }// backward filtering loops

        }// end of backward

    }
}
