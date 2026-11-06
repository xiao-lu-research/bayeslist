functions{
	
}

data {
	int N;
	int J; //number of non sensitive item
	array[N] int<lower = 0> Y1; // number of affirmative answers to the first list question
	array[N] int<lower = 0> Y2; // number of affirmative answers to the second list question
	int K;
        matrix[N,K] X;
        matrix[N,K] X2;
	array[N] int treat;
	array[N] int treat2;
        array[K] real mu_psi0;
        array[K] real mu_delta;
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
	vector[N] weights; // survey weights
}

parameters {
	vector[K] psi0; // coefficients of first list controls without treatment 1
	vector[K] delta; // coefficients of first list controls on Z
	vector[K] psi02; // coefficients of controls without treatment 2
	vector[K] delta2; // coefficients of second list controls on Z
}


transformed parameters {
  real<lower=0, upper=1> pas; // proportion affirmative on sensitive items
  pas = dot_product(weights, 0.5*(inv_logit(X * delta))+0.5*(inv_logit(X2 * delta2)));

  vector<lower=0, upper=1>[N] pz; //
  vector<lower=0, upper=1>[N] pz2; //
  vector<lower=0>[N] il; //
  vector<lower=0>[N] il2; //
  for (i in 1:N){
    pz[i]  = inv_logit(X[i,] * delta);
    pz2[i] = inv_logit(X2[i,] * delta2);
    il[i] = inv_logit(X[i,] * psi0);
    il2[i] = inv_logit(X2[i,] * psi02);
  }
}

model{

        delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
        delta2 ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi02  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
	
	for (i in 1:N){
          // First those that received ttt 1 (and not ttt 2)
	  if (treat[i] == 1 && Y1[i] == 0) {
	    target += log(1 - pz[i]) + binomial_lpmf( 0 | J, il[i]);
	  } else if (treat[i] == 1 && Y1[i] == J + 1){
	    target += log(    pz[i]) + binomial_lpmf(J | J,il[i]);
	  } else if (treat[i] == 1){
            target += log(pz[i] * exp(binomial_lpmf(Y1[i] - 1 | J,il[i])) + (1 - pz[i]) * exp(binomial_lpmf(Y1[i] | J,il[i])));
          } else { // no treatment
            target += log(pz[i] * exp(binomial_lpmf(Y1[i] | J,il[i])) + (1 - pz[i]) * exp(binomial_lpmf(Y1[i] | J,il[i])) );
          }
          // Now those that received ttt 2 (and not ttt 1)
	  if (treat2[i] == 1 && Y2[i] == 0) {
	    target += log(1 - pz2[i]) + binomial_lpmf( 0 | J, il2[i]);
	  } else if (treat2[i] == 1 && Y2[i] == J + 1){
	    target += log(    pz2[i]) + binomial_lpmf(J | J,il2[i]);
	  } else if (treat2[i] == 1){
            target += log(pz2[i] * exp(binomial_lpmf(Y2[i] - 1 | J,il2[i])) + (1 - pz2[i]) * exp(binomial_lpmf(Y2[i] | J,il2[i])));
	  } else { // no second treatment
            target += log(pz2[i] * exp(binomial_lpmf(Y2[i] | J,il2[i])) + (1 - pz2[i]) * exp(binomial_lpmf(Y2[i] | J,il2[i])) );
	  }
	}

}

generated quantities {
    vector[N] log_lik; 
    for (i in 1:N){
// First those that received ttt 1 (and not ttt 2)
	  if (treat[i] == 1 && Y1[i] == 0) {
	    log_lik[i] = log(1 - pz[i]) + binomial_lpmf( 0 | J, il[i]);
	  } else if (treat[i] == 1 && Y1[i] == J + 1){
	    log_lik[i] = log(    pz[i]) + binomial_lpmf(J | J,il[i]);
	  } else if (treat[i] == 1){
            log_lik[i] = log(pz[i] * exp(binomial_lpmf(Y1[i] - 1 | J,il[i])) + (1 - pz[i]) * exp(binomial_lpmf(Y1[i] | J,il[i])));
          } else { // no treatment
            log_lik[i] = log(pz[i] * exp(binomial_lpmf(Y1[i] | J,il[i])) + (1 - pz[i]) * exp(binomial_lpmf(Y1[i] | J,il[i])) );
          }
          // Now those that received ttt 2 (and not ttt 1)
	  if (treat2[i] == 1 && Y2[i] == 0) {
	    log_lik[i] += log(1 - pz2[i]) + binomial_lpmf( 0 | J, il2[i]);
	  } else if (treat2[i] == 1 && Y2[i] == J + 1){
	    log_lik[i] += log(    pz2[i]) + binomial_lpmf(J | J,il2[i]);
	  } else if (treat2[i] == 1){
            log_lik[i] += log(pz2[i] * exp(binomial_lpmf(Y2[i] - 1 | J,il2[i])) + (1 - pz2[i]) * exp(binomial_lpmf(Y2[i] | J,il2[i])));
	  } else { // no second treatment
            log_lik[i] += log(pz2[i] * exp(binomial_lpmf(Y2[i] | J,il2[i])) + (1 - pz2[i]) * exp(binomial_lpmf(Y2[i] | J,il2[i])) );
	  }
    }

}
