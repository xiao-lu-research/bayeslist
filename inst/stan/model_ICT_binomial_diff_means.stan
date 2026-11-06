functions{

}

data {
	int N;
	int J; //number of non sensitive item
	array[N] int<lower = 0> Y; // number of affirmative answers
	int K; // The number of control predictors
	array[N] int treat;
	vector[N] weights; // survey weights
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
}

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	vector[K] delta; // coefficients of controls on Z
}

transformed parameters {
  real<lower=0, upper=1> pas;
  pas = (J+1)*inv_logit(delta[1])-J*inv_logit(psi0[1]);
//pas = inv_logit(delta[1]);
}

model{

	psi0  ~ normal(0,sigma_psi0); // prior for the coefficients of controls (with or without treatment)
  	delta ~ normal(0,sigma_delta); // prior for the coefficients of controls on Z

	for (i in 1:N){
 //       if (treat[i] == 1 && Y[i] == 0) {
 //         target += log(1 - inv_logit(delta[1])) + binomial_lpmf(0|J,inv_logit(psi0[1]));
 //       } else if (treat[i] == 1 && Y[i] == J + 1){
 //         target += log(inv_logit(delta[1])) + binomial_lpmf(J|J,inv_logit(psi0[1]));
 //       } else if (treat[i] == 1){
          if (treat[i] == 1){
	    target += binomial_lpmf(Y[i]|J+1,inv_logit(delta[1]));
	  } else { // no treatment
	    target += binomial_lpmf(Y[i]|J,inv_logit(psi0[1]));
	  }
	}

}

generated quantities {
    vector[N] log_lik;
    for (i in 1:N){
	  if (treat[i] == 1) {
	    log_lik[i] = binomial_lpmf(Y[i]|J+1,inv_logit(delta[1]));
	  } else { // no treatment
	    log_lik[i] = binomial_lpmf(Y[i]|J,inv_logit(psi0[1]));
	  }
    }
} 
