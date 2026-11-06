functions{
	
}

data {
	int N;
	int J; //number of non sensitive item
	array[N] int<lower = 0> Y1; // number of affirmative answers to the first list question
	array[N] int<lower = 0> Y2; // number of affirmative answers to the second list question
	int<lower=1, upper=1> K; // this model has only an intercept
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
  real pas; // average difference in expected list counts
  pas = 0.5 * ((J + 1) * inv_logit(delta[1]) - J * inv_logit(psi0[1]))
      + 0.5 * ((J + 1) * inv_logit(delta2[1]) - J * inv_logit(psi02[1]));
}

model{

        delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
        delta2 ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi02  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
	
	for (i in 1:N){
// First those that received ttt 1 (and not ttt 2)
		if (treat[i] == 1){ // 0 < Y1[i] < J+1
			target += binomial_lpmf(Y1[i]|J+1,inv_logit(delta[1]));
		} else {
			target += binomial_lpmf(Y1[i]|J,inv_logit(psi0[1]));
		}

// Now those that received ttt 2 (and not ttt 1)

		if (treat2[i] == 1){
			target += binomial_lpmf(Y2[i]|J+1,inv_logit(delta2[1]));
		} else {
			target += binomial_lpmf(Y2[i]|J,inv_logit(psi02[1]));
		}
	}

}

generated quantities {
    vector[N] log_lik; 
    for (i in 1:N){
// First those that received ttt 1 (and not ttt 2)
		if (treat[i] == 1){ // 0 < Y1[i] < J+1
			log_lik[i] = binomial_lpmf(Y1[i]|J+1,inv_logit(delta[1]));
		} else {
			log_lik[i] = binomial_lpmf(Y1[i]|J,inv_logit(psi0[1]));
		}

// Now those that received ttt 2 (and not ttt 1)

		if (treat2[i] == 1){
			log_lik[i] += binomial_lpmf(Y2[i]|J+1,inv_logit(delta2[1]));
		} else {
			log_lik[i] += binomial_lpmf(Y2[i]|J,inv_logit(psi02[1]));
		}
	}

}
