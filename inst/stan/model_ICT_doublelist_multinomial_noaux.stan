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
	vector[N] weights; // survey weights
	array[K] real mu_psi0;
	array[K] real mu_delta;
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
        array[J-1] real<lower = 0> sigma_psiM; // prior psiM standard deviation
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
}

parameters {
	vector[K] psi0; // coefficients of first list controls without treatment 1
	vector[K] delta; // coefficients of first list controls on Z
	vector[K] psi02; // coefficients of controls without treatment 2
	vector[K] delta2; // coefficients of second list controls on Z
	ordered[J-1] psiM; // coefficients of controls without treatment
	ordered[J-1] psiM2; // coefficients of controls without treatment
}


transformed parameters {
  real<lower=0, upper=1> pas; // proportion affirmative on sensitive items
  pas = dot_product(weights, 0.5*(inv_logit(X * delta))+0.5*(inv_logit(X2 * delta2)));

  ordered[J] opsiM; // coefficients of controls without treatment
  ordered[J] opsiM2; // coefficients of controls without treatment
//opsiM[1] = psi0[1];
  opsiM[1] = 0;
  opsiM2[1] = 0;
  for (k in 2:J){
//  opsiM[k] = psiM[k-1] + psi0[1];
    opsiM[k] = opsiM[k-1] + exp(psiM[k-1]);
    opsiM2[k] = opsiM2[k-1] + exp(psiM2[k-1]);
  }
  vector<lower=0, upper=1>[N] pz; // 
  vector<lower=0, upper=1>[N] pz2; // 
  for (i in 1:N){
    pz[i]  = inv_logit(X[i,] * delta);
    pz2[i] = inv_logit(X2[i,] * delta2);
  }
}

model{

        delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
        delta2 ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi02  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment

	psiM  ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)
	psiM2 ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)
	
	for (i in 1:N){
          // First those that received ttt 1 (and not ttt 2)
	  if (treat[i] == 1 && Y1[i] == 0) {
	    target += log(1 - pz[i]) + ordered_logistic_lupmf( 1 | X[i,] * psi0, opsiM);
	  } else if (treat[i] == 1 && Y1[i] == J + 1){
	    target += log(    pz[i]) + ordered_logistic_lupmf(J+1| X[i,] * psi0, opsiM);
	  } else if (treat[i] == 1){
	    target += log(pz[i]*exp(ordered_logistic_lupmf(Y1[i]   | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lupmf(Y1[i]+1| X[i,] * psi0, opsiM)) );
	  } else { // no treatment
	    target += log(pz[i]*exp(ordered_logistic_lupmf(Y1[i]+1 | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lupmf(Y1[i]+1 | X[i,] * psi0, opsiM)) );
	  }
          // Now those that received ttt 2 (and not ttt 1)
	  if (treat2[i] == 1 && Y2[i] == 0) {
	    target += log(1 - pz2[i]) + ordered_logistic_lupmf( 1 | X2[i,] * psi02, opsiM2);
	  } else if (treat2[i] == 1 && Y2[i] == J + 1){
	    target += log(    pz2[i]) + ordered_logistic_lupmf(J+1| X2[i,] * psi02, opsiM2);
	  } else if (treat2[i] == 1){
	    target += log(pz2[i]*exp(ordered_logistic_lupmf(Y2[i]   | X2[i,] * psi02, opsiM2)) + (1 - pz2[i]) * exp(ordered_logistic_lupmf(Y2[i]+1| X2[i,] * psi02, opsiM2)) );
	  } else { // no second treatment
	    target += log(pz2[i]*exp(ordered_logistic_lupmf(Y2[i]+1 | X2[i,] * psi02, opsiM2)) + (1 - pz2[i]) * exp(ordered_logistic_lupmf(Y2[i]+1 | X2[i,] * psi02, opsiM2)) );
	  }
	}

}

generated quantities {
    vector[N] log_lik;
	for (i in 1:N){
          // First those that received ttt 1 (and not ttt 2)
	  if (treat[i] == 1 && Y1[i] == 0) {
	    log_lik[i] = log(1 - pz[i]) + ordered_logistic_lpmf( 1 | X[i,] * psi0, opsiM);
	  } else if (treat[i] == 1 && Y1[i] == J + 1){
	    log_lik[i] = log(    pz[i]) + ordered_logistic_lpmf(J+1| X[i,] * psi0, opsiM);
	  } else if (treat[i] == 1){
	    log_lik[i] = log(pz[i]*exp(ordered_logistic_lpmf(Y1[i]   | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lpmf(Y1[i]+1| X[i,] * psi0, opsiM)) );
	  } else { // no treatment
	    log_lik[i] = log(pz[i]*exp(ordered_logistic_lpmf(Y1[i]+1 | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lpmf(Y1[i]+1 | X[i,] * psi0, opsiM)) );
	  }
        // Now those that received ttt 2 (and not ttt 1)
	  if (treat2[i] == 1 && Y2[i] == 0) {
	    log_lik[i] += log(1 - pz2[i]) + ordered_logistic_lpmf( 1 | X2[i,] * psi02, opsiM2);
	  } else if (treat2[i] == 1 && Y2[i] == J + 1){
	    log_lik[i] += log(    pz2[i]) + ordered_logistic_lpmf(J+1| X2[i,] * psi02, opsiM2);
	  } else if (treat2[i] == 1){
	    log_lik[i] += log(pz2[i]*exp(ordered_logistic_lpmf(Y2[i]   | X2[i,] * psi02, opsiM2)) + (1 - pz2[i]) * exp(ordered_logistic_lpmf(Y2[i]+1| X2[i,] * psi02, opsiM2)) );
	  } else { // no second treatment
	    log_lik[i] += log(pz2[i]*exp(ordered_logistic_lpmf(Y2[i]+1 | X2[i,] * psi02, opsiM2)) + (1 - pz2[i]) * exp(ordered_logistic_lpmf(Y2[i]+1 | X2[i,] * psi02, opsiM2)) );
	  }
	}
} 
