functions{
 // Function to calculate the expected value of Y for given eta and cutpoints c
  real ordered_logistic_expectation(real eta, vector c) {
    int K = num_elements(c) + 1;
    vector[K] probs;
    vector[K] cum_probs_upper;

    // Calculate P(Y <= k) for k = 1 to K-1
    for (k in 1:(K-1)) {
      cum_probs_upper[k] = inv_logit(c[k] - eta);
    }
    cum_probs_upper[K] = 1.0; // P(Y <= K) is 1

    // Calculate individual category probabilities P(Y = k)
    probs[1] = cum_probs_upper[1]; // P(Y=1) = P(Y<=1) - P(Y<=0) = P(Y<=1) - 0
    for (k in 2:K) {
      probs[k] = cum_probs_upper[k] - cum_probs_upper[k-1];
    }

    // Compute the expectation E[Y] = sum(k * P(Y=k))
    return dot_product(probs, linspaced_vector(K, 1.0, 1.0*K));
  }
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
        array[J] real<lower = 0> sigma_psiMT; // prior psiM standard deviation
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
}

parameters {
	vector[K] psi0; // coefficients of first list controls without treatment 1
	vector[K] delta; // coefficients of first list controls on Z
	vector[K] psi02; // coefficients of controls without treatment 2
	vector[K] delta2; // coefficients of second list controls on Z
	ordered[J-1] psiM; // coefficients of controls without treatment
	ordered[J-1] psiM2; // coefficients of controls without treatment
	ordered[J] psiMT; // coefficients of controls without treatment
	ordered[J] psiMT2; // coefficients of controls without treatment
}


transformed parameters {
  ordered[J] opsiM;  // coefficients of controls without treatment
  ordered[J] opsiM2; // coefficients of controls without treatment
  ordered[J+1] opsiMT;  // coefficients of controls without treatment
  ordered[J+1] opsiMT2; // coefficients of controls without treatment
  opsiM[1] = 0;
  opsiM2[1] = 0;
  for (k in 2:J){
    opsiM[k] = opsiM[k-1] + exp(psiM[k-1]);
    opsiM2[k] = opsiM2[k-1] + exp(psiM2[k-1]);
  }
  opsiMT[1] = 0;
  opsiMT2[1] = 0;
  for (k in 2:(J+1)){
    opsiMT[k] = opsiMT[k-1] + exp(psiMT[k-1]);
    opsiMT2[k] = opsiMT2[k-1] + exp(psiMT2[k-1]);
  }

  real pas; // weighted average difference in expected list counts
  pas = 0;
  for (i in 1:N){
    pas += weights[i] * (
        0.5*(ordered_logistic_expectation(X[i,] * delta, opsiMT) - ordered_logistic_expectation(X[i,] * psi0, opsiM))
      + 0.5*(ordered_logistic_expectation(X2[i,] * delta2, opsiMT2) - ordered_logistic_expectation(X2[i,] * psi02, opsiM2)));
  }
  pas /= sum(weights);
}

model{

        delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
        delta2 ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi02  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment

	psiM  ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)
	psiM2 ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)
	psiMT ~ normal(0,sigma_psiMT); // prior for the coefficients of controls (with or without treatment)
	psiMT2 ~ normal(0,sigma_psiMT); // prior for the coefficients of controls (with or without treatment)
	
	for (i in 1:N){
          // First those that received ttt 1 (and not ttt 2)
	  if (treat[i] == 1){
	    target += ordered_logistic_lupmf(Y1[i]+1 | X[i,] * delta, opsiMT);
	  } else { // no treatment
	    target += ordered_logistic_lupmf(Y1[i]+1 | X[i,] * psi0, opsiM);
	  }
          // Now those that received ttt 2 (and not ttt 1)
	  if (treat2[i] == 1){
	    target += ordered_logistic_lupmf(Y2[i]+1 | X2[i,] * delta2, opsiMT2);
	  } else { // no second treatment
	    target += ordered_logistic_lupmf(Y2[i]+1 | X2[i,] * psi02, opsiM2);
	  }
	}

}

generated quantities {
    vector[N] log_lik;
	for (i in 1:N){
          // First those that received ttt 1 (and not ttt 2)
	  if (treat[i] == 1){
	    log_lik[i] = ordered_logistic_lpmf(Y1[i]+1 | X[i,] * delta, opsiMT);
	  } else { // no treatment
	    log_lik[i] = ordered_logistic_lpmf(Y1[i]+1 | X[i,] * psi0, opsiM);
	  }
        // Now those that received ttt 2 (and not ttt 1)
	  if (treat2[i] == 1){
	    log_lik[i] += ordered_logistic_lpmf(Y2[i]+1 | X2[i,] * delta2, opsiMT2);
	  } else { // no second treatment
	    log_lik[i] += ordered_logistic_lpmf(Y2[i]+1 | X2[i,] * psi02, opsiM2);
	  }
	}
} 
