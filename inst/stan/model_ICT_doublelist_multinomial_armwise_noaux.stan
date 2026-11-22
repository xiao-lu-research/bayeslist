// Separate proportional-odds models for each list and treatment/control arm.
// With intercept-only covariates, each arm spans all positive count probabilities.
// Additional covariates impose proportional odds within each list and arm.
// No cross-arm constraints enforce a common control count plus a binary item.
functions{
 // Return E[Y + 1]; ordered_logistic uses one-based categories for zero-based counts.
  real ordered_logistic_expectation(real eta, vector c) {
    int K = num_elements(c) + 1;
    vector[K] probs;
    vector[K] cum_probs_upper;

    // Calculate P(Y + 1 <= k) for k = 1,...,K-1.
    for (k in 1:(K-1)) {
      cum_probs_upper[k] = inv_logit(c[k] - eta);
    }
    cum_probs_upper[K] = 1.0; // P(Y + 1 <= K) is 1.

    // Calculate the probabilities of the one-based categories.
    probs[1] = cum_probs_upper[1]; // P(Y = 0).
    for (k in 2:K) {
      probs[k] = cum_probs_upper[k] - cum_probs_upper[k-1];
    }

    // The common +1 offset cancels in each treatment-minus-control contrast.
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
	vector[N] weights; // Weights for the mean-count contrast, not the likelihood.
	array[K] real mu_psi0;
	array[K] real mu_delta;
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
        array[J-1] real<lower = 0> sigma_psiM; // Prior SDs for control log cutpoint gaps.
        array[J] real<lower = 0> sigma_psiMT; // Prior SDs for treatment log cutpoint gaps.
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
}

parameters {
	vector[K] psi0; // List 1 control-count regression coefficients.
	vector[K] delta; // List 1 treatment-count coefficients, not sensitive-item log odds.
	vector[K] psi02; // List 2 control-count regression coefficients.
	vector[K] delta2; // List 2 treatment-count coefficients, not sensitive-item log odds.
	vector[J-1] psiM; // Unrestricted log cutpoint gaps for list 1 control counts.
	vector[J-1] psiM2; // Unrestricted log cutpoint gaps for list 2 control counts.
	vector[J] psiMT; // Unrestricted log cutpoint gaps for list 1 treatment counts.
	vector[J] psiMT2; // Unrestricted log cutpoint gaps for list 2 treatment counts.
}


transformed parameters {
  // Fix the first cutpoint at zero to identify the regression intercept.
  // Cumulative exp(log gaps) orders the cutpoints; the log gaps need no ordering.
  ordered[J] opsiM; // Ordered cutpoints for list 1 control counts.
  ordered[J] opsiM2; // Ordered cutpoints for list 2 control counts.
  ordered[J+1] opsiMT; // Ordered cutpoints for list 1 treatment counts.
  ordered[J+1] opsiMT2; // Ordered cutpoints for list 2 treatment counts.
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

  real pas; // Weighted average of the two mean-count contrasts; may lie outside [0,1].
  pas = 0;
  for (i in 1:N){
    pas += weights[i] * (
        0.5*(ordered_logistic_expectation(X[i,] * delta, opsiMT) - ordered_logistic_expectation(X[i,] * psi0, opsiM))
      + 0.5*(ordered_logistic_expectation(X2[i,] * delta2, opsiMT2) - ordered_logistic_expectation(X2[i,] * psi02, opsiM2)));
  }
  pas /= sum(weights);
}

model{

        delta ~ normal(mu_delta,sigma_delta); // Normal prior with supplied mean and SD.
        psi0  ~ normal(mu_psi0,sigma_psi0); // Normal prior with supplied mean and SD.
        delta2 ~ normal(mu_delta,sigma_delta); // Normal prior with supplied mean and SD.
        psi02  ~ normal(mu_psi0,sigma_psi0); // Normal prior with supplied mean and SD.

	psiM  ~ normal(0,sigma_psiM); // Normal prior on unrestricted log gaps for list 1 control counts.
	psiM2 ~ normal(0,sigma_psiM); // Normal prior on unrestricted log gaps for list 2 control counts.
	psiMT ~ normal(0,sigma_psiMT); // Normal prior on unrestricted log gaps for list 1 treatment counts.
	psiMT2 ~ normal(0,sigma_psiMT); // Normal prior on unrestricted log gaps for list 2 treatment counts.
	
	for (i in 1:N){
          // First-list contribution, using its observed treatment assignment.
	  if (treat[i] == 1){
	    target += ordered_logistic_lupmf(Y1[i]+1 | X[i,] * delta, opsiMT);
	  } else { // no treatment
	    target += ordered_logistic_lupmf(Y1[i]+1 | X[i,] * psi0, opsiM);
	  }
          // Second-list contribution, using its observed treatment assignment.
	  if (treat2[i] == 1){
	    target += ordered_logistic_lupmf(Y2[i]+1 | X2[i,] * delta2, opsiMT2);
	  } else { // no second treatment
	    target += ordered_logistic_lupmf(Y2[i]+1 | X2[i,] * psi02, opsiM2);
	  }
	}

}

generated quantities {
    vector[N] log_lik; // Joint pointwise log likelihood, summing both list contributions.
	for (i in 1:N){
          // First-list contribution, using its observed treatment assignment.
	  if (treat[i] == 1){
	    log_lik[i] = ordered_logistic_lpmf(Y1[i]+1 | X[i,] * delta, opsiMT);
	  } else { // no treatment
	    log_lik[i] = ordered_logistic_lpmf(Y1[i]+1 | X[i,] * psi0, opsiM);
	  }
        // Second-list contribution, using its observed treatment assignment.
	  if (treat2[i] == 1){
	    log_lik[i] += ordered_logistic_lpmf(Y2[i]+1 | X2[i,] * delta2, opsiMT2);
	  } else { // no second treatment
	    log_lik[i] += ordered_logistic_lpmf(Y2[i]+1 | X2[i,] * psi02, opsiM2);
	  }
	}
} 
