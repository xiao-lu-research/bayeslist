// Separate proportional-odds models for observed treatment and control counts.
// With an intercept-only X, each arm spans all positive category probabilities.
// With additional covariates, proportional odds restricts effects within each arm.
// No cross-arm constraints enforce a common latent control count plus binary item.
functions{
 // Return E[Y + 1]: ordered_logistic uses categories 1,...,K for counts 0,...,K-1.
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

    // The common +1 offset cancels in the treatment-minus-control contrast.
    return dot_product(probs, linspaced_vector(K, 1.0, 1.0*K));
  }
}

data {
	int N;
	int J; // Number of non-sensitive items; treatment lists contain J+1 items.
	array[N] int<lower = 0> Y; // Observed count, from 0 to J+treat[i].
	int K; // Number of predictors, including the intercept, in each arm's model.
	matrix[N,K] X; // Covariates used for both arm-specific regressions.
	array[N] int treat; // 1 for the treatment list; 0 for the control list.
	vector[N] weights; // Weights for the standardized mean contrast, not the likelihood.
        array[K] real<lower = 0> sigma_delta; // Prior SDs for treatment-count coefficients.
        array[J-1] real<lower = 0> sigma_psiM; // Prior SDs for control log cutpoint gaps.
        array[J] real<lower = 0> sigma_psiMT; // Prior SDs for treatment log cutpoint gaps.
        array[K] real<lower = 0> sigma_psi0; // Prior SDs for control-count coefficients.
}

parameters {
	vector[J-1] psiM; // Unrestricted log cutpoint gaps for control counts 0,...,J.
	vector[J] psiMT; // Unrestricted log cutpoint gaps for treatment counts 0,...,J+1.
	vector[K] psi0; // Control-count regression coefficients.
	vector[K] delta; // Treatment-count coefficients, not sensitive-item log odds.
}

transformed parameters {
  // Fix the first cutpoint in each arm to identify its regression intercept.
  // Positive increments order the cutpoints without ordering their gaps.
  ordered[J] opsiM; // Cutpoints for the J+1 control-count categories.
  opsiM[1] = 0;
  for (k in 2:J){
    opsiM[k] = opsiM[k-1] + exp(psiM[k-1]);
  }
  ordered[J+1] opsiMT; // Cutpoints for the J+2 treatment-count categories.
  opsiMT[1] = 0;
  for (k in 2:(J+1)){
    opsiMT[k] = opsiMT[k-1] + exp(psiMT[k-1]);
  }

  // Standardize both expected counts over the same weighted covariate distribution.
  // This contrast estimates prevalence under the ICT assumptions, but posterior
  // draws can be outside [0,1] because no ICT compatibility constraints couple the arms.
  real pas;
  pas = 0;
  for (i in 1:N){
    pas += weights[i] * (ordered_logistic_expectation(X[i,] * delta, opsiMT) - ordered_logistic_expectation(X[i,] * psi0, opsiM));
  }
  pas /= sum(weights);
}

model{

        // Independent zero-centered normal priors on coefficients and log gaps.
        // These induce priors on the cell probabilities; no prior truncates pas.
	psiM  ~ normal(0,sigma_psiM);
	psiMT ~ normal(0,sigma_psiMT);
	psi0  ~ normal(0,sigma_psi0);
	delta ~ normal(0,sigma_delta);

	// Independent respondent likelihood contributions, with separate arm parameters.
	// Add one to the observed count for Stan's one-based ordered-logistic categories.
	for (i in 1:N){
	  if (treat[i] == 1){
	    target += ordered_logistic_lupmf(Y[i]+1 | X[i,] * delta, opsiMT);
	  } else { // no treatment
	    target += ordered_logistic_lupmf(Y[i]+1 | X[i,] * psi0, opsiM);
	  }
	}

}

generated quantities {
    vector[N] log_lik; // Normalized pointwise log likelihood for LOO/WAIC.
    for (i in 1:N){
      if (treat[i] == 1){
        log_lik[i] = ordered_logistic_lpmf(Y[i]+1 | X[i,] * delta, opsiMT);
      } else { // no treatment
        log_lik[i] = ordered_logistic_lpmf(Y[i]+1 | X[i,] * psi0, opsiM);
      }
    }
} 
