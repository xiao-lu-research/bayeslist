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
	array[N] int<lower = 0> Y; // number of affirmative answers
	int K; // The number of control predictors
	matrix[N,K] X;
	array[N] int treat;
	vector[N] weights; // survey weights
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
        array[J-1] real<lower = 0> sigma_psiM; // prior psiM standard deviation
        array[J] real<lower = 0> sigma_psiMT; // prior psiM standard deviation
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
}

parameters {
	ordered[J-1] psiM; // coefficients of controls without treatment
	ordered[J] psiMT; // coefficients of controls without treatment
	vector[K] psi0; // coefficients of controls without treatment
	vector[K] delta; // coefficients of controls on Z
}

transformed parameters {
  ordered[J] opsiM; // coefficients of controls without treatment
  opsiM[1] = 0;
  for (k in 2:J){
    opsiM[k] = opsiM[k-1] + exp(psiM[k-1]);
  }
  ordered[J+1] opsiMT; // coefficients of controls without treatment
  opsiMT[1] = 0;
  for (k in 2:(J+1)){
    opsiMT[k] = opsiMT[k-1] + exp(psiMT[k-1]);
  }

//real<lower=0, upper=1> pas;
  real pas; // weighted average difference in expected list counts
  pas = 0;
  for (i in 1:N){
    pas += weights[i] * (ordered_logistic_expectation(X[i,] * delta, opsiMT) - ordered_logistic_expectation(X[i,] * psi0, opsiM));
  }
  pas /= sum(weights);
}

model{

	psiM  ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)
	psiMT ~ normal(0,sigma_psiMT); // prior for the coefficients of controls (with or without treatment)
	psi0  ~ normal(0,sigma_psi0); // prior for the coefficients of controls (with or without treatment)

// One or the other of these priors
  	delta ~ normal(0,sigma_delta); // prior for the coefficients of controls on Z
//	target += log((N+1)/(N-1)) + log(1-(pow(pas, N) + pow(1-pas, N))); // ND prior on pas

	for (i in 1:N){
	  if (treat[i] == 1){
	    target += ordered_logistic_lupmf(Y[i]+1 | X[i,] * delta, opsiMT);
	  } else { // no treatment
	    target += ordered_logistic_lupmf(Y[i]+1 | X[i,] * psi0, opsiM);
	  }
	}

}

generated quantities {
    vector[N] log_lik;
    for (i in 1:N){
      if (treat[i] == 1){
        log_lik[i] = ordered_logistic_lpmf(Y[i]+1 | X[i,] * delta, opsiMT);
      } else { // no treatment
        log_lik[i] = ordered_logistic_lpmf(Y[i]+1 | X[i,] * psi0, opsiM);
      }
    }
} 
