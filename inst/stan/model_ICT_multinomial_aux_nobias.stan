functions{
	
	int group_size(array[] int ref, int value) {
    	  int count;
    	  count = 0;
    	  for (ii in 1:size(ref))
      	  if (ref[ii]==value)
        	  count = count + 1;
      	  return count;
  	}

	// returns the subset of y[ref==value,1:k]
  	matrix subset_matrix(matrix y, int k, array[] int ref, int value) {
 	   int jj;
	    matrix[group_size(ref, value),k] res;
	    jj = 1;
	    for(ii in 1:size(ref)) {
	      if (ref[ii] == value) {
	    	for (kk in 1:k){
	    		res[jj,kk] = y[ii,kk];	
	    	}
	        jj = jj+1;
	      }
	  }
	  return res;
	}

 // 	returns the subset of y[ref==value]
   	vector subset_vector(vector y, array[] int ref, int value) {
     	  int jj;
     	  vector[group_size(ref, value)] res;
     	  jj = 1;
     	  for(ii in 1:size(ref)) {
       	    if (ref[ii] == value) {
       		res[jj] = y[ii];	
         	jj = jj+1;
       	    }
     	  }
     	return res;
   	}

}

data {
	int N;
	int J; //number of non sensitive item
	array[N] int<lower = 0> Y; // number of affirmative answers
	int K; // The number of control predictors
	matrix[N,K] X;
	array[N] int treat;
	int G; // number of subgroups
	array[G] real h; // auxiliary information for subgroups
	array[N] int g; // subgroup indicators
	vector[N] weights; // survey weights
	array[G] real sigma_h; // subgroup standard deviations
	array[K] real mu_psi0;
	array[K] real mu_delta;
        array[K] real<lower = 0> sigma_delta; // prior delta standard deviation
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
        array[J-1] real<lower = 0> sigma_psiM; // prior psiM standard deviation
}

transformed data{
  vector<lower=0>[G] a;
  vector<lower=0>[G] b;
  for (gg in 1:G){
    a[gg] = h[gg] * ((h[gg] * (1 - h[gg]) / (sigma_h[gg]*sigma_h[gg])) - 1);
    b[gg] = (1 - h[gg]) * ((h[gg] * (1 - h[gg]) / (sigma_h[gg]*sigma_h[gg])) - 1);
  }
}

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	// vector[K] psi1; // coefficients of controls with treatment
	vector[K] delta; // coefficients of controls on Z
	ordered[J-1] psiM; // coefficients of controls without treatment
}

transformed parameters {
  vector<lower=0, upper=1>[G] pas;
  for (gg in 1:G){
// returns X[g==gg,1:K]*delta
    pas[gg] = dot_product(subset_vector(weights,g,gg),
                (inv_logit(subset_matrix(X,K,g,gg) * delta)) ) /
                sum(subset_vector(weights,g,gg));
  }
  ordered[J] opsiM; // coefficients of controls without treatment
  opsiM[1] = 0;
  for (k in 2:J){
    opsiM[k] = opsiM[k-1] + exp(psiM[k-1]);
  }
  vector<lower=0, upper=1>[N] pz; // 
  for (i in 1:N){
    pz[i] = inv_logit(X[i,] * delta);
  }
}

model{

	delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
	psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment

	psiM  ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)

	for (gg in 1:G){
	  pas[gg] ~ beta(a[gg], b[gg]);
	}

	for (i in 1:N){
	  if (treat[i] == 1 && Y[i] == 0) {
	    target += log(1 - pz[i]) + ordered_logistic_lupmf( 1 | X[i,] * psi0, opsiM);
	  } else if (treat[i] == 1 && Y[i] == J + 1){
	    target += log(    pz[i]) + ordered_logistic_lupmf(J+1| X[i,] * psi0, opsiM);
	  } else if (treat[i] == 1){
	    target += log(pz[i]*exp(ordered_logistic_lupmf(Y[i]   | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lupmf(Y[i]+1| X[i,] * psi0, opsiM)) );
	  } else { // no treatment
	    target += log(pz[i]*exp(ordered_logistic_lupmf(Y[i]+1 | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lupmf(Y[i]+1 | X[i,] * psi0, opsiM)) );
	  }
	}

}

generated quantities {
    vector[N] log_lik;
    for (i in 1:N){
      if (treat[i] == 1 && Y[i] == 0) {
        log_lik[i] = log(1 - pz[i]) + ordered_logistic_lpmf( 1 | X[i,] * psi0, opsiM);
      } else if (treat[i] == 1 && Y[i] == J + 1){
        log_lik[i] = log(    pz[i]) + ordered_logistic_lpmf(J+1| X[i,] * psi0, opsiM);
      } else if (treat[i] == 1){
        log_lik[i] = log(pz[i]*exp(ordered_logistic_lpmf(Y[i]   | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lpmf(Y[i]+1| X[i,] * psi0, opsiM)) );
      } else { // no treatment
        log_lik[i] = log(pz[i]*exp(ordered_logistic_lpmf(Y[i]+1 | X[i,] * psi0, opsiM)) + (1 - pz[i]) * exp(ordered_logistic_lpmf(Y[i]+1 | X[i,] * psi0, opsiM)) );
      }
    }
} 
