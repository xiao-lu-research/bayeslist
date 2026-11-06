functions{

	int group_size(array[] int ref, int value) {
    	  int count;
    	  count = 0;
    	  for (ii in 1:size(ref))
      	  if (ref[ii]==value)
        	  count = count + 1;
      	  return count;
  	}

// 	returns the subset of y[ref==value,1:k]
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
        real<lower = 0> sigma_bias; // prior bias standard deviation
        array[J-1] real<lower = 0> sigma_psiM; // prior psiM standard deviation
	vector[N] weights; // survey weights
	int G; // number of subgroups
	array[G] real h; // auxiliary information for subgroups
	array[N] int g; // subgroup indicators
	array[G] real sigma_h; // subgroup standard deviations
}

//transformed data {
//	real logith[G]; // auxiliary information for subgroups (on the logit scale)
//        for (gg in 1:G){
//          logith[gg] = logit(h[gg]);
//        }
//}

transformed data{
  vector<lower=0>[G] a;
  vector<lower=0>[G] b;
  for (gg in 1:G){
    a[gg] = h[gg] * ((h[gg] * (1 - h[gg]) / (sigma_h[gg]*sigma_h[gg])) - 1);
    b[gg] = (1 - h[gg]) * ((h[gg] * (1 - h[gg]) / (sigma_h[gg]*sigma_h[gg])) - 1);
  }
}

parameters {
	vector[K] psi0; // coefficients of first list controls without treatment 1
	vector[K] delta; // coefficients of first list controls on Z
	vector[K] psi02; // coefficients of controls without treatment 2
	vector[K] delta2; // coefficients of second list controls on Z
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
	real<lower = 0, upper = 1> rho02; //variance parameter without treatment
	ordered[J-1] psiM; // coefficients of controls without treatment
	ordered[J-1] psiM2; // coefficients of controls without treatment
}


//  real<lower=0, upper=1> pas; // proportion affirmative on sensitive items
//  pas = dot_product(weights, 0.5*(inv_logit(X * delta))+0.5*(inv_logit(X2 * delta2)));

transformed parameters {
  vector<lower=0, upper=1>[G] pas; // proportion affirmative on sensitive items
  for (gg in 1:G){
    pas[gg] = dot_product(subset_vector(weights,g,gg),
     0.5*(inv_logit(subset_matrix(X,      K,g,gg) * delta )) +
     0.5*(inv_logit(subset_matrix(X2,K,g,gg) * delta2)) ) /
     sum(subset_vector(weights,g,gg));
  }
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
    pz[i] = inv_logit(X[i,] * delta);
   pz2[i] = inv_logit(X2[i,] * delta2);
  }
}

model{

        delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
        delta2 ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi02  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
	
	psiM  ~ normal(0,sigma_psiM); // prior for the coefficients of controls (with or without treatment)

	for (gg in 1:G){
  	  pas[gg] ~ beta(a[gg], b[gg]);
	}

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
