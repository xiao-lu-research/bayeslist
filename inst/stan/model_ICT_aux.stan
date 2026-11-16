functions{
	

	real gd(real p){
		return(binomial_lpmf(1|1,p));
	}

	real h0d(int y, int J, real alpha, real beta){
		return(beta_binomial_lpmf(y|J,alpha,beta));
	}

	real h1d(int y, int J, real alpha, real beta){
		return(beta_binomial_lpmf(y|J,alpha,beta));
	}

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
        real<lower = 0> sigma_bias; // prior bias standard deviation
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
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
	// real<lower = 0, upper = 1> rho1; //variance parameter with treatment
        real bias; // bias
}

transformed parameters {
  vector<lower=0, upper=1>[G] pas;
  for (gg in 1:G){
// returns X[g==gg,1:K]*delta
    pas[gg] = dot_product(subset_vector(weights,g,gg),
                (inv_logit(subset_matrix(X,K,g,gg) * delta)) ) /
                sum(subset_vector(weights,g,gg));
  }
}

model{

  	bias ~ normal(0,sigma_bias);
	delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
	psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
	// psi1  ~ normal(0,10); // priors subject to modification

	for (gg in 1:G){
		pas[gg] ~ beta(a[gg], b[gg]);
	}

	for (i in 1:N){
		if (treat[i] == 1 && Y[i] == 0) {
			target += log(1 - exp(gd(inv_logit(X[i,] * delta + bias)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && Y[i] == J + 1){
			target += gd(inv_logit(X[i,] * delta + bias)) + h1d(J,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1){
			target += log(exp(gd(inv_logit(X[i,] * delta + bias))) * exp(h1d(Y[i] - 1,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta + bias)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else {
			target += log(exp(gd(inv_logit(X[i,] * delta + bias))) * exp(h1d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta + bias)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		}
		//target += normal_lpdf(logit(h[g[i]])| X[i,] * delta + bias, 1);
	}

}

generated quantities {
    vector[N] log_lik;
    for (i in 1:N){
	if (treat[i] == 1 && Y[i] == 0) {
		log_lik[i] = log(1 - exp(gd(inv_logit(X[i,] * delta + bias)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
	} else if (treat[i] == 1 && Y[i] == J + 1){
		log_lik[i] = gd(inv_logit(X[i,] * delta + bias)) + h1d(J,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
	} else if (treat[i] == 1){
		log_lik[i] = log(exp(gd(inv_logit(X[i,] * delta + bias))) * exp(h1d(Y[i] - 1,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta + bias)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
	} else {
		log_lik[i] = log(exp(gd(inv_logit(X[i,] * delta + bias))) * exp(h1d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta + bias)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
	}
    }
} 
