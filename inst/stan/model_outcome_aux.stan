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

	int group_size(int[] ref, int value) {
    int count;
    count = 0;
    for (ii in 1:size(ref))
      if (ref[ii]==value)
        count = count + 1;
      return count;
  	}

  matrix subset_matrix(matrix y, int k, int[] ref, int value) {
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

}

data {
	int N;
	int J; //number of non sensitive item
	int<lower = 0> Y[N]; // number of affirmative answers
	int K;
	matrix[N,K] X;
	int treat[N];
  // int N1; // length of the indicator
	int G; // number of subgroups
	real h[G]; // auxiliary information for subgroups
	int g[N]; // subgroup indicators
}

transformed data{

}

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	// vector[K] psi1; // coefficients of controls with treatment
	vector[K] delta; // coefficients of controls on Z
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
	// real<lower = 0, upper = 1> rho1; //variance parameter with treatment
	real<lower = 0> sigma; //sd for the moment conditions of the auxiliary information
}



model{


	sigma ~ cauchy(0,2.5);
	delta ~ normal(0,10); // priors subject to modification
	psi0  ~ normal(0,10); // priors subject to modification
	// psi1  ~ normal(0,10); // priors subject to modification

	for (gg in 1:G){
		logit(h[gg]) ~ normal(mean(subset_matrix(X,K,g,gg) * delta), sigma); 
	}

	for (i in 1:N){
		if (treat[i] == 1 && Y[i] == 0) {
			target += log(1 - exp(gd(inv_logit(X[i,] * delta)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && Y[i] == J + 1){
			target += gd(inv_logit(X[i,] * delta)) + h1d(J,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1){
			target += log(exp(gd(inv_logit(X[i,] * delta))) * exp(h1d(Y[i] - 1,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else {
			target += log(exp(gd(inv_logit(X[i,] * delta))) * exp(h1d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		}
		//target += normal_lpdf(logit(h[g[i]])| X[i,] * delta, 1);
	}




}