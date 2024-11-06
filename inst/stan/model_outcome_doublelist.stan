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


}

data {
	int N;
	int J; //number of non sensitive item
	int<lower = 0> Y[N]; // number of affirmative answers
	int K;
	matrix[N,K] X;
	int treat[N];
	real mu_psi0[K];
	real<lower = 0> sigma_psi0[K];
	real mu_delta[K];
	real<lower = 0> sigma_delta[K];
	real mu_rho0;
	real<lower = 0> sigma_rho0;
}


parameters {
	vector[K] psi0; // coefficients of controls without treatment
	vector[K] delta; // coefficients of controls on Z
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
}



model{



	

	for (i in 1:K){
		psi0[i]  ~ normal(mu_psi0[i],sigma_psi0[i]); // priors subject to modification
		delta[i] ~ normal(mu_delta[i],sigma_delta[i]); // priors subject to modification
	}
	
	rho0 ~ normal(mu_rho0, sigma_rho0);

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
	}




}