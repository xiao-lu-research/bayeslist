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

	// normal link for continuous dv
	real fd_l(int y, real p){
		return(binomial_lpmf(y|1,p));
	}


}

data {
	int N;
	int J; //number of non sensitive item
	int<lower = 0> Y[N]; // number of affirmative answers
	int K;
	matrix[N,K] X;
	int treat[N];
	int outcome[N];
	real mu_delta[K];
	real<lower = 0> sigma_delta[K];
	real mu_psi0[K];
	real<lower = 0> sigma_psi0[K];
	real mu_psi2[K];
	real<lower = 0> sigma_psi2[K];
	real mu_rho0;
	real<lower = 0> sigma_rho0;
	real mu_phi;
	real<lower = 0> sigma_phi;
	real mu_gamma0;
	real<lower = 0> sigma_gamma0;
}

transformed data{

}

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	vector[K] psi2; // coefficients of controls on outcome
	vector[K] delta; // coefficients of controls on Z
	real gamma0; // coefficient of the sensitive item on outcome
	real phi; // coefficients of latent number of affirmative answers to controls on outcome
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
}



model{

	for (i in 1:K){
		delta[i]  ~ normal(mu_delta[i],sigma_delta[i]); // priors subject to modification
		psi0[i]  ~ normal(mu_psi0[i],sigma_psi0[i]); // priors subject to modification
		psi2[i]  ~ normal(mu_psi2[i],sigma_psi2[i]); // priors subject to modification
	}
	
	rho0 ~ normal(mu_rho0, sigma_rho0);
	phi ~ normal(mu_phi, sigma_phi);
	gamma0 ~ normal(mu_gamma0,sigma_gamma0);


	for (i in 1:N){
		if (treat[i] == 1 && Y[i] == 0) {
			target += fd_l(outcome[i], inv_logit(X[i,] * psi2 + Y[i] * phi + 0 * gamma0)) +  log(1 - exp(gd(inv_logit(X[i,] * delta)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && Y[i] == J + 1){
			target += fd_l(outcome[i], inv_logit(X[i,] * psi2 + (Y[i] - 1) * phi + 1 * gamma0)) + gd(inv_logit(X[i,] * delta)) + h1d(J,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1){
			target += log(exp(fd_l(outcome[i], inv_logit(X[i,] * psi2 + (Y[i] - 1) * phi + 1 * gamma0))) * exp(gd(inv_logit(X[i,] * delta))) * exp(h1d(Y[i] - 1,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + exp(fd_l(outcome[i], inv_logit(X[i,] * psi2 + Y[i] * phi + 0 * gamma0)))  * (1 - exp(gd(inv_logit(X[i,] * delta)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else {
			target += log(exp(fd_l(outcome[i], inv_logit(X[i,] * psi2 + Y[i] * phi + 1 * gamma0))) * exp(gd(inv_logit(X[i,] * delta))) * exp(h1d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + exp(fd_l(outcome[i], inv_logit(X[i,] * psi2 + Y[i] * phi + 0 * gamma0))) * (1 - exp(gd(inv_logit(X[i,] * delta)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		}
	}




}