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
	real<lower = 0> a[K];
	real<lower = 0> b[K];
}

transformed data{

}

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	// vector[K] psi1; // coefficients of controls with treatment
	vector[K] psi2; // coefficients of controls on outcome
	vector<lower = 0, upper = 1>[K] delta00; // auxiliary
	real gamma0; // coefficient of the sensitive item on outcome
	real phi; // coefficients of latent number of affirmative answers to controls on outcome
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
	// real<lower = 0, upper = 1> rho1; //variance parameter with treatment
	//real<lower = 0> sigma; // scale parameter for normal outcome
}

transformed parameters {
	vector[K] delta; // coefficients of controls of sensitive item 
	for (i in 1:K){
		delta[i] = log(1/(1/delta00[i] - 1));
	}
}


model{

	psi0  ~ normal(0,10); // priors subject to modification
	psi2  ~ normal(0,10); // priors subject to modification
	gamma0 ~ normal(0,10);
	for (i in 1:K){
		delta00[i]  ~ beta(a[i]+1,b[i]+1); // priors subject to modification
	}


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