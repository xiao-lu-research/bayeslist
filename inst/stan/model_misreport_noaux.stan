functions{
	

	real gd(real p){
		return(binomial_lpmf(1|1,p));
	}

	real jd(real p){
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
	int direct[N];
}

transformed data{

}

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	//vector[K] psi1; // coefficients of controls with treatment
	vector[K] delta; // coefficients of controls on Z
	vector[K] gamma0; // coefficients of controls on U
	real treat_e; // effect of treatment on U
	real U_e; // effect of U on Y*
	real Z_e; // effect of Z on Y*
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
	//real<lower = 0, upper = 1> rho1; //variance parameter with treatment
}



model{

	gamma0 ~ normal(0,10); // priors subject to modification
	delta ~ normal(0,10); // priors subject to modification
	psi0  ~ normal(0,10); // priors subject to modification
	//psi1  ~ normal(0,10); // priors subject to modification
	treat_e  ~ normal(0,10); // priors subject to modification
	U_e  ~ normal(0,10); // priors subject to modification
	Z_e  ~ normal(0,10); // priors subject to modification

	for (i in 1:N){

		if (treat[i] == 0 && direct[i] == 0){
			target += log(exp(jd(inv_logit(X[i,] * gamma0))) * exp(gd(inv_logit(X[i,] * delta))) * exp(h1d(Y[i],J,inv_logit(X[i,]*psi0 + Z_e + U_e)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0 + Z_e + U_e))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else if (treat[i] == 0 && direct[i] == 1) {
			target += log(1 - exp(jd(inv_logit(X[i,] * gamma0)))) + gd(inv_logit(X[i,] * delta)) + h1d(Y[i],J,inv_logit(X[i,]*psi0 + Z_e)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0 + Z_e))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && direct[i] == 0 && Y[i] == 0) {
			target += log(1 - exp(gd(inv_logit(X[i,] * delta)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && direct[i] == 0 && Y[i] == J + 1) {
			target += jd(inv_logit(X[i,] * gamma0 + treat_e)) + gd(inv_logit(X[i,] * delta)) + h1d(J,J,inv_logit(X[i,]*psi0 + Z_e + U_e)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0 + Z_e + U_e))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && direct[i] == 0) {
			target += log(exp(jd(inv_logit(X[i,] * gamma0 + treat_e ))) * exp(gd(inv_logit(X[i,] * delta))) * exp(h1d(Y[i] - 1,J,inv_logit(X[i,]*psi0 + Z_e + U_e)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0 + Z_e + U_e))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit(X[i,] * delta)))) * exp(h0d(Y[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else if (treat[i] == 1 && direct[i] == 1 && Y[i] > 0){
			target += log(1 - exp(jd(inv_logit(X[i,] * gamma0 + treat_e)))) + gd(inv_logit(X[i,] * delta)) + h1d(Y[i] - 1,J,inv_logit(X[i,]*psi0 + Z_e)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0 + Z_e))*(1 - rho0)/rho0);
		} else {
			target += 1;
		}
			
			
	}




}