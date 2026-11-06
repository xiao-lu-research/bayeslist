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
	vector[N] weights; // survey weights
}

parameters {
	vector[K] psi0; // coefficients of first list controls without treatment 1
	vector[K] delta; // coefficients of first list controls on Z
	vector[K] psi02; // coefficients of controls without treatment 2
	vector[K] delta2; // coefficients of second list controls on Z
	real<lower = 0, upper = 1> rho0; //variance parameter without treatment
	real<lower = 0, upper = 1> rho02; //variance parameter without treatment
}


transformed parameters {
  real<lower=0, upper=1> pas; // proportion affirmative on sensitive items
  pas = dot_product(weights, 0.5*(inv_logit(X * delta))+0.5*(inv_logit(X2 * delta2)));
}

model{

        delta ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi0  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
        delta2 ~ normal(mu_delta,sigma_delta); // prior for the coefficients of controls on Z
        psi02  ~ normal(mu_psi0,sigma_psi0); // prior for the coefficients of controls without treatment
	
//      rho0 ~ normal(mu_rho0, sigma_rho0);
  	rho0 ~ uniform(0, 1);
//	rho02 ~ uniform(0, 1);

	for (i in 1:N){
// First those that received ttt 1 (and not ttt 2)
		if (treat[i] == 1 && Y1[i] == 0) { 
			target += log(1 - exp(gd(inv_logit(X[i,] * delta)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && Y1[i] == (J + 1)){
			target += gd(inv_logit(X[i,] * delta)) + h1d(J,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1){ // 0 < Y1[i] < J+1
			target += log(exp(gd(inv_logit((X[i,] * delta)))) * exp(h1d(Y1[i] - 1,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit((X[i,] * delta))))) * exp(h0d(Y1[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else {
			target += log(exp(gd(inv_logit((X[i,] * delta)))) * exp(h1d(Y1[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit((X[i,] * delta))))) * exp(h0d(Y1[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		}

// Now those that received ttt 2 (and not ttt 1)

		if (treat2[i] == 1 && Y2[i] == 0) {
			target += log(1 - exp(gd(inv_logit((X2[i,] * delta2))))) + h0d(0,J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02);
		} else if (treat2[i] == 1 && Y2[i] == (J + 1)){
			target += gd(inv_logit((X2[i,] * delta2))) + h1d(J,J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02);
		} else if (treat2[i] == 1){
			target += log(exp(gd(inv_logit((X2[i,] * delta2)))) * exp(h1d(Y2[i] - 1,J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)) + (1 - exp(gd(inv_logit((X2[i,] * delta2))))) * exp(h0d(Y2[i],J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)));
		} else {
			target += log(exp(gd(inv_logit((X2[i,] * delta2)))) * exp(h1d(Y2[i],J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)) + (1 - exp(gd(inv_logit((X2[i,] * delta2))))) * exp(h0d(Y2[i],J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)));
		}
	}

}

generated quantities {
    vector[N] log_lik; 
    for (i in 1:N){
// First those that received ttt 1 (and not ttt 2)
		if (treat[i] == 1 && Y1[i] == 0) { 
			log_lik[i] = log(1 - exp(gd(inv_logit(X[i,] * delta)))) + h0d(0,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1 && Y1[i] == (J + 1)){
			log_lik[i] = gd(inv_logit(X[i,] * delta)) + h1d(J,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0);
		} else if (treat[i] == 1){ // 0 < Y1[i] < J+1
			log_lik[i] = log(exp(gd(inv_logit((X[i,] * delta)))) * exp(h1d(Y1[i] - 1,J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit((X[i,] * delta))))) * exp(h0d(Y1[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		} else {
			log_lik[i] = log(exp(gd(inv_logit((X[i,] * delta)))) * exp(h1d(Y1[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)) + (1 - exp(gd(inv_logit((X[i,] * delta))))) * exp(h0d(Y1[i],J,inv_logit(X[i,]*psi0)*(1 - rho0)/rho0, (1 - inv_logit(X[i,]*psi0))*(1 - rho0)/rho0)));
		}

// Now those that received ttt 2 (and not ttt 1)

		if (treat2[i] == 1 && Y2[i] == 0) {
			log_lik[i] += log(1 - exp(gd(inv_logit((X2[i,] * delta2))))) + h0d(0,J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02);
		} else if (treat2[i] == 1 && Y2[i] == (J + 1)){
			log_lik[i] += gd(inv_logit((X2[i,] * delta2))) + h1d(J,J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02);
		} else if (treat2[i] == 1){
			log_lik[i] += log(exp(gd(inv_logit((X2[i,] * delta2)))) * exp(h1d(Y2[i] - 1,J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)) + (1 - exp(gd(inv_logit((X2[i,] * delta2))))) * exp(h0d(Y2[i],J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)));
		} else {
			log_lik[i] += log(exp(gd(inv_logit((X2[i,] * delta2)))) * exp(h1d(Y2[i],J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)) + (1 - exp(gd(inv_logit((X2[i,] * delta2))))) * exp(h0d(Y2[i],J,inv_logit(X2[i,]*psi02)*(1 - rho02)/rho02, (1 - inv_logit(X2[i,]*psi02))*(1 - rho02)/rho02)));
		}
	}

}
