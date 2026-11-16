functions{
	

	real h0d(int y, int J, real alpha, real beta){
//print("h0d y ",y," J ",J," alpha ",alpha," beta ",beta);
		return(beta_binomial_lpmf(y|J,alpha,beta));
	}

	real h1d(int y, int J, real alpha, real beta){
//print("h1d y ",y," J ",J," alpha ",alpha," beta ",beta);
		return(beta_binomial_lpmf(y|J,alpha,beta));
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
        array[K] real<lower = 0> sigma_psi0; // prior psi0 standard deviation
        real<lower = 0, upper = 1> mu_rho0; //variance parameter without treatment
        real<lower = 0> sigma_rho0; //prior SD of variance parameter without treatment
}

// transformed data{
//   real<lower=0> a;
//   real<lower=0> b;
//   a = mu_rho0 * ((mu_rho0 * (1 - mu_rho0) / (sigma_rho0*sigma_rho0)) - 1);
//   b = (1 - mu_rho0) * ((mu_rho0 * (1 - mu_rho0) / (sigma_rho0*sigma_rho0)) - 1);
// }

parameters {
	vector[K] psi0; // coefficients of controls without treatment
	vector[K] delta; // coefficients of controls on Z
        real lrho0; // logit variance parameter without treatment
}

transformed parameters {
  real<lower=0, upper=1> pas;
  pas = dot_product(weights, (inv_logit(X * delta)));

  real<lower = 0, upper = 1> rho0; //variance parameter without treatment
  rho0 = inv_logit(lrho0);

  vector<lower=0, upper=1>[N] pz; // 
  vector<lower=0>[N] il; // 
  vector<lower=0>[N] ilc; // 
  for (i in 1:N){
    pz[i] = inv_logit(X[i,] * delta);
    il[i] = inv_logit(X[i,] * psi0)/exp(lrho0);
    ilc[i] = (1 - inv_logit(X[i,] * psi0))/exp(lrho0);
  }
}

model{
  //    Note that the next line is the standard prior
    	delta ~ normal(0,sigma_delta); // prior for the coefficients of controls on Z
	psi0  ~ normal(0,sigma_psi0); // prior for the coefficients of controls (with or without treatment)

    //  rho0 ~ beta(a,b);
    //  rho0 ~ uniform(0.0,1.0);
        lrho0 ~ normal(logit(mu_rho0), sigma_rho0);
    //  lrho0 ~ normal(0, 2);

    //  Non-degeneracy prior on PAS
    //  target += log((N+1)/(N-1)) + log(1-(pow(pas, N) + pow(1-pas, N))); // ND prior on pas
    //  remember the log absolute determinant of the Jacobian
    //	target += log(dot_product(weights, inv_logit(X * delta) .* (1-inv_logit(X * delta)) ));
    //  for (j in 1:K){
    //    target += log(dot_product(weights .* col(X,j), (inv_logit(X * delta) .* (1-inv_logit(X * delta))) ));
    //  }

  	for (i in 1:N){
//print("model i ",i," J ",J," il[i] ",il[i]," pz[i] ",pz[i]);
		if (treat[i] == 1 && Y[i] == 0) {
			target += log(1 - pz[i]) + h0d(0,J,il[i],ilc[i]);
		} else if (treat[i] == 1 && Y[i] == J + 1){
			target += log(pz[i]) + h1d(J,J,il[i],ilc[i]);
		} else if (treat[i] == 1){
			target += log(pz[i] * exp(h1d(Y[i] - 1,J,il[i],ilc[i])) + (1 - pz[i]) * exp(h0d(Y[i],J,il[i], ilc[i])));
		} else { // no treatment
			target += log(pz[i] * exp(h1d(Y[i],J,il[i], ilc[i])) + (1 - pz[i]) * exp(h0d(Y[i],J,il[i],ilc[i])) );
		}
  	}

}

generated quantities {
    vector[N] log_lik;  
    for (i in 1:N){
//print(" generated i ",i,"J ",J,"il[i] ",il[i],"pz[i]",pz[i]);
 	if (treat[i] == 1 && Y[i] == 0) {
 		log_lik[i]= log(1 - pz[i]) + h0d(0,J,il[i], ilc[i]);
 	} else if (treat[i] == 1 && Y[i] == J + 1){
 		log_lik[i]= log(pz[i]) + h1d(J,J,il[i], ilc[i]);
 	} else if (treat[i] == 1){
 		log_lik[i]= log(pz[i] * exp(h1d(Y[i] - 1,J,il[i], ilc[i])) + (1 - pz[i]) * exp(h0d(Y[i],J,il[i], ilc[i])));
 	} else { // no treatment
 		log_lik[i]= log(pz[i] * exp(h1d(Y[i],J,il[i], ilc[i])) + (1 - pz[i]) * exp(h0d(Y[i],J,il[i], ilc[i])) );
}
     }
}
