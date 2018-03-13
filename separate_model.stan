//analysis of longitudinal data
data{
  // qol data
  int<lower=0> N; //number of total observations
  int<lower=0> I; //number of subjects
  vector[N] KCCQ; //quality of life for each patient at each time
  vector[N] Timepoint; //timepoints
  vector[I] Treatment; //treatment indicators
  int<lower = 1, upper = I> Subjects[N]; //indicator for subjects
  
  // survival data
  vector<lower=0>[I] Y; //survival time
  int Cen[I]; //censoring indicator  
}


parameters{
  //qol parameters
  real mu_theta_0; //intercept for control group (baseline)
  real mu_theta_1; //slope for control group (baseline)
  //real mu_theta_1_t; //slope for treatment group
  real<lower = 0> sigma_theta_0; //standard deviation for intercepts
  real<lower = 0> sigma_theta_1; //standard deviation for slopes
  vector[I] theta_0_raw; //patient specific intercepts
  vector[I] theta_1_raw; //patient specific slopes
  real gamma_0; //population level treatment effect (intercept)
  real gamma_1; //population level treatment slope
  real<lower = 0> sigma_kccq; //standard deviation for kccq score
  
  //survival parameters
  real<lower=0> alpha; //shape parameter
  real beta_0; //intercept (baseline treatment = BMS)
  real beta_t; //treatment coefficient
}

transformed parameters{
  vector[N] linear_predictor; //fitted values
  vector[I] theta_0;
  vector[I] theta_1;
  theta_0 = theta_0_raw * sigma_theta_0; //theta_0 is scaled by sigma_theta_0 (its variance)
  theta_1 = theta_0_raw * sigma_theta_1; //theta_1 is scaled by sigma_theta_1 (its variance)
  
  for(n in 1:N){
    linear_predictor[n] = mu_theta_0 + theta_0[Subjects[n]] + gamma_0 * Treatment[Subjects[n]] + gamma_1 * Treatment[Subjects[n]] * Timepoint[n] +  mu_theta_1 * Timepoint[n] + theta_1[Subjects[n]] * Timepoint[n]; 
  }
}
  

model{
  //longitudinal model priors
  theta_0_raw ~ normal(0, 1);
  theta_1_raw ~ normal(0, 1); //sigma_theta_1 tends to collapse to zero
  sigma_theta_0 ~ student_t(3, 0, 1);
  sigma_theta_1 ~ student_t(3, 0, 1);
  sigma_kccq ~ student_t(3, 0, 1);
  mu_theta_0 ~ normal(0, 1);
  mu_theta_1 ~ normal(0, 1);
  gamma_0 ~ normal(0, 1);
  gamma_1 ~ normal(0, 1);
  
  //survival model priors
  alpha ~ normal(0, 10);
  beta_0 ~ normal(0, 10);
  beta_t ~ normal(0, 10);

  //longitudinal model
  KCCQ ~ normal(linear_predictor, sigma_kccq);
  //survival model
  for(i in 1:I){
    if(Cen[i] == 0){
      Y[i] ~ weibull(alpha, exp(-(beta_0 + Treatment[i] * beta_t) / alpha));
  } else {
      target += weibull_lccdf(Y[i] | alpha, exp(-(beta_0 + Treatment[i] * beta_t) / alpha));
   }
 }
}

generated quantities{
  vector[N] bounded_predictor; //variable on the original KCCQOS scale
  bounded_predictor = Phi(linear_predictor);
}
  
