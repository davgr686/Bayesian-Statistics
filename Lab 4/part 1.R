library("rstan") # observe startup messages
set.seed(1234)


#a
u = 10 
sigma = 2
t = 200
phi = 1


ARSim = function (u, sigma, t, phi){
  x = u
  x_vec = c()
  for( j in 1:t ){
    err = rnorm(1,0,sqrt(sigma))
    x_t = u + phi*(x - u)+ err
    x = x_t
    x_vec[j] = x_t
  }
  return(x_vec)
}
plot(ARSim(u, sigma, t, -1), type="l")
plot(ARSim(u, sigma, t, 0), type="l")
plot(ARSim(u, sigma, t, 0.5), type="l")
plot(ARSim(u, sigma, t, 1), type="l")


#b

X <- ARSim(u, sigma, t, 0.3)
Y <- ARSim(u, sigma, t, 0.95)
plot(X, type="l")
plot(Y, type="l")



X_list <- list(N=t,y=X)
Y_list <- list(N=t,y=Y)

stan_model <- 'data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=-1,upper=1> phi;
}
model {
  for (n in 2:N)
    y[n] ~ normal(mu + phi * (y[n-1] - mu), sqrt(sigma));
}
'
fit_X <- stan(model_code = stan_model, data = X_list, iter = 5000, warmup = 1000)
fit_Y <- stan(model_code = stan_model, data = Y_list, iter = 5000, warmup = 1000)

X_draws = extract(fit_X)
Y_draws = extract(fit_Y)



summary(X_draws$mu)
summary(X_draws$sigma)
summary(X_draws$phi)

summary(Y_draws$mu)
summary(Y_draws$sigma)
summary(Y_draws$phi)


mu_cred_interval_x <- quantile(X_draws$mu, c(0.025, 0.975))
sigma_cred_interval_x <- quantile(X_draws$sigma, c(0.025, 0.975))
phi_cred_interval_x <- quantile(X_draws$phi, c(0.025, 0.975))

mu_cred_interval_y <- quantile(Y_draws$mu, c(0.025, 0.975))
sigma_cred_interval_y <- quantile(Y_draws$sigma, c(0.025, 0.975))
phi_cred_interval_y <- quantile(Y_draws$phi, c(0.025, 0.975))



traceplot(fit_X)
traceplot(fit_Y)

plot(X_draws$mu, X_draws$phi)
plot(Y_draws$mu, Y_draws$phi)

#fit_Y <- stan(model_code = stan_model, data = Y)


