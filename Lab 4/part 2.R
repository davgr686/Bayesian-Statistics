library("rstan") # observe startup messages
set.seed(1234)


data <- read.table("https://raw.githubusercontent.com/STIMALiU/BayesLearnCourse/master/Labs/campy.dat", header = T)

stan_model <- 'data {
  int<lower=0> N;
  int c[N];
}
parameters {
  real x[N];
  real mu;
  real<lower=0> sigma;
  real<lower=-1,upper=1> phi;
}

model {
  sigma ~ scaled_inv_chi_square(1,1);
  
  for (n in 2:N) {
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sqrt(sigma));
    c[n] ~ poisson(exp(x[n]));
  }
}
'

X_list <- list(N=dim(data)[1], c=data$c)

fit <- stan(model_code = stan_model, data = X_list, iter = 5000, warmup = 1000)
draws = extract(fit)
traceplot(fit)

plot(data$c)

x_means <- colMeans(draws$x)
lines(exp(x_means), col="red", lwd = 2)

x_cred_interval_25 <- summary(fit)$summary[,"2.5%"][1:140]
x_cred_interval_97 <- summary(fit)$summary[,"97.5%"][1:140]

lines(exp(x_cred_interval_25), col="green", lwd = 1)
lines(exp(x_cred_interval_97), col="green", lwd = 1)

# d

stan_model2 <- 'data {
  int<lower=0> N;
  int c[N];
}
parameters {
  real x[N];
  real mu;
  real<lower=0> sigma;
  real<lower=-1,upper=1> phi;
}

model {
  sigma ~ scaled_inv_chi_square(N,0.1);
  
  for (n in 2:N) {
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sqrt(sigma));
    c[n] ~ poisson(exp(x[n]));
  }
}
'

fit <- stan(model_code = stan_model2, data = X_list, iter = 5000, warmup = 1000)
draws = extract(fit)
traceplot(fit)

plot(data$c)

x_means <- colMeans(draws$x)
lines(exp(x_means), col="red", lwd = 2)

x_cred_interval_25 <- summary(fit)$summary[,"2.5%"][1:140]
x_cred_interval_97 <- summary(fit)$summary[,"97.5%"][1:140]

lines(exp(x_cred_interval_25), col="green", lwd = 1)
lines(exp(x_cred_interval_97), col="green", lwd = 1)
