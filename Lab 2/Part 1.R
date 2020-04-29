# a

library(mvtnorm)
set.seed(1234)
setwd('C:/Users/daavi/Documents/Bayesian Learning - 732A91/Lab 2')

data <- read.delim("TempLinkoping.txt", sep='\t')
data$temp <- ifelse(!is.na(data$X), data$X, data$temp)
data <- subset(data, select = -c(2))


mu_0 <- c(-10, 130, -130)
omega_0 <- 1 * diag(3)
var_0 <- 1
v_0 <- 4

predict_sim <- function(v_0, var_0, omega_0, x_vals) {
  
  var <- (v_0*var_0)/rchisq(1,v_0)
  betas <- rmvnorm(1, mean = mu_0, sigma = var * solve(omega_0))
  
  preds <- array(1:length(x_vals))
  for (i in 1:length(x_vals)) {
    preds[i] <- betas[,1] + betas[,2] * x_vals[i] + betas[,3] * x_vals[i]**2
  }
  return (preds)
}

plot(data$time, data$temp)
for (i in 1:100) {
  lines(data$time, predict_sim(v_0, var_0, omega_0, data$time), type="l", col="red")
}

# b
X <- as.matrix(cbind(1, data$time, data$time**2))
y <- as.matrix(data$temp)

beta_hat <- solve(t(X)%*%X)%*%(t(X)%*%y)
mu_n <- solve(t(X)%*%X + omega_0) %*% (t(X)%*%X%*%beta_hat + omega_0%*%mu_0)
omega_n <- t(X)%*%X + omega_0
v_n <- v_0 + length(data$time)
var_n <- (v_0 * var_0 + (t(y)%*%y + t(mu_0)%*%omega_0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n))/v_n

posterior_draw <- function(v_n, var_n, mu_n, omega_n) {
  
  var <- (v_n*var_n)/rchisq(1,v_n)
  betas <- rmvnorm(1, mean = mu_n, sigma = var[1] * solve(omega_n))
  list(betas=betas, var = var)
}

n_sim <- 10000

posterior_betas <- matrix(0,n_sim,3)
posterior_vars <- c()


median_regression <- matrix(0,n_sim,length(data$temp))

for (i in 1:n_sim) {
  posterior_sim <- posterior_draw(v_n, var_n, mu_n, omega_n )
  posterior_betas[i,] <- posterior_sim$betas
  posterior_vars[i] <- posterior_sim$var
  
  for (x in 1:length(data$temp)) {
    median_regression[i,x] <- posterior_betas[i,1] + posterior_betas[i,2] * data$time[x] + posterior_betas[i,3] * data$time[x]**2
  }
}

medians <- apply(median_regression, 2, median)

cred_interval_lower <- c()
cred_interval_upper <- c()


for (u in 1:length(data$temp)) {
  cred_interval <- quantile(median_regression[,u],probs=c(0.025,0.975))
  cred_interval_lower[u] <- cred_interval[1]
  cred_interval_upper[u] <- cred_interval[2]
}

plot(data$time, data$temp, col="blue", xlab="Time", ylab="Temperature")
lines(data$time, medians, col="black")
lines(data$time, cred_interval_lower, col="red")
lines(data$time, cred_interval_upper, col="green")

par(mfrow=c(2,2))
hist(posterior_betas[,1], 100, main="Posterior Distibution of Beta 0", xlab="beta 0")
hist(posterior_betas[,2], 100, main="Posterior Distibution of Beta 1", xlab="beta 1")
hist(posterior_betas[,3], 100, main="Posterior Distibution of Beta 2", xlab="beta 2")
hist(posterior_vars, 100, main="Posterior Distibution of sigma^2", xlab="sigma^2")
par(mfrow=c(1,1))
# c

x_hat <- -posterior_betas[,2] / (2*posterior_betas[,3])
hist(x_hat, 100)

# d






