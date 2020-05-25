set.seed(1234)
data <- read.table("eBayNumberOfBidderData.dat", header = T)

# a

model <- glm(nBids ~ . - Const, family="poisson", data=data)
barplot(model$coefficients, ylab = "Beta", xlab="Covariate", main="Covariates with their corresponding beta values")
summary(model)

# b

library(mvtnorm)
y <- data$nBids
X <- as.matrix(data[,-c(1)])

Sigma <- 100 * solve(t(X)%*%X)
mu <- matrix(0, nrow=1, ncol=9)


LogPoisson <- function(betas,y,X,mu,Sigma){
  linPred <- betas %*% t(X)
  
  logLik <- sum(linPred * y - exp(linPred))
  if (abs(logLik) == Inf) logLik = -20000
  
  logPrior <- dmvnorm(betas, matrix(0, length(betas), 1), Sigma, log=TRUE) ## zellners prior
  return(logLik + logPrior)
}
initVal <- as.vector(rep(0,dim(X)[2]));

OptimResults<-optim(initVal,LogPoisson,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

beta_tilde <- OptimResults$par
j_inv <- -solve(OptimResults$hessian)
beta_draws <- rmvnorm(10000,beta_tilde, j_inv )

par(mfrow=c(3,3))

hist(beta_draws[, 1], 100, xlab = "Intercept", main = "")
hist(beta_draws[, 2], 100, xlab = "PowerSeller", main = "")
hist(beta_draws[, 3], 100, xlab = "VerifyID", main = "")
hist(beta_draws[, 4], 100, xlab = "Sealed", main = "")
hist(beta_draws[, 5], 100, xlab = "Minblem", main = "")
hist(beta_draws[, 6], 100, xlab = "MajBlem", main = "")
hist(beta_draws[, 7], 100, xlab = "LargNeg", main = "")
hist(beta_draws[, 8], 100, xlab = "LogBook", main = "")
hist(beta_draws[, 9], 100, xlab = "MinBidShare", main = "")

par(mfrow=c(1,1))

# c

Metropolis <- function(c, N, summa, postFunction, ...) {
  proposal_last <- as.vector(rep(0, dim(X)[2]))
  thetas <- matrix(0, N, dim(X)[2])
  acceptance <- 0
  for (x in 1:N) {
    proposal <- as.vector(rmvnorm(1, mean = proposal_last, sigma = c * summa))
    p_theta <- postFunction(proposal, ...)
    p_theta_last <- postFunction(proposal_last, ...)
    alpha <- min(1, exp(p_theta - p_theta_last))
    r <- runif(1)
    if (alpha > r) {
      proposal_last <- proposal
      acceptance <- acceptance + 1
    }
    thetas[x, ] <- proposal_last
    
  }
  cat("Acceptance Ratio: ", acceptance / N)
  return (thetas)
}
n <- 10000

draws <- Metropolis(0.5, n, j_inv, LogPoisson, y, X, mu, Sigma)

cum_mean <- matrix(0, n, dim(X)[2])

for (i in 1:dim(draws)[2]) {
  cum_mean[,i] <- cumsum(draws[,i]) / seq_along(draws[,i]) 
}

par(mfrow=c(3,3))

plot(draws[,1], main="Intercept", xlab="iteration", ylab="Value", type="l")    
plot(draws[,2], main="PowerSeller", xlab="iteration", ylab="Value", type="l")  
plot(draws[,3], main="VerifyID", xlab="iteration", ylab="Value", type="l") 
plot(draws[,4], main="Sealed", xlab="iteration",ylab="Value", type="l") 
plot(draws[,5], main="Minblem", xlab="iteration",ylab="Value", type="l") 
plot(draws[,6], main="MajBlem", xlab="iteration", ylab="Value", type="l") 
plot(draws[,7], main="LargNeg", xlab="iteration", ylab="Value", type="l") 
plot(draws[,8], main="LogBook", xlab="iteration", ylab="Value", type="l") 
plot(draws[,9], main="MinBidShare", xlab="iteration", ylab="Value", type="l") 

plot(cum_mean[,1],  main="Intercept", xlab="iteration", ylab="Cumulative Mean", type="l")  
abline(h=model$coefficients[1], col="red")
plot(cum_mean[,2], main="PowerSeller", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[2], col="red")
plot(cum_mean[,3], main="VerifyID", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[3], col="red")
plot(cum_mean[,4], main="Sealed", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[4], col="red")
plot(cum_mean[,5], main="Minblem", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[5], col="red")
plot(cum_mean[,6], main="MajBlem", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[6], col="red")
plot(cum_mean[,7], main="LargNeg", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[7], col="red")
plot(cum_mean[,8], main="LogBook", xlab="iteration", ylab="Cumulative Mean", type="l")
abline(h=model$coefficients[8], col="red")
plot(cum_mean[,9], main="MinBidShare", xlab="iteration", ylab="Cumulative Mean", type="l") 
abline(h=model$coefficients[9], col="red")

par(mfrow=c(1,1))

betas_mean <- as.vector(draws[n-1,])
#d
cov = as.vector(c(1,1,1,1,0,0,0,1,0.5))

lambda = exp(cov%*%betas_mean)
nBids_sim = rpois(10000, lambda)

barplot(table(nBids_sim)/n, ylab="Probability", xlab="nBids")


