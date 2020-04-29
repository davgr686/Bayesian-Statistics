library("mvtnorm")
set.seed(1234)
Data<-read.table('WomenWork.dat', header=T)  
y <- as.vector(Data[,1])
X <- as.matrix(Data[,2:9])
covNames <- names(Data)[2:length(names(Data))]
nPara <- dim(X)[2];

tau <- 10
mu <- as.vector(rep(0,nPara)) 
Sigma <- tau^2*diag(nPara)


LogPostLogistic <- function(betas,y,X,mu,Sigma){
  nPara <- length(betas)
  linPred <- X%*%betas
  
  logLik <- sum( linPred*y -log(1 + exp(linPred)))
  if (abs(logLik) == Inf) logLik = -20000
  
  logPrior <- dmvnorm(betas, matrix(0,nPara,1), Sigma, log=TRUE)
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,dim(X)[2])); 

OptimResults<-optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)


beta_tilde <- OptimResults$par
j_inv <- -solve(OptimResults$hessian)

betas <- rmvnorm(10000, beta_tilde, j_inv)

hist(betas[,7], 100)
cred_interval <- quantile(betas[,7], c(0.025, 0.975))
abline(v=cred_interval[1])
abline(v=cred_interval[2])


glmModel <- glm(Work ~ 0 + ., data = Data, family="binomial")
summary(glmModel)

#b

covariates <- matrix(c(1, 10, 8, 10, 1, 40, 1, 1), 1,8)
prob <- c()

for (i in 1:10000) {
  betas <- rmvnorm(1, beta_tilde, j_inv)
  prob[i] <- exp(covariates%*%t(betas))/(1+exp(covariates%*%t(betas)))
}

hist(prob, 50)

# c
working_women <- c()

for (i in 1:length(prob)) {
  working_women[i] <- rbinom(1, 10, prob[i])
}


hist(working_women, 50)
