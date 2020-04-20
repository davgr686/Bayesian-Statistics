# A
s <- 5
n <- 20
f <- n - s
a0 <- 2
b0 <- 2
set.seed(12345)

posterior <- rbeta(100000, a0 + s, b0 + f)

means <- c()
vars <- c()
indx <- 0
for (i in seq(1,100000, 100))
{
  post <- sample( posterior, i, replace=TRUE )
  means[indx] <- mean(post)
  vars[indx] <- sd(post)
  indx <- indx + 1
  
}

truemean <- (a0 + s) / (a0 + s + b0 + f)
truevar <- (a0 + s) * (b0 + f) / ((b0 + f + a0 + s)**2 * (b0 + f + a0 + s + 1))


plot(means, type="l")
abline(h=truemean, col = 'red')
plot(vars**2, type="l")
abline(h=truevar, col = 'red')

## B
post <- sample( posterior, 10000, replace=TRUE )
sum(post > 0.3) / length(post)


prob_theta_bigger_than_03 <- pbeta(0.3, a0 + s, b0 + f, lower.tail = FALSE)


## C

phis <- density(log(post/(1-post)))
plot(phis)



