## part 3

# a
set.seed(12345)
Y <- c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)

u <- 2.39

k <- seq(0.001,10,0.1)

posterior <- function(k) {
  return (prod(exp(k * cos(Y - u))/(2 * pi * besselI(k, 0))) * dexp(k))
}

post <- sapply(k, posterior)

plot(k, post, type='l')
abline(v=2.15, col="red")
