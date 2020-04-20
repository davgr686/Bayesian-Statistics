# a

y <- c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
n <- 10 

tau2 = sum((log(y)-3.7)^2)/n


sig2_posterior  = (n-1) * tau2 /rchisq(10000, n-1)


hist(sig2_posterior, 200, xlim = c(0,1))


# b

G <- 2 * pnorm(sqrt(sig2_posterior)/sqrt(2), 0, 1) - 1
CI = quantile(G,probs=c(0.05,0.95))


h <- hist(G, breaks=100, plot=FALSE)
cuts <- cut(h$breaks, c(-Inf,CI[1],CI[2],Inf))
plot(h, col=c("white","green","white")[cuts])
abline(v=CI[1], col="green")
abline(v=CI[2], col="green")


kdens <- density(G)


sorted_kdens <- sort(kdens$y, index.return=TRUE, decreasing = TRUE)

target_dens <- sum(sorted_kdens$x) * 0.9

sum = 0

for (i in 1:length(sorted_kdens$x))
{
  sum = sum + sorted_kdens$x[i]
  if (sum > target_dens) {
    #breakpoint = kdens$y[sorted_kdens$ix[i]]
    breakpoint1 = kdens$y[sorted_kdens$ix[i-1]]
    break
  }
}
plot(kdens)
abline(h=breakpoint1)
