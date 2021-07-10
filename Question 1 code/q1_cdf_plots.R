library(readxl)
setwd("C:/Users/rianh/Dropbox/MVA880/EXAM/EXAM DATA")

x = read_xlsx("dist.xlsx")
x = as.numeric(unlist(x))

pi1n = 0.3
pi2n = 0.25
pi3n = 0.45

mu1n = 140
mu2n = 150
mu3n = 200

sigma1n = 3
sigma2n = 4
sigma3n = 20

pi1e = 0.27
pi2e = 0.31
pi3e = 0.42

mu1e = 137.6
mu2e = 153.2
mu3e = 203.2

sigma1e = 3
sigma2e = 6.5
sigma3e = 16.1


x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=5000)


y_n_cdf <- pi1n * pnorm(x1, mean=mu1n, sd=sigma1n) + pi2n * pnorm(x1, mean=mu2n, sd=sigma2n) + pi3n * pnorm(x1, mean=mu3n, sd=sigma3n)


y_e_cdf <- pi1e * pnorm(x1, mean=mu1e, sd=sigma1e) + pi2e * pnorm(x1, mean=mu2e, sd=sigma2e) + pi3e * pnorm(x1, mean=mu3e, sd=sigma3e)

y_n_cdf

plot(ecdf(x), main='Empirical CDF plot of question 1 data, with Null and GMM CDFs',
     cex.main=1.5)
lines(x1, y_n_cdf, col="brown1", lwd=2, lty=2)
lines(x1, y_e_cdf, col='chartreuse1', lwd=2)
legend('bottomright', col=c('black', 'brown 1', 'chartreuse1'), lwd=2, lty = c(4, 2, 1), 
       legend=c("Empirical CDF", "Null CDF", "EM GMM CDF"), cex=0.7, bty = 'n')
