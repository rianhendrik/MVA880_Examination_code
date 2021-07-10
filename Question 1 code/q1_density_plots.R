#GMM Plots for Animations
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


hist(x, prob=T, breaks=15, xlim=c(range(x)[1], range(x)[2]), 
     ylim = c(0, max(y_n)),
     main='A histogram of the data for question 1 \n with the null and 3 component EM GMM solutions overlayed',
     cex.main=1.5,
     col = ('cornflowerblue'))
x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
y_n <- pi1n * dnorm(x1, mean=mu1n, sd=sigma1n) + pi2n * dnorm(x1, mean=mu2n, sd=sigma2n) + pi3n * dnorm(x1, mean=mu3n, sd=sigma3n)
y_e <- pi1e * dnorm(x1, mean=mu1e, sd=sigma1e) + pi2e * dnorm(x1, mean=mu2e, sd=sigma2e) + pi3e * dnorm(x1, mean=mu3e, sd=sigma3e)
lines(x1, y_n, col="brown1", lwd=2, lty=2, ylim = max(y_n))
lines(x1, y_e, col='chartreuse1', lwd=2)
legend('topright', col=c('brown 1', 'chartreuse1'), lwd=2, lty = c(2, 1), 
       legend=c("k=3 Null GMM density", "k=3 EM GMM density"), cex=0.6, bty = 'n')
