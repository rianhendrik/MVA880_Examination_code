library(readxl)
library(mixtools)
setwd("C:/Users/rianh/Dropbox/MVA880/EXAM/EXAM DATA")


#Problems 
#1) How do I know which mean belongs to which cluster - labelling issue

sample_sizes = c(200, 2000, 20000)


nsim = 10000
n = 3000
k = 2


#Generate a Gaussian mixture density

pis_t = c(1/2, 1/2)
mus_t = c(80, 80) #clearly distinguished
stds_t = c(3, 7)


delta = c(1:k)
deltas = sample(delta, n, replace = T, prob = pis_t)

data = matrix(NA, nrow=n, ncol=2)
for (i in 1:n){
  delta_i = deltas[i]
  data[i,1] = rnorm(1, mus_t[delta_i], stds_t[delta_i])
  data[i,2] = delta_i
}



hist(data[,1], nclass=100)


#Determining inital values with k-means
#The problem with this is that my inital values
#Will almost always be the same - and so by GMM results will always be the same
#Conditional upon the fact that I am using the same dataset.

# res = kmeans(data[,1],k)$cluster
# 
# mus = matrix(NA,nrow=k, ncol=d)
# sigmas = matrix(nrow=k, ncol=d)
# pis = matrix(nrow=k, ncol=d)
# for (i in 1:k){
#   mus[i] = mean(data[res == i,1])
#   sigmas[i] = sd(data[res == i,1])
#   pis[i] = sum(res == i)/length(res)
# }


lambdas = matrix(NA, nrow=nsim, ncol=k)
means = matrix(NA, nrow=nsim, ncol=k)
deviations = matrix(NA, nrow=nsim, ncol=k)
prop_missclass = matrix(NA, nrow=nsim, ncol=1)

for (iter in 1:nsim){
  

  conv = F
  while (conv == F){
  
  res = kmeans(data[,1],k)$cluster
  d = 1
  mus = matrix(NA,nrow=k, ncol=d)
  sigmas = matrix(nrow=k, ncol=d)
  pis = matrix(nrow=k, ncol=d)
  for (i in 1:k){
    mus[i] = mean(data[res == i,1])
    sigmas[i] = sd(data[res == i,1])
    pis[i] = sum(res == i)/length(res)
  }
  
# 
#     pis = matrix(NA, nrow=k, ncol=1)
#     unis = runif(k)
#     for (i in 1:k){
#       pis[i] = unis[i]/sum(unis)
#     }
# 
#     mus = matrix(NA, nrow=k, ncol=1)
#     for(i in 1:k){
#       mus[i] = runif(1, min(data[,1]), max(data[,1]))
#     }
# 
#     sigmas = matrix(NA, nrow=k, ncol=1)
#     for (i in 1:k){
#       sigmas[i] = runif(1, 1, 10)
#     }

    
    mus = sort(mus)
    sigmas = sort(sigmas)
    pis = sort(pis)
    
    s = Sys.time()
    gmm = normalmixEM(data[,1], lambda = pis , mu = mus, sigma = sigmas)
    f = Sys.time()
    td = as.numeric(f-s)
    
    if (td<0.5){
      conv = T
    }
  
  
  
  hard_clusts = matrix(NA, nrow=n, ncol=1)
  for (i in 1:n){
    hard_clusts[i] = which(gmm$posterior[i,] == max(gmm$posterior[i,]))
  }
  
  
  
  lambdas[iter,] = gmm$lambda
  means[iter,] = gmm$mu
  deviations[iter,] = gmm$sigma

  #testing classification ability
  sorted = sort(gmm$sigma)
  gmm_labels = matrix(NA, nrow=k, ncol=1)
  for (i in 1:k){
    compi = sorted[i]
    gmm_labels[i] = which(gmm$sigma == compi)
  }
  
  #now, kth entry of gmm_label, is the label that was given for the kth component. 
  #if gmm_labe[2] = 1, the second component was lablled as the first. 
  #we need to change that
  
  hard_clusts_true = matrix(NA, nrow=n, ncol=1)
  for (i in 1:n){
    gmm_comp = hard_clusts[i]
    hard_clusts_true[i] = which(gmm_labels == gmm_comp)
  }
  
  cbind(hard_clusts, hard_clusts_true, data[,2])
  
  missclass = hard_clusts_true-data[,2]
  prop_missclass[iter] = length(which(missclass != 0))/n
  }
}

#Sorting means and deviations into components
means= cbind(t(apply(means[,1:ncol(means)], 1, sort)))
deviations = cbind(t(apply(deviations[,1:ncol(deviations)], 1, sort)))

#Calculating expected value of each component
colMeans(means)
mus_t
colMeans(deviations)
stds_t
colMeans(lambdas)

mean_bias = mus_t - colMeans(means)
deviation_bias = stds_t - colMeans(deviations)
pis_bias = pis_t - colMeans(lambdas)
mean_bias
deviation_bias
pis_bias

mean(prop_missclass)
median(prop_missclass)
