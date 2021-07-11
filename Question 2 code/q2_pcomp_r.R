setwd("C:/Users/rianh/Dropbox/MVA880/EXAM/EXAM DATA")

d = read.csv("digitq2.csv")
digits = d[,2:257]



scores50 = prcomp(digits, scale = T)$x


write.csv(scores50)