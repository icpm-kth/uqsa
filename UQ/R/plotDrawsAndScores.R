#load("/Users/fedmil/Documents/uqsa/Regularization/Experiment3_ns_5000_npc_5000.RData")
#load("/Users/fedmil/Documents/uqsa/Regularization/Experiment3_12_ns_5000_npc_5000.RData")
load("/Users/fedmil/Documents/uqsa/Regularization/Experiment3_12_18_ns_5000_npc_5000.RData")

par(mfrow=c(3,1))
plot(scores1)
plot(scores2)
plot(scores3)

draws<-draws3
par(mfrow=c(3,3))
for(j in 1:9){
  plot(draws[,j], main = paste("Par ", j, " ", parNames[j]))
  }
for(j in 10:18){
  plot(draws[,j], main = paste("Par ", j, " ", parNames[j]))
}
for(j in 19:27){
  plot(draws[,j], main = paste("Par ", j, " ", parNames[j]))
}

#HISTOGRAMS
for(j in 1:9){
  plot(hist(draws[,j]), main = paste("Par ", j, " ", parNames[j]))
}
for(j in 10:18){
  plot(hist(draws[,j]), main = paste("Par ", j, " ", parNames[j]))
}
for(j in 19:27){
  plot(hist(draws[,j]), main = paste("Par ", j, " ", parNames[j]))
}

# Number of unique samples
dim(unique(draws))

# Acceptance rate
dim(unique(draws))[1]/dim(draws)[1]

