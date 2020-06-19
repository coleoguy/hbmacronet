

library(BMhyd)
library(diversitree)
library(phytools)




##read control split mcmcs - nj
control_mcmc<-rep( list(list()), ntrees ) 
for(i in 1:ntrees){
  for(j in 1:1){
    control_mcmc[[i]][[j]]<-read.csv(paste("../outputs/control/nj_mcmc_split/tree",i,"/tree",i,"mcmc",j,".csv",sep=""))
  }
}
##read clade split mcmcs - nj
clade_mcmc<-rep( list(list()), ntrees ) 
for(i in 1:ntrees){
  for(j in 1:1){
    clade_mcmc[[i]][[j]]<-read.csv(paste("../outputs/clade/50_10/nj_mcmc_split/tree",i,"/tree",i,"mcmc",j,".csv",sep=""))
  }
}




i<-1
j<-i

##Plot the input tree with the hybridizations added 
PlotNetwork(input_trees[[3]],flow=hybrids[[i]],cex=0.5)
bfit[[i]]

##Plot two trees side by side to compare them
##With these types of plots we need the concensus trees
co<-cophylo(control_trees[[i]][[j]],clade_trees[[i]][[j]])
plot(co,cex=0.5)

plot(control_mcmc[[1]][[1]]$p)
i<-2
##profiles plot of the MCMCs
profiles.plot(data.frame(control_mcmc[[i]][[j]]$lambda.1-control_mcmc[[i]][[j]]$mu.1,control_mcmc[[i]][[j]]$lambda.2-control_mcmc[[i]][[j]]$mu.2),col.line = c("red","blue"),main="control")
legend("topright",legend = c("r1","r2"),fill=c("red","blue"))
profiles.plot(data.frame(clade_mcmc[[i]][[j]]$lambda.1-clade_mcmc[[i]][[j]]$mu.1,clade_mcmc[[i]][[j]]$lambda.2-clade_mcmc[[i]][[j]]$mu.2),col.line = c("red","blue"),main="50_10c")
legend("topright",legend = c("r1","r2"),fill=c("red","blue"))

##lineage through time plot
mltt.plot(control_trees[[i]],legend = FALSE) 
title(main = "control")

mltt.plot(clade_trees[[i]],legend = FALSE) 
title(main = "clade")
