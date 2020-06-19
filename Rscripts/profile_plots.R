##This file shows how we can use MCMCs and paste them together
##Also shows how we take burn-in
##File locations may need to be adjusted



library(dplyr)
ntrees<-5
nposterior<-2501
treeburn<-0.25 ## percentage of trees to toss as burn in
start<-ceiling(treeburn*nposterior)


##read control split mcmcs 
control_mcmc<-rep( list(list()), ntrees ) 
for(i in 1:ntrees){
  for(j in start:nposterior){
    control_mcmc[[i]][[j]]<-read.csv(paste("../outputs/control/mcmc_split/tree",i,"/tree",i,"mcmc",j,".csv",sep=""))
  }
}
##read clade split mcmcs 
clade_mcmc<-rep( list(list()), ntrees ) 
for(i in 1:ntrees){
  for(j in start:nposterior){
    clade_mcmc[[i]][[j]]<-read.csv(paste("../outputs/clade/50_10c/mcmc_split/tree",i,"/tree",i,"mcmc",j,".csv",sep=""))
  }
}


##bind  mcmcs together
bound_control<-list()
bound_clade<-list()
for(i in 1:ntrees){
  bound_control[[i]]<-bind_rows(control_mcmc[[i]])
  bound_clade[[i]]<-bind_rows(clade_mcmc[[i]])
}

i<-2

profiles.plot(data.frame(bound_control[[i]]$lambda.1-bound_control[[i]]$mu.1,bound_control[[i]]$lambda.2-bound_control[[i]]$mu.2),col.line = c("red","blue"),main="control")
profiles.plot(data.frame(bound_clade[[i]]$lambda.1-bound_clade[[i]]$mu.1,bound_clade[[i]]$lambda.2-bound_clade[[i]]$mu.2),col.line = c("red","blue"),main="50_10c")


control_mcmc
  