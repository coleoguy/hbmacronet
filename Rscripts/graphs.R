
load(file="graph.RData") ##This compares 50_10c to control
library(BMhyd)
library(diversitree)
library(phytools)


i<-1
j<-i

##Plot the input tree with the hybridizations added 
PlotNetwork(input_trees[[i]],flow=hybrids[[i]],cex=0.5)
bfit[[i]]

##Plot two trees side by side to compare them
##With these types of plots we need the concensus trees
co<-cophylo(rec_trees[[j]],con_trees[[j]])
plot(co,cex=0.5)


##Clade credability
##With these types of plots we need the concensus trees
{plot.phylo(con_trees[[i]],cex=0.5,main="control concensus")
  nodelabels(con_trees[[i]]$node.label,cex=0.5)}

{plot.phylo(con_trees[[i]],cex=0.5,main= "clade concensus")
  nodelabels(con_trees[[i]]$node.label,cex=0.5)}

##profiles plot of the MCMCs
profiles.plot(data.frame(bound_control[[i]]$lambda.1-bound_control[[i]]$mu.1,bound_control[[i]]$lambda.2-bound_control[[i]]$mu.2),col.line = c("red","blue"),main="control",xlim = c(-1,8))
legend("topright",legend = c("r1","r2"),fill=c("red","blue"))
profiles.plot(data.frame(bound_clade[[i]]$lambda.1-bound_clade[[i]]$mu.1,bound_clade[[i]]$lambda.2-bound_clade[[i]]$mu.2),col.line = c("red","blue"),main="50_10c",xlim = c(-1,8))
legend("topright",legend = c("r1","r2"),fill=c("red","blue"))

##lineage through time plot
mltt.plot(control_trees[[i]],legend = FALSE) 
title(main = "control")

mltt.plot(clade_trees[[i]],legend = FALSE) 
title(main = "clade")
