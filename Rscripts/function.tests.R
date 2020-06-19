# This script just shows how we can slice trees and get relatedness out of
# data

library(diversitree)
# make some trees
trees <- trees(pars = c(1, .1), max.taxa = 10, type = "bd", n = 100)

# here we convert a tree to MS and return the MS string a modified edge
# matrix with deme names for every branch
MSdata <- getMS(tree = trees[[1]], samps = 20, report = "T")

# next we can pull out a tree that represent the relatedness of lineages 
# at any point in the past
t <- .5*max(branching.times(trees[[1]]))
deme.tree <- HybProbs(edges=MSdata[[2]], time=t)

# lets look at our trees
par(mfcol=c(1,2))
plot(trees[[1]])
nodelabels(frame="circle", cex=.5)
t <- .5*max(branching.times(trees[[1]]))
abline(v=max(branching.times(trees[[1]]))-t,col="red",lwd=2)

plot(deme.tree)
nodelabels(frame="circle", cex=.5)

# lets draw 1000 hybridization events and make sure that 
# the distribution makes sense
result <- vector()
for(i in 1:1000){
  result[i] <- GetHybrid(deme.tree)
}

table(result)
# thats kinda beautiful

pair <- GetHybrid(deme.tree)

# now we draw an MS component to add to our string
MakeMsHyb(pair=pair, time=t, edges=MSdata[[2]], strength=.5)




