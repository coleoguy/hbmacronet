library(BMhyd)
library(diversitree)
library(phytools)
source("../Rscripts/getMS.R")
source("../Rscripts/macroevolve.R")

##Here we are taking trees that were made using seqgen sequences and trying to find thier speciation and extinction rates


##Read everything we need to do a downstream analysis
input_trees<-read_input_trees(location = "../input/trees",ntrees)
ntrees<-5 ##number of trees to read in

##Where the MCMCs are going to get saved
save_location<-"../outputs/clade/50_10/nj_mcmc_split"  
dir.create(save_location,showWarnings = FALSE)
for(i in 1:ntrees){
  dir.create(paste(save_location,"/tree",i,sep = ""),showWarnings = F)
}

# This is for reading in posterior sets
clade_trees<-readnfix(location = "../outputs/clade/50_10/mb",ntrees)


# # This is for reading in nj trees
# clade_trees<-read_input_trees(location = "../outputs/clade/50_10/nj/trees",ntrees,rename = F)
# #Now we need to make put this in a funky format so it works with the rest of the code - we're going to nest each element of this list in another list and then make a list of those lists so that way it is in a similar format to the posterior trees
# for(i in 1:ntrees){
#   clade_trees[[i]]<-list(clade_trees[[i]])
# }



#Read Hybrid records
##TODO fix the way the records are saved so I don't need to format them so much upon reading
hybrids<-list()
for(i in 1:ntrees){
  hy<-read.csv(paste("../hybridizations/clade/50_10/hybrid",i,".csv",sep = ""))[-1]
  hy[,5]<-max(branching.times(input_trees[[i]]))-hy[,3]
  hy[,3]<-hy$str
  hy[,4]<-hy[,5]
  hy[,1]<-as.character(hy[,1])
  hy[,2]<-as.character(hy[,2])
  hybrids[[i]]<-hy
  colnames(hybrids[[i]])<-c("donor","recipient","m","time.from.root.donor","time.from.root.recipient")
}


#Find out which node we need to do a split bd model
clades<-find_hybrid_clade(input_tree = input_tree,
                   clade_trees = clade_trees,
                   hybrids =hybrids,
                   tree_amount = ntrees) 
clades<-clades[[4]] ## we only care about the nodes that were identified to contain the hybrid clade



##Split MCMC
##bd mcmc for clade trees
for(i in 1:ntrees){ ##run through each set of trees 
  
  for(j in 1:length(clade_trees[[i]])){ ## Run through each posterior tree
    tree<-clade_trees[[i]][[j]]
    tree$edge.length<-tree$edge.length/max(branching.times(tree)) ## make tree unit length
    lk<-make.bd.split(tree = tree, nodes= clades[[i]][[j]],split.t=Inf)
    y<-mcmc(lik=lk,c(1,0.1,1,0.1),nsteps = 1000,w=0.1,print.every = 0,save.every=1000,save.file=paste(save_location,"/tree",i,"/tree",i,"mcmc",j,".csv",sep=""))
    w <- diff(sapply(y[2:3], range)) ##recalculate w to get a better runtime
    print(paste(i,j))
  }
}

