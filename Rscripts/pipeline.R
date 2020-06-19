
# This script will slowly be built up to become our primary
# pipeline for the production of posterior distribtuions of 
# trees and the associated DNA alignement they were inferred from
#
#
# This is just a reporting flag when turned on the script will save
# intermediate results for record keeping or troubleshooting
reporting <- F

output_dir<-"../outputs/clade/50_10" ## Directory to upload records from pipeline
dir.create(output_dir)
dir.create(paste(output_dir,"ms",sep=""))
dir.create(paste(output_dir,"mb",sep=""))
dir.create(paste(output_dir,"seq.gen",sep=""))
#
#############
# 1 Diversitree - produces phylo trees
#############
library(diversitree)
# # generate new trees
# # here we generate 100 trees each with 50 taxa 
# # using a speciation rate of 1 and extinction rate of .1
taxa <- 50   # the number of taxa on each tree
ntrees <- 5  # the number of simulated trees
# 
# trees <- trees(pars = c(1, .1), max.taxa = taxa, type = "bd", n = ntrees)
# for (i in 1:ntrees){
#   trees[[i]]$tip.label<-1:50
# }
# rm(taxa)


# read previously generated trees 
trees<-read.nexus('../input/trees')

###############
# 2 phylo -> MS - produces strings for MS
###############
source("getMS.R")
source("addHybrid.R")

ngenes <- 20  # number of gene tree samples to produce
ntrees<-100
ms.strings <- vector()
for(i in 1:ntrees){
  ms_tree<- getMS(tree = trees[[i]], samps = ngenes, report = "T") ##get phylo trees into ms strings
  
  tree_length<-diag(vcv.phylo(trees[[i]]))[1] ## get tree length - used in introgression
  
  ## control
  # ms.strings[i]<-ms_tree
  
  ## clade specific introgression - smallest clade of atleast 25% of species
  clade<-find_clade(get_ms_frame(ms_tree)[[2]],0.25)
  out<-addHybrid(ms_tree,strength = 0.50,
                 Nevents = 10,
                 clade=clade,
                 report=TRUE,
                 type = 'point')
  ms.strings[i]<-out[[1]]
  write.csv(out[[2]],file=paste('../input/hybridizations/hybrid',i,'.csv',sep=''))
  print(i)
  
}
if(reporting == T) write(ms.strings, file = "ms.strings.txt")

#############
# 3 MS - produces samples of possible gene trees
#############
# when we start calling other programs you will need to make
# sure that the executable is in your system path
bp <- 400# the number of nuceleotides we want to simulate in seqGen

for(j in 1:ntrees){
  z <-system(command = paste("ms ", ms.strings[j],
                             "|tail -n +5 | grep -v // ",
                             sep = ""), intern=T)
  # this next loop adds a bit of text in front of tree to
  # define the number of bp to simulate based on it
  for(i in 1:length(z)){ 
    if(z[i] != "") z[i] <- paste("[", bp, "]", z[i], sep = "")
  }
  write(z, file=paste(output_dir,"/ms/gene.trees.tree.", j, ".txt", sep = ""))
}

#############
# 4 SeqGen
#############
alignment.length <- 8000 # total length of desired alignment
# rather than reformatting the data it is now saved in nexus format
# MrBayes should handle this without reformatting
for(i in 1:ntrees){
  system(command = paste("seq-gen -mHKY -l ",alignment.length,
                         " -s 0.01 -on -p ",
                         ngenes, 
                         " < ",output_dir,"/ms/gene.trees.tree.", i,".txt > ../output/seq.gen/dna_tree.", i, ".nex", sep=""))
}


############
# 5 MrBayes - produces a posterior distribution of trees
############
# we will create mb block for each run that will have all settings 
# for a given run
gens <- 25000 # generations for mcmc

for(i in 1:ntrees){
  mb.block <- paste("begin mrbayes;\n",                 # identifies thsi as control block
                    "set autoclose=yes nowarn=yes;\n",  # close MB when done
                    "execute ../output/seq.gen/dna_tree.", i, ".nex;\n",   # this is our DNA from seqgen
                    "lset nst=2 rates=gamma;\n",        # sets of the DNA model
                    "prset brlenspr=clock:uniform;\n",  # make our tree ultrametric
                    "mcmc ngen=", gens, " samplefreq=10 file=../output/mb/tree", i, ";\n",  # sets up the mcmc and save files
                    "sumt Conformat=simple;\n",                          # summarizes trees
                    "end;\n", sep = "")
  write(mb.block, file=paste(output_dir,"/mb/mb.block",sep = ""))
  system(command = paste("mb ",output_dir,"/mb/mb.block",sep = ""))
}

####################
# process trees
####################

# need to still write some code to take final trees
# downsample to 100 trees for trait analysis

