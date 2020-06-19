library(ape)
source("./macroevolve.R")

readnjoin<-function(location,seq_amount,print.every=0){ ##read sequences and then use neighbor-joining to infer trees. Then save trees - will make new folder for trees
  ##Reads in sets of sequences and set them as DNAbin objects
  trees<-list()
  for(i in 1:seq_amount){
    seqs<-as.DNAbin(read.nexus.data(file = paste(location,"/dna_tree.",i,".nex",sep = "")))
    dist_matrix<-dist.dna(seqs) ## get distance matrix of sequences 
    t<-nj(dist_matrix) ## Do neighbor-joining to get tree
    t<-check_and_fix_ultrametric(t)#make ultrametric -fix rounding errors
    trees[[i]]<-t
    
    if(print.every!=0){ ##because I'm inpatient we can monitor the progress even though its pretty quick relative to other things that run in this project
      if(i%%print.every==0){
        print(i)
      }
    }
    
  }
  new_location<-paste(sub("/seq.gen","",location),"/nj",sep = "") ##new location to save data
  dir.create(new_location,showWarnings = FALSE) ##create new diectory if not there already
  write.nexus(trees,file=paste(new_location,"/trees",sep = "")) ## save trees in new folder
}

##Now we can can get neighbor-joining trees for the sequences we've already generated 
readnjoin("../outputs/control/seq.gen",100,print.every=1)
readnjoin("../outputs/clade/50_25c/seq.gen",2)
readnjoin("../outputs/clade/50_10/seq.gen",5)
readnjoin("../outputs/test/seq.gen",9)




## We can then work with trees as we would from MrBayes and do the things that are outlined in pipeline2
