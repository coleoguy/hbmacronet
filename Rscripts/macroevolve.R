check_and_fix_ultrametric <- function(phy){
  # This function checks trees to see if they pass ape ultrametricity test.
  # If not, it computes the differential root-to-tip distance across all tips.
  # It adds the appropriate quantity to each terminal branch length to ensure that 
  # tree passes ultrametric test.
  # Note: this is only a valid method of making trees ultrametric when the 
  # 	non-ultrametricity is due to small numerical discrepancies, e.g., 
  #   rounding or other floating point issues during phylogeny construction.
  # 
  # written by Dan Rabosky
  
  if (!is.ultrametric(phy)){
    
    vv <- vcv.phylo(phy)
    dx <- diag(vv)
    mxx <- max(dx) - dx
    for (i in 1:length(mxx)){
      phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
    }
    if (!is.ultrametric(phy)){
      stop("Ultrametric fix failed\n")
    }	
  }
  
  return(phy)
}

readnfix<-function(location,tree_amount){
  ##Reads in a set of posterior trees and attempts to make them ultrametric
  trees<-list()
  for(i in 1:tree_amount){
    trees[[i]]<-read.nexus(file = paste(location,"/tree",i,".run1.T",sep = "")) ##NOTE: currently only reads one of two chains produced by mrBayes
    for(j in 1:length(trees[[i]])){
      trees[[i]][[j]]<-check_and_fix_ultrametric(trees[[i]][[j]])
    }
    print(i)
  }
  return(trees)
}

read_input_trees<-function(location,tree_amount,rename=T){
  ##Reads in input trees and attempts to rename them to match the names given by getMS()
  input<-read.nexus(location)
  input_trees<-list()
  for( i in 1:tree_amount){
    input_trees[[i]]<-input[[i]]
    if(rename){
    input_trees[[i]]$tip.label<-getMS(input_trees[[i]],1,"T",new.names = T)[[2]][,2] ##rename the input trees to match how they are changed in getMS() ##TODO- rename trees when saving them as oposed to when reading them
    }
    input_trees[[i]]$tip.label<-as.character(input_trees[[i]]$tip.label)
  }
  return(input_trees)
}

find_hybrid_clade<-function(input_tree,hybrids,clade_trees,tree_amount){
  ##attempts to identify the hybridized clade from the input tree and hybrid records
  ##Then tries to find this clade on one of the posterior trees 

  #find smallest node that fits with all hybrid events 
  bfit<-list()
  bnode<-list()
  for(i in 1:tree_amount){ ##Look at every input tree
    hy<-hybrids[[i]]
    tree<-input_trees[[i]]
    hyspecies<-as.numeric(unique(c(hy$donor,hy$recipient)))##get all recorded tips that are involved in events
    bfit[[i]]<-rep(NA,tree$Nnode*2)
    for(j in (tree$Nnode+2):(tree$Nnode+50)){ ## look at each node in the tree
      desc_indices<-as.numeric(get.descendants(j,tree)) ##Find all descendents
      temp<-as.numeric(tree$tip.label[desc_indices[desc_indices<=length(tree$tip.label)]]) ## find all tip labels that are a descendent
      if(all(hyspecies %in% temp) & (length(temp)<length(bfit[[i]])) ){ ##if node has all species in the record and is smaller than current best node
        bfit[[i]]<-temp
        bnode[[i]]<-j
      }
    }
  }
  #Tries to find a clade that is close to the one identified from the input hybrid events
  #Nodes will be scored with a +1 for having a correct species and -1 for having an incorrect species
  #if two nodes tie then the node with more species will be taken
  bclade_tips<-rep(list(list()),length(clade_trees))
  bclade_node<-rep(list(list()),length(clade_trees))
  for(i in 1:tree_amount){ ## Do this for each group of output trees
    for(j in 1:length(clade_trees[[i]])){ ##Do this for each posterior output tree of a group
      tree<-clade_trees[[i]][[j]]
      bscore<-(-Inf)
      for(k in (length(tree$tip.label)+1):(2*tree$Nnode+(length(tree$tip.label)-tree$Nnode))){##Look at each node in the tree
        ##get all tip labels for a node 
        desc_indices<-as.numeric(get.descendants(k,tree))
        desc<-as.numeric(tree$tip.label[desc_indices[desc_indices<=length(tree$tip.label)]])
        
        ##score tips based off of input hybridizations
        addition<-sum(desc %in% bfit[[i]]) ## tips that were identified by the input
        subtract<-sum(!(desc %in% bfit[[i]])) ##tips that were not identified by the input
        score<-addition-subtract ## score of node 
        if(score>=bscore){
          if(score==bscore && length(bclade_node[[i]][[j]]>desc)){
          }
          else{
            bclade_tips[[i]][[j]]<-desc
            bclade_node[[i]][[j]]<-k
            bscore<-score
          }
        }
      }
    }
  }
  ##return a list of:
    ##bfit - all species as identified by the input tree - NOTE: The numbers here are the tiplabels NOT the indicies of the tips
    ##bnode- The nodeof the input tree that results in bfit
    ##bclade_tips - all species found in the posterior tree that result in a clade close to bfit - NOTE: The numbers here are the tiplabels NOT the indicies of the tips
      ##this is a nested list with bclade_tips[[i]][[j]] - "j" refers to a specific posterior tree while "i" refers to a set of posterior trees, each set corresponding to an input tree of the same index
    ##bclade_node - the node that results in bclade_tips
      ##This is a nested list in the same format as bclade_tips- bclade_node[[i]][[j]]
  return(list(bfit,bnode,bclade_tips,bclade_node))
}