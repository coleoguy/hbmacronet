

addHybrid<-function(ms_input,timespan=NA,Nevents,strength,clade=NA,type="point",report=FALSE,doubles=TRUE,tolerance=NA){
  ##ms_input - a String in the format of ms - usually the output of getMS()
  ##timespan - a vector of times in the format of c(start,end)
  ##Nevents  - a natural number for how many hybridization events you want 
  ##Strength - a number where 0< x< 1 that corresponds to the fraciton of gene trees that follow the hybridized history
  ##Clade    - Used in Clade specific hybridization. This is the node that defines the clade for hybridization. NOTE: the node number given corresponds to the node of the phylo object version of the tree
  ##type     -Type of hybridization. "Point" or "continuous." Point specifies hybridization at a given time point whereas continuous specifies continual hybridization from the corresponding node up until the time point
  
  
  ms<-get_ms_frame(ms_input) ##put data into a more workable format
  ##TODO Make it possible to do a time period AND a clade rather than OR
  if(!is.na(clade)){ ##add hybrids based on the clade rather than by a time period
    if(!is.na(timespan)){
      message("ignoring given timespan and using clade birth to tips")
    }
    ms_frame<-clade_slice(ms[[2]],node = clade)
    begin<-0
    end<-max(ms_frame$time)
  }
  else{ ## add hybrids based on a time period rather than a clade
    begin<-timespan[1]
    end<-timespan[2]
    ms_frame<-ms[[2]]
  }
  
  ##prepare variables and vectors for storing hybrid information
  ej<-c(rep("-ej",Nevents))
  es<-c(rep("-es",Nevents))
  em<-c(rep("-em",Nevents))
  time<-c()
  sp1<-c()
  sp2<-c()
  str<-rep(strength,Nevents)
  temp_sp<-c()
  npop<-max(c(max(ms[[2]]$sp1),max(ms[[2]]$sp2)))
  
  time<-sort(runif(Nevents,begin,end)) ##find times of events, specified by start and end time 
  i<-1
  tol<-1
  ##add in hybridization events up to intensity amount
  while(i<=Nevents){
    include_pair<-TRUE ##denotess whether to include the given pair when they are found. Will only change to false if doubles=F and there's a repeat
    
    ##cut the tree to be at a time found earlier
    sliced<-timeslice(ms_frame,time[i])
    
    ##weight pairs for hybridiszation and find a pair
    joins<-find_progeny(sliced) ##find all progeny for each node
    joins<-new_join(sliced,joins) ##find when to progeny have a common ancestor, used in weighting hybridization events
    ##find a pair for hybrid
    pair<-joins[sample(seq(nrow(joins)),1,replace=TRUE,prob=joins$weight),]
    if(!doubles){
      ##check to see if pairing already exists 
      p<-pair[c(1,2)]
      
      if((sp1 %in% p) && (sp2 %in% p)){
        ##pair already exists
        include_pair<-FALSE
        tol<-tol+1 
      }
      ##If attempts exceed tolerance then move on to another hybridization
      if(tol>=tolerance){
        tol<-1 ##resest tol for next time
        i<-i+1 ## Move on to next event
      }
    }
    if (include_pair){
      ##randomly order the hybrid pair 
      pair<-sample(pair[c(1,2)],2)
      temp_sp<-c(temp_sp,npop+i)
      sp1<-c(sp1,pair[1])
      sp2<-c(sp2,pair[2])
      i<-i+1
    }
  }
  ms_start<-paste(ms[[1]],collapse=" ")
  if(type=="point"){
    splits<-data.frame(es,time,unlist(sp1),str,stringsAsFactors = FALSE) ## make df of es statements
    joins<-data.frame(ej,time, temp_sp, unlist(sp2), stringsAsFactors = FALSE) ## make df of ej statements
    merged<-merge_in_time(ms[[2]],splits=splits,joins=joins)
  }
  if(type=="continuous"){
    hyb<-data.frame(em,time,unlist(sp1),unlist(sp2),str,stringsAsFactors = FALSE) ## put into nice data frame 
    merged<-merge_in_time(ms[[2]],splits = hyb)
  }
  saved<-data.frame(as.character(unlist(sp1)),as.character(unlist(sp2)),str,time,time,stringsAsFactors = FALSE)
  colnames(saved)<-c("donor","recipient","m","time.from.root.donor","time.from.root.recipient") ##saved in a goofy data frame so it can be graphed later. NOTE the type of hybridization is NOT recorded here, you'll need to do that elsewhere
  output<-paste(ms_start,merged,"-T" ,sep=" ")
  if(report){ ## if true then include a data frame that denotes the hybridizations
    output<-list(output,saved)
  }
  return(output)
}

timeslice<-function(ms,time){
  ##get rid of all events that happen past the given time
  slice<-ms[ms$time>=time,]
  ##reduce all other lengths by amount designated
  slice$time<-slice$time-time
  return(slice)
}

## finds all tips that are associated with a given node
## returns a list where the indexes relate to the rows of the given ms dataframe and the elements are the associated tips
##initiates recursive call of progeny function
find_progeny<-function(ms){
  ms_frame<-ms
  n<-ms_frame[which.max(ms_frame$time),] ## find node at root
  c<-list() ## make empty list of children at a given node. Will get filled in throughout recursion
  npop<-max(c(max(ms_frame$sp1),max(ms_frame$sp2))) ## find the number of species, used to get correct length of children
  c[[npop]]<-0
  node_relations<-progeny(ms = ms_frame,node=n,children=c) ## find progeny of the root 
  node_relations[[npop]]<-NULL
  return(node_relations)
}

##return a list that gives the tips associated with all nodes
##the index of the list corresponds with the row of the ms data frame
progeny<-function(ms,node,children){
  ##find left and right children for node of interest
  ms_frame<-ms
  left<-node$sp1
  right<-node$sp2
  
  l_node<-ms_frame[((ms_frame$sp1==left | ms_frame$sp2==left) & ms_frame$time<node$time),] ## find all nodes that contain left species that occur after our input node
  l_node<-l_node[which.max(l_node$time),]  ## the node with the greatest time is the one that occurs next
  l_pos<-which((ms_frame$sp2==l_node$sp2) & (ms_frame$sp1 == l_node$sp1)) ## row position of left node on dataframe
  
  
  r_node<-ms_frame[((ms_frame$sp1==right | ms_frame$sp2==right)& ms_frame$time<node$time),]  ## find all nodes that contain right species that occur after our input node
  r_node<-r_node[which.max(r_node$time),]  ## the node with the greatest time is the one that occurs next
  r_pos<-which(ms_frame$sp2==r_node$sp2 & ms_frame$sp1 == r_node$sp1) ## row position of right node on dataframe
  
  
  ##if the left branch has no child node
  if(length(l_pos)==0){
    l_child<-node$sp1 ## set left branch as the species 
  }
  ##if the left branch has child nodes then find progeny of the left node
  else{
    if(is.null(children[[l_pos]])){
      #children[[l_pos]]<-child(ms_frame,l_node,children)[[l_pos]]
      children<-cmbnlst(children,progeny(ms_frame,l_node,children)) ##combine children list with an updated children list that has gone further down in recursion
    }
    l_child<-children[[l_pos]]
  }
  
  ##if the right branch has no child node
  if(length(r_pos)==0){
    r_child<-node$sp2 ## set right branch as the species
  }
  ##if the right branch has child nodes then find progeny of the right node
  else{
    if(is.null(children[[r_pos]])){
      #children[[r_pos]]<-child(ms_frame,r_node,children)[[r_pos]]
      children<-cmbnlst(children,progeny(ms_frame,r_node,children)) ##combine children list with an updated children list that has gone further down in recursion
    }
    r_child<-children[[r_pos]]
  }
  pos<-which(ms_frame$sp2==node$sp2 & ms_frame$sp1 == node$sp1) ##find position of current node 
  children[[pos]]<-c(r_child,l_child) ## progeny of current node is the sum of progeny of left and right branches
  return(children)
}

##Combine two lists of same info in most elements but one list may have data where the other has NULL space
cmbnlst<-function(x,y){
  i<-1
  output<-list()
  while(i<=length(x)){
    if(is.null(x[[i]])){
      output[[i]]<-y[[i]]
    } 
    else{
      output[[i]]<-x[[i]]
    }
    i<-i+1
  }
  return(output)
}

##identify which node that two species first join and assign a weight based on the time of the node
new_join<-function(ms,node_relations){
  ##elements for the output
  sp1<-c()
  sp2<-c()
  time<-c()
  exist<-unique(c(ms$sp1,ms$sp2)) ## find which species currently exist
  ##iterate thru all pair combinations and find the minimum node that includes them both
  npop<-max(c(max(ms$sp1),max(ms$sp2)))
  first<-1
  while(first<npop){
    second<-first+1
    while(second<=npop){
      i<-1
      t<-c()
      s<-c()
      while(i<=length(node_relations)){
        ##look thru all node relations to find nodes that include the two species of interest
        ##get time of node
        if((first %in% node_relations[[i]]) & (second %in% node_relations[[i]])){
          t<-c(t,ms$time[i])
          s<-c(s,length(node_relations[[i]]))
        }
        else{
          t<-c(t,max(ms$time))
          s<-c(s,node_relations[[which.max(ms$time)]])
        }
        i<-i+1
      }
      ##only record pairs if the two species currently exist
      if((first %in% exist) & (second %in% exist)){
        sp1<-c(sp1,first)
        sp2<-c(sp2,second)
        time<-c(time,min(t))
      }
      second<-second+1
    }
    first<-first+1
  }
  size<-as.list(table(time))
  weight<-c()
  for(i in 1:length(time)){
    weight<-c(weight,1/time[i]/size[[as.character(time[i])]])
  }
  joins<-data.frame(sp1,sp2,time,weight)
  return(joins)
}

##Put information from the ms string into a more workable format
##ouputs a list where the first element is the first part of ms string while the second element is the workable data table
## TODO edit this list so there is a 3rd element that contains the "-T" or whatever is at the end of a ms string - currently the last bit is hacked in there at the end of addHybrid()
get_ms_frame<-function(tree){
  ms_net<-strsplit(tree, split=" ")[[1]] ##split the string by spaces 
  start<-min(grep("-ej",ms_net)) ## The starting things are the things before the first "-ej" statement, we don't do anything with this section until we want to put everything back together
  ms_start<-ms_net[1:start-1]
  ms_net<-ms_net[start:length(ms_net)] ## this is a vector of strings that we are interested in
  ms_net
  ##TODO Make more efficient, currently rewrites each column with each iteration
  i<-1
  type<-c()
  time<-c()
  sp1<-c()
  sp2<-c()
  while (i<length(ms_net)){
    type<-c(type,ms_net[i]) ## stores all "ej"
    time<-c(time,ms_net[i+1]) ## stores all times for each node
    sp1<-c(sp1,ms_net[i+2])## stores the first population in a node
    sp2<-c(sp2,ms_net[i+3]) ## stores the second population in a node
    i<-i+4
  }
  ms_frame<-data.frame(type,as.numeric(time),as.numeric(sp1),as.numeric(sp2),stringsAsFactors = FALSE) ## put it all nicely in a data frame 
  names(ms_frame)<-c("type", "time", "sp1",  "sp2")
  return(list(ms_start,ms_frame))
}


##takes the ms data frame and merges it with hybrid data frames
## The ms  string requires ms commands to be in the order the they happen on the tree so we do that here
##when doing "point" type hybridizations there are two ms commands that are used: es and ej. As a result two hybridization data frames are given, one for es(splits) and one for ej(joins)
merge_in_time<-function(original,splits,joins=NA){
  output<-""
  while((nrow(original)!=0)  | nrow(splits)!=0){ ##while our original data frame and our hybrid frame are not empty
    min_original<-original[which.min(original$time),]$time ##find the earliest original event
    min_split<-splits[which.min(splits$time),]$time ## find the earliest hybridization event
    
    if(length(min_original)==0){ ## if no more original events left
      min_original<-Inf
    }
    if(length(min_split)==0){ ## if no more hybridization events left
      min_split<-Inf
    }
    if(min_split<min_original){ ## the next earliest event is a hybridization event
      new_splits<-paste(splits[which.min(splits$time),],collapse=" ") ##Get the detials of that event
      splits<-splits[-which.min(splits$time),] ## delete that event from the data frame
      
      
      if(!is.na(joins)){ ## if there is a join data frame too then we we want that information as well
        new_joins<-paste(joins[which.min(joins$time),],collapse=" ")
        output<-paste(output,new_splits,new_joins)
        joins<-joins[-which.min(joins$time),]
      }
      else{
        output<-paste(output,new_splits) ## add earliest hybrid event to our output
      }
      
    }
    else{## the next earliest event is an original event
      new_original<-paste(original[which.min(original$time),],collapse=" ") ## get the info of the earliest original event
      output<-paste(output,new_original) ## add that information to our output
      original<-original[-which.min(original$time),] ##Delete that event from the data frame
    }
  }
  return(output)
}
clade_slice<-function(ms,node){## Slice the tree such that the only nodes and species that remain are those that are within the specified node
  relations<-find_progeny(ms) ## find the progeny of each node 
  i<-1
  sub_node<-c()
  ##get rid of all nodes that are not within the clade of interest
  while(i<=length(ms[,1])){
    if(all(relations[[i]] %in% relations[[node]])){
      sub_node<-c(sub_node,i)
    }
    i<-i+1
  }
  return(ms[sub_node,])
}


##find the smallest node that contains at least the proportion of total species as specified by the clade.size variable
find_clade<-function(ms,clade.size){
  ##clade.size is a proportion of the total species- so it should range from 0 to 1
  relations<-find_progeny(ms)
  j<-1
  sizes<-c()
  while(j<=length(relations)){ ##Look at each node and find the proportion of species that that node contains
    sizes<-c(sizes,length(relations[[j]]))
    j<-j+1
  }
  npop<-max(sizes)
  prop<-sizes/npop
  return(which(prop == min(prop[prop>=clade.size]))[1])
}
