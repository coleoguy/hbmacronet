source("getMS.R")
source("addHybrid.R")

##Lets make sure that the components of addHybrid work by using them in a somewhat trival exmaple
##First generate a tree of say, 10 species 
set.seed(1)
tree<-tree.bd(c(1,0.1),max.taxa = 10)

##we're going to want to  rename the tiplabels so they are similar to how the mstring will name them
##And for viewing ease we will do the same with how the nodes are labeled during the addHybrid process
##When we have a list of the nodes in the ms format they are ordered such that most recent nodes from the tips appear first
##The species are named in a simiar fashion where the species from the most recent nodes are given the lowest number names

ms<-getMS(tree,1,"T",new.names = T) ## We've got an ms string in [[1]] and the new names in [[2]]
tree$tip.label<-ms[[2]][,2]
tree$node.label<-as.character(tree$Nnode:1)
ms<-ms[[1]] ## for ease we dn't care about the names anymore

##Lets look to see if the renaming worked
plot(tree) ##Great!
nodelabels(text=as.character(c(9,7,8,4,6,5,3,2,1))) ##Now I could probably write something that does this but we're using this just as a test so renaming them by hand will suffice for now

##Now we've got a string in the ms format
ms

##Since we got the string in the ms format, let's first reformat it so it is more workable 
## this returns a list, the second list item is the real workhorse that contains all tree information
##The first item is just the stuff at the start of the ms string
ms_frame<-get_ms_frame(ms)[[2]]
ms_frame

##The next thing we would do in addHybrid is slice the clade we want to add hybrid events
##However, this requires a a node to slice by. In the pipeline we find a clade by using find_clade()
##Let's say that we want the smallest clade that has at least 40% of the species (4 species),
##Looking at the plot, we can easily see that we want node  8 - the node that contains the four species 1,2,7, and 10.
##Let's try find_clade() to see if we get node 8
node<-find_clade(ms_frame,0.4)
node ## great!

##The next thing to do is slice the tree to just give us nodes in that clade 
##We should reduce our data frame to nodes 1,4, and 8
clade_ms<-clade_slice(ms_frame,node) ##Great!
clade_ms

##We can see that only species 1,2,7 and 10 are listed in the data frame


##Since we say events can happen anywhere throughout the clade we want from 0 to the node that defines the clade
start<-0
end<-max(clade_ms$time)

##Now that we've got our clade we can run through the process of adding an event

##First, we pick a time 
time<-runif(1,start,end) ## if the seed is still there we should get 0.1498328

##For the purposes of weighting event probabilities properly we need to slice the tree at that time.
##We are saying that things after this time haven't happened yet so we shouldn't factor them into our weightings
##Node 1 occurs after this event so we shouldn't see this node after we slice
time_ms<-timeslice(clade_ms,time)
time_ms

##You should also notice that the times of the remaining nodes got reduced by the amount of time in the event. This puts the branch lengths to what they would be at the time of the event
##Also, even though species 2 doesn't really exist yet at this time point, the branch that  would have led up to node 1 is now treated as the branch species 2. It's an ms thing that I'm not going to get into right now.

##Now we want to find the descendents at each node
desc<-find_progeny(time_ms)

##although the output here refers [[1]] and [[2]], we can really think of these as the descendents for nodes 4 and 8. The indices line up with nodes 4 and 8 in time_ms.
##Looks pretty good

##We then use this to find out which species are introduced at a node. I don't have this step in its own function but you can try to trust that new_join() in addHybrid.R does this.
##Even looking at desc, it should be easy to see that species to see the new species for each node:
##Node 4: 2,7
##Node 8: 10

##We then use this information to weight each event
##Events  between two species are wieghted by the inverse of the branch length from one species to the MRCA of the other species.
## They are then wieghted by how many "species" are shared among the common ancestor
##We should see that here:
joins<-new_join(time_ms,desc)
joins

##to try and explain that second point better, let's look at events between 2,7 and 2,10.
## They both have a branch length of 1.238 which should give a weight of about .81. 
##However, we divide both of those by two since there are two species at that time. Effectively we are chosing the branch by using the branch length of 1.63 and then saying that picking either 2 or 7 is equally likely.

##All that's left to do now is select a species pair based off of the weights and randomly order the species pair
##The random ordering of the species pair is used to say that the directionality of the event can go either way
pair<-joins[sample(seq(nrow(joins)),1,replace=TRUE,prob=joins$weight),]
pair<-pair<-sample(pair[c(1,2)],2)
pair

##Cool! We've successfully picked an introgression event. You can try this with your own tree if you want!


##We can imagine to do this however many times we want, recording the species pairs and times while doing so to generate a list of hybrid events

##All we need to do then is get this record of events into the ms string
##I'm not going to do that here but you can see how it is done in lines 78-90 in addHybrids.R
##NOTE: if you're going to look at lines 78-90 you may also want to look at 28-35 too. These lines are the vectors that get populated when doing a bunch of events and go into creating data frames that eventually merge into the ms string