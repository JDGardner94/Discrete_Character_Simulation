# Discrete correlation script - Control Tests
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/18/2017

library("ape")
library("corHMM")
library("geiger")
library("MCMCglmm")
library("phytools")

NSIMS <- 1000 # Number of iterations in simulation
NTAXA <- 100 # number of taxa
NODES <- NTAXA-1 #number of nodes in tree
pvec <- numeric(NSIMS) # vector for pvalues
node.prob <- numeric(NODES) # vector for highest probability nodes in each tree
nodes.meanprob <- numeric(NSIMS) # mean of highest probability for all nodes for each tree
nodes.varprob <- numeric(NSIMS) # variance of highest probability for all nodes in each iteration
trees <- character(NSIMS) # store trees
character.data <- table(NSIMS) # Table for all data
DarwinCharacters <- table(NSIMS) # Table for Darwin characters
UnrepBurstCharacters <- table(NSIMS) # Table for UnrepBurst characters
UnrepBurstCharacters.mean <- numeric(NSIMS) # Table for frequency of 1's in UnrepBurst data
Data.BothDarwinCharacters <- table(NSIMS) # Table for Binded Darwin characters
Data.UnrepBurst.UnrepBurst2Characters <- table(NSIMS) # Table for Binded Darwin and UnrepBurst characters
Data.BothUnrepBurstCharacters <- table(NSIMS) # Table for Positive Control characters
Data.BothUnrepBurstCharacters.mean <- numeric(NSIMS) # Table for frequency of 1's in Positive Control data
Data.UnrepBurst.UnrepBurst2 <- table(NSIMS) # Table for Negative Control characters
Data.UnrepBurst.UnrepBurst2.mean <- numeric(NSIMS) # Table for frequency of 1's in Negative Control data
cladesize <- numeric(NSIMS) # size of clade under study

# Pagel (1994) Method
model.BothUnrepBurst.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.UnrepBurst.UnrepBurst2.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication


i<-0
while (i<NSIMS){
    i<-i+1
cat("iteration", i, '\n')

## Simulate tree in Geiger
phy <- sim.bdtree(b=0.1, d=0, stop="taxa", n=(NTAXA))

## Generate Data
# Sets all data to 0
darwindata <-rep(0, NTAXA)

#Find clade to change values from 0 to 1
cladetip <- tips(phy, mrca(phy)["s40", "s20"])
cladetip <- substr(cladetip,start=2,stop=4) #removes "s" from tip name
tips <- as.numeric(cladetip)

# Keeps the clade size constrained between 40 and 60 taxa
	while (length(tips) > 60 | length(tips) < 40){
	phy <- sim.bdtree(b=0.1, d=0, stop="taxa", n=(NTAXA))

	# Sets all data to 0
	darwindata <-rep(0, NTAXA)

	#Find clade to change values from 0 to 1
	cladetip <- tips(phy, mrca(phy)["s40", "s20"])
	cladetip <- substr(cladetip,start=2,stop=4) #removes "s" from tip name
	tips <- as.numeric(cladetip)
	}

# Sets value for clade to 1 to make Darwin data
darwindata[c(tips)]<-1

# Randomly change 1s back to zeros to make Unreplicated Burst data
unrepburstdata <- rep(0, NTAXA)
unrepbursttips <- sample(NTAXA, size=length(tips), replace = FALSE)
unrepburstdata[c(unrepbursttips)]<-1

# This is for our negative control: UnrepBurst vs UnrepBurst2 should have high p-val
unrepburstdata2 <-rep(0, NTAXA)
unrepbursttips2 <- sample(NTAXA, size=length(tips), replace = FALSE)
unrepburstdata2[c(unrepbursttips2)]<-1

Data.BothDarwinCharacters <- factor(paste(as.factor(darwindata),as.factor(darwindata),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurstCharacters <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.BothUnrepBurstCharacters <- factor(paste(as.factor(unrepburstdata),as.factor(unrepburstdata),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.UnrepBurst.UnrepBurst2 <- factor(paste(as.factor(unrepburstdata),as.factor(unrepburstdata2),sep=","),levels=c("0,0","0,1","1,0","1,1"))

# ##Modelling of discrete correlated evolution.

# ##Replicating Maddison and FitzJohn (2015) - Pagel (1994) Method
iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE) # independent matrix
dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE) # dependent matrix

# #Models use try command to detect and trap error message. If error message, 'next' halts the processing of the current iteration and advances the looping index

#Positive Control: Should have low p-val
model.BothUnrepBurst.iQ <- try(ace(Data.BothUnrepBurstCharacters, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.BothUnrepBurst.iQ), "try-error"))) {
	writeLines("model.BothUnrepBurst.iQ failed")
	i<-i-1
	next
	} else {
	}

model.BothUnrepBurst.dQ <- try(ace(Data.BothUnrepBurstCharacters, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.BothUnrepBurst.dQ), "try-error"))) {
	writeLines("model.BothUnrepBurst.dQ failed")
	i<-i-1
	next
	} else {
	}

#Negative Control: Should have high p-val
model.UnrepBurst.UnrepBurst2.iQ <- try(ace(Data.UnrepBurst.UnrepBurst2, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.UnrepBurst.UnrepBurst2.iQ), "try-error"))) {
	writeLines("model.UnrepBurst.UnrepBurst2.iQ failed")
	i<-i-1
	next
	} else {
	}

model.UnrepBurst.UnrepBurst2.dQ <- try(ace(Data.UnrepBurst.UnrepBurst2, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.UnrepBurst.UnrepBurst2.dQ), "try-error"))) {
	writeLines("model.UnrepBurst.UnrepBurst2.dQ failed")
	i<-i-1
	next
	} else {
	}

#Likelihood Ratio Tests
if(model.BothUnrepBurst.dQ$loglik > model.BothUnrepBurst.iQ$loglik){	
model.BothUnrepBurst.lrt <- 2*(model.BothUnrepBurst.dQ$loglik - model.BothUnrepBurst.iQ$loglik) #low p-values favor dependent model
	} else {
	model.BothUnrepBurst.lrt <- 2*(model.BothUnrepBurst.iQ$loglik - model.BothUnrepBurst.dQ$loglik)
	} 

if(model.UnrepBurst.UnrepBurst2.dQ$loglik > model.UnrepBurst.UnrepBurst2.iQ$loglik) {
model.UnrepBurst.UnrepBurst2.lrt <- 2*(model.UnrepBurst.UnrepBurst2.dQ$loglik - model.UnrepBurst.UnrepBurst2.iQ$loglik) #low p-values favor dependent model
	} else {
	model.UnrepBurst.UnrepBurst2.lrt <- 2*(model.UnrepBurst.UnrepBurst2.iQ$loglik - model.UnrepBurst.UnrepBurst2.dQ$loglik)
	} 

#Saves data in a table for each iteration. For explanations of data structure see Maddison and Fitzjohn. Syst. Biol. 64(1):127–136, 2015.
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)
model.BothUnrepBurst.pvec[i] <- 1 - pchisq(model.BothUnrepBurst.lrt, 4)
model.UnrepBurst.UnrepBurst2.pvec[i] <- 1 - pchisq(model.UnrepBurst.UnrepBurst2.lrt, 4)


}

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##Maddison Results
length(trees)
length(cladesize)
length(DarwinCharacters)
length(UnrepBurstCharacters)
length(UnrepBurstCharacters.mean)
model.BothUnrepBurst.pvec #list of all p-values
hist(model.BothUnrepBurst.pvec) #positive control histogram
model.UnrepBurst.UnrepBurst2.pvec #list of all p-values
hist(model.UnrepBurst.UnrepBurst2.pvec) #negative control histogram
hist(UnrepBurstCharacters.mean)




