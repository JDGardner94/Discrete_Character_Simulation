# Discrete correlation script
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/24/2017

library(ape)
library(corHMM)
library(geiger)
library(MCMCglmm)
library(phytools)

NSIMS <- 1000 # Number of iterations in simulation
pvec <- numeric(NSIMS) # vector for pvalues
nodes.meanprob <- numeric(NSIMS) # mean of highest probability for all nodes for each tree
nodes.varprob <- numeric(NSIMS) # variance of highest probability for all nodes in each iteration
trees <- character(NSIMS) # store trees
treesize <- numeric(NSIMS) # vector for number of taxa each iteration
character.data <- table(NSIMS) # Table for all data
DarwinCharacters <- table(NSIMS) # Table for Darwin characters
UnrepBurstCharacters <- table(NSIMS) # Table for UnrepBurst characters
UnrepBurstCharacters.mean <- numeric(NSIMS) # Table for frequency of 1's in UnrepBurst data
cladesize <- numeric(NSIMS) # size of clade under study
Data.BothDarwinCharacters <- table(NSIMS) # Table for Darwin characters
Data.Darwin.UnrepBurstCharacters <- table(NSIMS) # Table for UnrepBurst characters
model.BothDarwin.lrt <- numeric(NSIMS) # vector for Darwin's scenario likelihood ratio test statistic
model.Darwin.UnrepBurst.lrt <- numeric(NSIMS) # vector for unreplicated burst scenario likelihood ratio test statistic

# Maddison and FitzJohn Replication (Pagel's Method)
model.BothDarwin.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.Darwin.UnrepBurst.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
UnrepBurstCharacters.mean <- numeric(NSIMS) # vector for frequency of 1's in UnrepBurst data


i<-0
while (i<NSIMS){
    i<-i+1
cat("iteration", i, '\n')

NTAXA <- sample(50:1000, 1) # number of taxa
NODES <- NTAXA-1 #number of nodes in tree
node.prob <- numeric(NODES) # vector for highest probability nodes in each tree

## Simulate tree in Geiger
phy <- sim.bdtree(b=0.1, d=0, stop="taxa", n=(NTAXA))

## Generate Data
# Sets all data to 0
darwindata <-rep(0, NTAXA)

#Find clade to change values from 0 to 1
cladetip <- tips(phy, mrca(phy)["s40", "s20"])
cladetip <- substr(cladetip,start=2,stop=4) #removes "s" from tip name
tips <- as.numeric(cladetip)

# Keeps the clade size constrained to 40 taxa
	while (length(tips) > 40 | length(tips) < 40){
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
unrepburstdata <-rep(0, NTAXA)
unrepbursttips <- sample(tips, size=length(tips)/2, replace = FALSE)
unrepburstdata[c(unrepbursttips)]<-1

Data.BothDarwinCharacters <- factor(paste(as.factor(darwindata),as.factor(darwindata),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurstCharacters <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata),sep=","),levels=c("0,0","0,1","1,0","1,1"))


###Modelling of discrete correlated evolution.


# ##Replicating Maddison and FitzJohn (2015) - Pagel (1994) Method
iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE) # independent matrix
dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE) # dependent matrix

# #Models use try command to detect and trap error message. If error message, 'next' halts the processing of the current iteration and advances the looping index

#Darwin's Scenario
model.BothDarwin.iQ <- try(ace(Data.BothDarwinCharacters, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.BothDarwin.iQ), "try-error"))) {
	writeLines("model.BothDarwin.iQ failed")
	i<-i-1
	next
	} else {
	}

model.BothDarwin.dQ <- try(ace(Data.BothDarwinCharacters, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.BothDarwin.dQ), "try-error"))) {
	writeLines("model.BothDarwin.dQ failed")
	i<-i-1
	next
	} else {
	}

#Unreplicated Burst Scenario
model.Darwin.UnrepBurst.iQ  <- try(ace(Data.Darwin.UnrepBurstCharacters, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst.iQ ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst.iQ failed")
	i<-i-1
	next
	} else {
	}

model.Darwin.UnrepBurst.dQ <- try(ace(Data.Darwin.UnrepBurstCharacters, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst.dQ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst.dQ failed")
	i<-i-1
	next
	} else {
	}

#Likelihood Ratio Tests
if(model.BothDarwin.dQ$loglik > model.BothDarwin.iQ$loglik){
model.BothDarwin.lrt[i] <- 2*(model.BothDarwin.dQ$loglik - model.BothDarwin.iQ$loglik) #low p-values favor dependent model
	} else {
	model.BothDarwin.lrt[i] <- 2*(model.BothDarwin.iQ$loglik - model.BothDarwin.dQ$loglik)
	}

if(model.Darwin.UnrepBurst.dQ$loglik > model.Darwin.UnrepBurst.iQ$loglik){
model.Darwin.UnrepBurst.lrt[i] <- 2*(model.Darwin.UnrepBurst.dQ$loglik - model.Darwin.UnrepBurst.iQ$loglik) #low p-values favor dependent model
	} else {
	model.Darwin.UnrepBurst.lrt[i] <- 2*(model.Darwin.UnrepBurst.iQ$loglik - model.Darwin.UnrepBurst.dQ$loglik)
	}

##Saves data in a table for each iteration. For explanations of data structure see Maddison and Fitzjohn. Syst. Biol. 64(1):127�136, 2015.
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
treesize[i] <- NTAXA
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)
model.BothDarwin.pvec[i] <- 1 - pchisq(model.BothDarwin.lrt[i], 4)
model.Darwin.UnrepBurst.pvec[i] <- 1 - pchisq(model.Darwin.UnrepBurst.lrt[i], 4)

}

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##General Results
length(trees)
length(cladesize)
cladesize
length(treesize)
hist(treesize)
treesize
length(DarwinCharacters)
length(UnrepBurstCharacters)
length(UnrepBurstCharacters.mean)
length(model.BothDarwin.pvec)
length(model.Darwin.UnrepBurst.pvec)
sum(model.BothDarwin.pvec <= 0.05)
sum(model.Darwin.UnrepBurst.pvec <= 0.05)

##Pagel Darwin Scenario Results
model.BothDarwin.lrt
model.BothDarwin.pvec
hist(model.BothDarwin.pvec)
min(model.BothDarwin.pvec)
max(model.BothDarwin.pvec)
median(model.BothDarwin.pvec)
plot(cladesize, model.BothDarwin.pvec)
PagelDarwin.lm <- lm(model.BothDarwin.pvec ~ cladesize)
summary(PagelDarwin.lm)
plot(treesize, model.BothDarwin.lrt)
PagelDarwin.lm2 <- lm(model.BothDarwin.lrt ~ treesize)
summary(PagelDarwin.lm2)
plot(treesize, model.BothDarwin.pvec)
PagelDarwin.lm3 <- lm(model.BothDarwin.pvec ~ treesize)
summary(PagelDarwin.lm3)

##Pagel Unreplicated Burst Scenario Results
model.Darwin.UnrepBurst.lrt
model.Darwin.UnrepBurst.pvec
hist(model.Darwin.UnrepBurst.pvec)
min(model.Darwin.UnrepBurst.pvec)
max(model.Darwin.UnrepBurst.pvec)
median(model.Darwin.UnrepBurst.pvec)
plot(cladesize, model.Darwin.UnrepBurst.pvec)
PagelUnrepburst.lm <- lm(model.Darwin.UnrepBurst.pvec ~ cladesize)
summary(PagelUnrepburst.lm)
plot(treesize, model.Darwin.UnrepBurst.lrt)
PagelUnrepburst.lm2 <- lm(model.Darwin.UnrepBurst.lrt ~ treesize)
summary(PagelUnrepburst.lm2)
plot(treesize, model.Darwin.UnrepBurst.pvec)
PagelUnrepburst.lm3 <- lm(model.Darwin.UnrepBurst.pvec ~ treesize)
summary(PagelUnrepburst.lm3)
