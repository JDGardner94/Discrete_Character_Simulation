# Discrete correlation script
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/18/2017

library(ape)
library(corHMM)
library(geiger)
library(MCMCglmm)
library(phytools)

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
cladesize <- numeric(NSIMS) # size of clade under study
Data.BothDarwinCharacters <- table(NSIMS) # Table for Darwin characters
Data.Darwin.UnrepBurstCharacters <- table(NSIMS) # Table for UnrepBurst characters
Data.Darwin.UnrepBurst2 <- table(NSIMS) # Table for UnrepBurst characters 2
Data.Darwin.UnrepBurst3 <- table(NSIMS) # Table for UnrepBurst characters 3
Data.Darwin.UnrepBurst4 <- table(NSIMS) # Table for UnrepBurst characters 4
Data.Darwin.UnrepBurst5 <- table(NSIMS) # Table for UnrepBurst characters 5

# Maddison and FitzJohn Replication (Pagel's Method)
model.BothDarwin.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.Darwin.UnrepBurst.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
UnrepBurstCharacters.mean <- numeric(NSIMS) # vector for frequency of 1's in UnrepBurst data
model.Darwin.UnrepBurst2.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.Darwin.UnrepBurst3.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.Darwin.UnrepBurst4.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.Darwin.UnrepBurst5.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication


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
unrepburstdata <-rep(0, NTAXA)
unrepbursttips <- sample(tips, size=length(tips)/2, replace = FALSE)
unrepburstdata[c(unrepbursttips)]<-1

# Progressive departure from Darwin to Unreplicated Burst (predict progressive decrease in p-value)
unrepburstdata2 <- darwindata
unrepburstsample2 <- sample(tips, size=length(tips)/8, replace = FALSE)
unrepburstdata2[c(unrepburstsample2)] <-0 

unrepburstdata3 <- darwindata
unrepburstsample3 <- sample(tips, size=length(tips)/6, replace = FALSE)
unrepburstdata3[c(unrepburstsample3)] <-0

unrepburstdata4 <- darwindata
unrepburstsample4 <- sample(tips, size=length(tips)/4, replace = FALSE)
unrepburstdata4[c(unrepburstsample4)] <-0 

unrepburstdata5 <- darwindata
unrepburstsample5 <- sample(tips, size=length(tips)/2, replace = FALSE)
unrepburstdata5[c(unrepburstsample5)] <-0 

Data.BothDarwinCharacters <- factor(paste(as.factor(darwindata),as.factor(darwindata),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurstCharacters <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata),sep=","),levels=c("0,0","0,1","1,0","1,1"))

Data.Darwin.UnrepBurst2 <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata2),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurst3 <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata3),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurst4 <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata4),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurst5 <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata5),sep=","),levels=c("0,0","0,1","1,0","1,1"))


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

#Unrep Burst 2 (divide tips by 6)
model.Darwin.UnrepBurst2.iQ  <- try(ace(Data.Darwin.UnrepBurst2, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst2.iQ ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst2.iQ failed")
	i<-i-1
	next
	} else {
	}

model.Darwin.UnrepBurst2.dQ <- try(ace(Data.Darwin.UnrepBurst2, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst2.dQ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst2.dQ failed")
	i<-i-1
	next
	} else {
	}

#Unrep Burst 3 (divide tips by 5)
model.Darwin.UnrepBurst3.iQ  <- try(ace(Data.Darwin.UnrepBurst3, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst3.iQ ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst3.iQ failed")
	i<-i-1
	next
	} else {
	}

model.Darwin.UnrepBurst3.dQ <- try(ace(Data.Darwin.UnrepBurst3, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst3.dQ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst3.dQ failed")
	i<-i-1
	next
	} else {
	}

#Unrep Burst 4 (divide tips by 4)
	model.Darwin.UnrepBurst4.iQ  <- try(ace(Data.Darwin.UnrepBurst4, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst4.iQ ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst4.iQ failed")
	i<-i-1
	next
	} else {
	}

model.Darwin.UnrepBurst4.dQ <- try(ace(Data.Darwin.UnrepBurst4, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst4.dQ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst4.dQ failed")
	i<-i-1
	next
	} else {
	}

#Unrep Burst 5 (divide by 3)
	model.Darwin.UnrepBurst5.iQ  <- try(ace(Data.Darwin.UnrepBurst5, phy, type = 'd', model = iQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst5.iQ ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst5.iQ failed")
	i<-i-1
	next
	} else {
	}

model.Darwin.UnrepBurst5.dQ <- try(ace(Data.Darwin.UnrepBurst5, phy, type = 'd', model = dQ, use.expm = FALSE, use.eigen = TRUE))
if(isTRUE(all.equal(class(model.Darwin.UnrepBurst5.dQ), "try-error"))) {
	writeLines("model.Darwin.UnrepBurst5.dQ failed")
	i<-i-1
	next
	} else {
	}

#Likelihood Ratio Tests
if(model.BothDarwin.dQ$loglik > model.BothDarwin.iQ$loglik){	
model.BothDarwin.lrt <- 2*(model.BothDarwin.dQ$loglik - model.BothDarwin.iQ$loglik) #low p-values favor dependent model
	} else {
	model.BothUnrepBurst.lrt <- 2*(model.BothDarwin.iQ$loglik - model.BothDarwin.dQ$loglik)
	} 	

if(model.Darwin.UnrepBurst.dQ$loglik > model.Darwin.UnrepBurst.iQ$loglik){	
model.Darwin.UnrepBurst.lrt <- 2*(model.Darwin.UnrepBurst.dQ$loglik - model.Darwin.UnrepBurst.iQ$loglik) #low p-values favor dependent model
	} else {
	model.Darwin.UnrepBurst.lrt <- 2*(model.Darwin.UnrepBurst.iQ$loglik - model.Darwin.UnrepBurst.dQ$loglik)
	} 
		
if(model.Darwin.UnrepBurst2.dQ$loglik > model.Darwin.UnrepBurst2.iQ$loglik){	
model.Darwin.UnrepBurst2.lrt <- 2*(model.Darwin.UnrepBurst2.dQ$loglik - model.Darwin.UnrepBurst2.iQ$loglik) #low p-values favor dependent model
	} else {
	model.Darwin.UnrepBurst2.lrt <- 2*(model.Darwin.UnrepBurst2.iQ$loglik - model.Darwin.UnrepBurst2.dQ$loglik)
	} 

if(model.Darwin.UnrepBurst3.dQ$loglik > model.Darwin.UnrepBurst3.iQ$loglik){	
model.Darwin.UnrepBurst3.lrt <- 2*(model.Darwin.UnrepBurst3.dQ$loglik - model.Darwin.UnrepBurst3.iQ$loglik) #low p-values favor dependent model
	} else {
	model.Darwin.UnrepBurst3.lrt <- 2*(model.Darwin.UnrepBurst3.iQ$loglik - model.Darwin.UnrepBurst3.dQ$loglik)
	} 

if(model.Darwin.UnrepBurst4.dQ$loglik > model.Darwin.UnrepBurst4.iQ$loglik){	
model.Darwin.UnrepBurst4.lrt <- 2*(model.Darwin.UnrepBurst4.dQ$loglik - model.Darwin.UnrepBurst4.iQ$loglik) #low p-values favor dependent model
	} else {
	model.Darwin.UnrepBurst4.lrt <- 2*(model.Darwin.UnrepBurst4.iQ$loglik - model.Darwin.UnrepBurst4.dQ$loglik)
	} 
	
if(model.Darwin.UnrepBurst5.dQ$loglik > model.Darwin.UnrepBurst5.iQ$loglik){	
model.Darwin.UnrepBurst5.lrt <- 2*(model.Darwin.UnrepBurst5.dQ$loglik - model.Darwin.UnrepBurst5.iQ$loglik) #low p-values favor dependent model
	} else {
	model.Darwin.UnrepBurst5.lrt <- 2*(model.Darwin.UnrepBurst5.iQ$loglik - model.Darwin.UnrepBurst5.dQ$loglik)
	} 	

##Saves data in a table for each iteration. For explanations of data structure see Maddison and Fitzjohn. Syst. Biol. 64(1):127136, 2015.
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)
model.BothDarwin.pvec[i] <- 1 - pchisq(model.BothDarwin.lrt, 4)
model.Darwin.UnrepBurst.pvec[i] <- 1 - pchisq(model.Darwin.UnrepBurst.lrt, 4)
model.Darwin.UnrepBurst2.pvec[i] <- 1 - pchisq(model.Darwin.UnrepBurst2.lrt, 4)
model.Darwin.UnrepBurst3.pvec[i] <- 1 - pchisq(model.Darwin.UnrepBurst3.lrt, 4)
model.Darwin.UnrepBurst4.pvec[i] <- 1 - pchisq(model.Darwin.UnrepBurst4.lrt, 4)
model.Darwin.UnrepBurst5.pvec[i] <- 1 - pchisq(model.Darwin.UnrepBurst5.lrt, 4)


}

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##General Results
length(trees)
length(cladesize)
length(DarwinCharacters)
length(UnrepBurstCharacters)
length(UnrepBurstCharacters.mean)
length(model.BothDarwin.pvec)
length(model.Darwin.UnrepBurst.pvec)
sum(model.BothDarwin.pvec <= 0.05)
sum(model.Darwin.UnrepBurst.pvec <= 0.05)

##Pagel Darwin Scenario Results
model.BothDarwin.pvec
hist(model.BothDarwin.pvec)
min(model.BothDarwin.pvec)
max(model.BothDarwin.pvec)	
median(model.BothDarwin.pvec)
plot(cladesize, model.BothDarwin.pvec)
PagelDarwin.lm <- lm(model.BothDarwin.pvec ~ cladesize)
summary(PagelDarwin.lm)

##Pagel Unreplicated Burst Scenario Results
model.Darwin.UnrepBurst.pvec
hist(model.Darwin.UnrepBurst.pvec)
min(model.Darwin.UnrepBurst.pvec)
max(model.Darwin.UnrepBurst.pvec)	
median(model.Darwin.UnrepBurst.pvec)
plot(cladesize, model.Darwin.UnrepBurst.pvec)
PagelUnrepburst.lm <- lm(model.Darwin.UnrepBurst.pvec ~ cladesize)
summary(PagelUnrepburst.lm)

##Pagel Progressive Departure Results
model.Darwin.UnrepBurst2.pvec
hist(model.Darwin.UnrepBurst2.pvec)
min(model.Darwin.UnrepBurst2.pvec)
max(model.Darwin.UnrepBurst2.pvec)
median(model.Darwin.UnrepBurst2.pvec)
model.Darwin.UnrepBurst3.pvec
hist(model.Darwin.UnrepBurst3.pvec)
min(model.Darwin.UnrepBurst3.pvec)
max(model.Darwin.UnrepBurst3.pvec)
median(model.Darwin.UnrepBurst3.pvec)
model.Darwin.UnrepBurst4.pvec
hist(model.Darwin.UnrepBurst4.pvec)
min(model.Darwin.UnrepBurst4.pvec)
max(model.Darwin.UnrepBurst4.pvec)
median(model.Darwin.UnrepBurst4.pvec)
model.Darwin.UnrepBurst5.pvec
hist(model.Darwin.UnrepBurst5.pvec)
min(model.Darwin.UnrepBurst5.pvec)
max(model.Darwin.UnrepBurst5.pvec)
median(model.Darwin.UnrepBurst5.pvec)

