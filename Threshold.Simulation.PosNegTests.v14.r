# Discrete correlation script - Control Tests
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/24/2017

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

# Threshold Model - Positive Scenario
mcmc.bothunrepburst.R <- numeric(NSIMS) # vector for R values - Darwin Scenario Threshold Model
mcmc.bothunrepburst.R2 <- numeric(NSIMS) # vector for R^2 values - Darwin Scenario Threshold Model

# Threshold Model - Negative Scenario
mcmc.UnrepBurst.UnrepBurst2.R <- numeric(NSIMS) # vector for R values - UnrepBurst Scenario Threshold Model
mcmc.UnrepBurst.UnrepBurst2.R2 <- numeric(NSIMS) # vector for R^2 values - UnrepBurst Scenario Threshold Model


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

# ##Threshold Model in Phytools
ngen <- 10000000
sample <- 1000
burnin <- ngen*0.25
	#ThreshBayes Darwin Scenario Data
darwindata.Matrix <- matrix(c(darwindata))
darwindata.MatrixBind <- cbind(darwindata.Matrix,darwindata.Matrix)
dimnames(darwindata.MatrixBind) <- list(as.list(phy$tip.label))

	#ThreshBayes Unreplicated Burst Scenario Data
UnrepBurst.Matrix <- matrix(c(unrepburstdata))
UnrepBurst.MatrixBind <- cbind(darwindata.Matrix,UnrepBurst.Matrix)
dimnames(UnrepBurst.MatrixBind) <- list(as.list(phy$tip.label))

	#ThreshBayes Positive Control Data
BothUnrepBurst.MatrixBind <- cbind(UnrepBurst.Matrix,UnrepBurst.Matrix)
dimnames(BothUnrepBurst.MatrixBind) <- list(as.list(phy$tip.label))

	#ThreshBayes Negative Control Data
UnrepBurst2.Matrix <- matrix(c(unrepburstdata2))
UnrepBurst.UnrepBurst2.MatrixBind <- cbind(UnrepBurst.Matrix,UnrepBurst2.Matrix)
dimnames(UnrepBurst.UnrepBurst2.MatrixBind) <- list(as.list(phy$tip.label))

	#ThreshBayes Positive Control Analysis
mcmc.bothunrepburst <- threshBayes(phy,BothUnrepBurst.MatrixBind,types=c("discrete"),ngen=ngen)
burnin.bothunrepburst <-colMeans(mcmc.bothunrepburst$par[10:nrow(mcmc.bothunrepburst$par),2:6])

	#ThreshBayes Negative Control Analysis
mcmc.UnrepBurst.UnrepBurst2 <- threshBayes(phy,UnrepBurst.UnrepBurst2.MatrixBind,types=c("discrete"),ngen=ngen)
burnin.UnrepBurst.UnrepBurst2 <-colMeans(mcmc.UnrepBurst.UnrepBurst2$par[10:nrow(mcmc.UnrepBurst.UnrepBurst2$par),2:6])

	#Saves data for the Threshold models
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)

	#Threshold Positive Control R Values 
mcmc.bothunrepburst.R[i] <- mcmc.bothunrepburst$par[2501:nrow(mcmc.bothunrepburst$par),"r"] #extract R value

	#Threshold Negative Control R Values
mcmc.UnrepBurst.UnrepBurst2.R[i] <- mcmc.UnrepBurst.UnrepBurst2$par[2501:nrow(mcmc.UnrepBurst.UnrepBurst2$par),"r"] #extract R value

}

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results

##Threshold Results
mcmc.bothunrepburst.R
hist(mcmc.bothunrepburst.R)
median(mcmc.bothunrepburst.R) #median r value
quantile(mcmc.bothunrepburst.R, c(0.025, 0.975)) #95% credible interval of posterior distribution of r
boxplot(mcmc.bothunrepburst.R) #box plot of posterior distribution of r
plot(mcmc.bothunrepburst$par[, "logL"], type = "l", xlab = "generation", ylab = "logL")

mcmc.UnrepBurst.UnrepBurst2.R
hist(mcmc.UnrepBurst.UnrepBurst2.R)
median(mcmc.UnrepBurst.UnrepBurst2.R) #median r value
quantile(mcmc.UnrepBurst.UnrepBurst2.R, c(0.025, 0.975)) #95% credible interval of posterior distribution of r
boxplot(mcmc.UnrepBurst.UnrepBurst2.R) #box plot of posterior distribution of r
plot(mcmc.UnrepBurst.UnrepBurst2$par[, "logL"], type = "l", xlab = "generation", ylab = "logL")

