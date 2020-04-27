# Discrete correlation script
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/24/2017

start.time <- Sys.time()

library(ape)
library(corHMM)
library(geiger)
library(MCMCglmm)
library(phytools)

NSIMS <- 100 # Number of iterations in simulation
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

# Threshold Model - Darwin Scenario
mcmcdarwin.r <- numeric(NSIMS) # vector for R values - Darwin Scenario Threshold Model
mcmcdarwin.r2 <- numeric(NSIMS) # vector for R^2 values - Darwin Scenario Threshold Model

# Threshold Model - Unreplicated Burst Scenario
mcmcunrepburst.r <- numeric(NSIMS) # vector for R values - UnrepBurst Scenario Threshold Model
mcmcunrepburst.r2 <- numeric(NSIMS) # vector for R^2 values - UnrepBurst Scenario Threshold Model


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

Data.BothDarwinCharacters <- factor(paste(as.factor(darwindata),as.factor(darwindata),sep=","),levels=c("0,0","0,1","1,0","1,1"))
Data.Darwin.UnrepBurstCharacters <- factor(paste(as.factor(darwindata),as.factor(unrepburstdata),sep=","),levels=c("0,0","0,1","1,0","1,1"))


###Modelling of discrete correlated evolution.

# ##Threshold Model in Phytools (Revell, 2012)
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

	#ThreshBayes Darwin Scenario Analysis
mcmc.darwin <- threshBayes(phy,darwindata.MatrixBind,types=c("disc", "disc"),ngen=ngen, control = list(sample = sample))
burnin.darwin <- colMeans(mcmc.darwin$par[10:nrow(mcmc.darwin$par),2:6])

	#ThreshBayes Unreplicated Burst Scenario Analysis
mcmc.unrepburst <- threshBayes(phy,UnrepBurst.MatrixBind,types=c("disc", "disc"),ngen=ngen, control = list(sample = sample))
burnin.unrepburst <- colMeans(mcmc.unrepburst$par[10:nrow(mcmc.unrepburst$par),2:6])

	#Saves data for the Threshold models
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)

	#Threshold Darwin Scenario P-Values
		#Extract R and T-test Approach
mcmcdarwin.r[i] <-mcmc.darwin$par[2501:nrow(mcmc.darwin$par),"r"]
mcmcdarwin.r2[i] <- (mcmcdarwin.r)^2

	#Threshold Unreplicated Burst Scenario P-Values
		#Extract R and T-test Approach
mcmcunrepburst.r[i] <-mcmc.unrepburst$par[2501:nrow(mcmc.unrepburst$par),"r"]
mcmcunrepburst.r2[i] <- (mcmcunrepburst.r)^2

}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##General Results
length(trees)
length(cladesize)
length(DarwinCharacters)
length(UnrepBurstCharacters)
length(UnrepBurstCharacters.mean)

##Threshold Darwin Scenario Results
mcmcdarwin.r
hist(mcmcdarwin.r)
boxplot(mcmcdarwin.r)
min(mcmcdarwin.r)
max(mcmcdarwin.r)
median(mcmcdarwin.r)
quantile(mcmcdarwin.r, c(0.025, 0.975)) #95% credible interval of posterior distribution of r
plot(cladesize, mcmcdarwin.r)
ThresholdDarwin.lm <- lm(mcmcdarwin.r ~ cladesize)
summary(ThresholdDarwin.lm)
plot(mcmc.darwin$par[, "logL"], type = "l", xlab = "generation", ylab = "logL")

mcmcdarwin.r^2
hist(mcmcdarwin.r^2)
boxplot(mcmcdarwin.r^2)
min(mcmcdarwin.r^2)
max(mcmcdarwin.r^2)
median(mcmcdarwin.r^2)
quantile(mcmcdarwin.r^2, c(0.025, 0.975)) #95% credible interval of posterior distribution of r

##Threshold Unreplicated Burst Scenario Results
mcmcunrepburst.r
hist(mcmcunrepburst.r)
boxplot(mcmcunrepburst.r)
min(mcmcunrepburst.r)
max(mcmcunrepburst.r)
median(mcmcunrepburst.r)
quantile(mcmcunrepburst.r, c(0.025, 0.975)) #95% credible interval of posterior distribution of r
plot(cladesize, mcmcunrepburst.r)
ThresholdUnrepburst.lm <- lm(mcmcunrepburst.r ~ cladesize)
summary(ThresholdUnrepburst.lm)
plot(mcmc.unrepburst$par[, "logL"], type = "l", xlab = "generation", ylab = "logL")

mcmcunrepburst.r^2
hist(mcmcunrepburst.r^2)
boxplot(mcmcunrepburst.r^2)
min(mcmcunrepburst.r^2)
max(mcmcunrepburst.r^2)
median(mcmcunrepburst.r^2)
quantile(mcmcunrepburst.r^2, c(0.025, 0.975)) #95% credible interval of posterior distribution of r
