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

# Pagel (1994) Method
model.BothUnrepBurst.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication
model.UnrepBurst.UnrepBurst2.pvec <- numeric(NSIMS) # vector for pvalues - Maddison and FitzJohn Replication

# Threshold Model - Positive Scenario
mcmc.bothunrepburst.R <- numeric(NSIMS) # vector for R values - Darwin Scenario Threshold Model
mcmc.bothunrepburst.R2 <- numeric(NSIMS) # vector for R^2 values - Darwin Scenario Threshold Model

# Threshold Model - Negative Scenario
mcmc.UnrepBurst.UnrepBurst2.R <- numeric(NSIMS) # vector for R values - UnrepBurst Scenario Threshold Model
mcmc.UnrepBurst.UnrepBurst2.R2 <- numeric(NSIMS) # vector for R^2 values - UnrepBurst Scenario Threshold Model

# MCMCglmm Model
model.mcmcglmm.BothUnrepBurst.pvec <- numeric(NSIMS) # vector for "pvalues" - MCMCglmm Model
model.mcmcglmm.UnrepBurst.UnrepBurst2.pvec <- numeric(NSIMS) # vector for "pvalues" - MCMCglmm Model
solutions <- numeric(NSIMS) #summary stats
HPD.Sol.Both <- numeric(NSIMS)
HPD.VCV.Both <- numeric(NSIMS)
HPD.Sol.Unrep <- numeric(NSIMS)
HPD.VCV.Unrep <- numeric(NSIMS)


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


# ##MCMCglmm Model in MCMCglmm
PA<-inverseA(phy, nodes="TIPS")$Ainv ###inverted tree matrix
prior1 <-list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) #set prior, expected covariances: G = random (phylogenetic) effects, R = residuals, V = variance, nu = , alpha.mu = , alpha.V = variance structure component

	#MCMCglmm Darwin Scenario Data
darwindata.MatrixRownames <- cbind(Taxa = rownames(darwindata.MatrixBind), darwindata.MatrixBind)
darwindata.dataframe <- as.data.frame(darwindata.MatrixRownames)

	#MCMCglmm Unreplicated Burst Scenario Data
UnrepBurst.MatrixRownames <- cbind(Taxa = rownames(UnrepBurst.MatrixBind), UnrepBurst.MatrixBind)
UnrepBurst.dataframe <- as.data.frame(UnrepBurst.MatrixRownames)

	#MCMCglmm Positive Control Data
BothUnrepBurst.MatrixRownames <- cbind(Taxa = rownames(BothUnrepBurst.MatrixBind), BothUnrepBurst.MatrixBind)
BothUnrepBurst.dataframe <- as.data.frame(BothUnrepBurst.MatrixRownames)

	#MCMCglmm Negative Control Data
UnrepBurst.UnrepBurst2.MatrixRownames <- cbind(Taxa = rownames(UnrepBurst.UnrepBurst2.MatrixBind), UnrepBurst.UnrepBurst2.MatrixBind)
UnrepBurst.UnrepBurst2.dataframe <- as.data.frame(UnrepBurst.UnrepBurst2.MatrixRownames)

	#MCMCglmm Positive Control Analysis
BothUnrepBurst.dataframe[,2]<- as.numeric(as.character(BothUnrepBurst.dataframe[,2])) #change second column from a factor to a numeric
BothUnrepBurst.dataframe[,3]<- as.numeric(as.character(BothUnrepBurst.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.BothUnrepBurst <-MCMCglmm(V2~V3, random =~Taxa, family="categorical", ginverse = list(Taxa = PA), prior = prior1, nitt=1000000, thin=1000, burnin=250000, data=BothUnrepBurst.dataframe)

	#MCMCglmm Negative Control Analysis
UnrepBurst.UnrepBurst2.dataframe[,2]<- as.numeric(as.character(UnrepBurst.UnrepBurst2.dataframe[,2])) #change second column from a factor to a numeric
UnrepBurst.UnrepBurst2.dataframe[,3]<- as.numeric(as.character(UnrepBurst.UnrepBurst2.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.UnrepBurst.UnrepBurst2 <-MCMCglmm(V2~V3, random =~Taxa, family="categorical", ginverse = list(Taxa = PA), prior = prior1, nitt=1000000, thin=1000, burnin=250000, data=UnrepBurst.UnrepBurst2.dataframe)

	#Saves data for the MCMCglmm models
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)

		##Positive P-vals
s <- summary(MCMCglmm.BothUnrepBurst)
model.mcmcglmm.BothUnrepBurst.pvec[i] <- s[["solutions"]][[10]] #P-Value for Positive Control - should have low p-val
HPD.Sol.Both [i] <- HPDinterval(MCMCglmm.BothUnrepBurst$Sol)
HPD.VCV.Both [i] <- HPDinterval(MCMCglmm.BothUnrepBurst$VCV)

		##Negative P-vals
s2 <- summary(MCMCglmm.UnrepBurst.UnrepBurst2)
model.mcmcglmm.UnrepBurst.UnrepBurst2.pvec[i] <- s2[["solutions"]][[10]] #P-Value for Negative Control - should have high p-val
HPD.Sol.Unrep [i] <- HPDinterval(MCMCglmm.UnrepBurst.UnrepBurst2$Sol)
HPD.VCV.Unrep [i] <- HPDinterval(MCMCglmm.UnrepBurst.UnrepBurst2$VCV)

}

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##Maddison Results
length(trees)
length(cladesize)
cladesize
length(DarwinCharacters)
length(UnrepBurstCharacters)
length(UnrepBurstCharacters.mean)
hist(UnrepBurstCharacters.mean)

##Pagel's Discrete Model Results
model.BothUnrepBurst.pvec #list of all p-values
hist(model.BothUnrepBurst.pvec) #positive control histogram
model.UnrepBurst.UnrepBurst2.pvec #list of all p-values
hist(model.UnrepBurst.UnrepBurst2.pvec) #negative control histogram

##Threshold Results
mcmc.bothunrepburst.R
hist(mcmc.bothunrepburst.R)
median(mcmc.bothunrepburst.R) #median r value
boxplot(mcmc.bothunrepburst.R) #box plot of posterior distribution of r
plot(mcmc.bothunrepburst$par[, "logL"], type = "l", xlab = "generation", ylab = "logL")

mcmc.UnrepBurst.UnrepBurst2.R
hist(mcmc.UnrepBurst.UnrepBurst2.R)
median(mcmc.UnrepBurst.UnrepBurst2.R) #median r value
boxplot(mcmc.UnrepBurst.UnrepBurst2.R) #box plot of posterior distribution of r
plot(mcmc.UnrepBurst.UnrepBurst2$par[, "logL"], type = "l", xlab = "generation", ylab = "logL")

##MCMCglmm Results
autocorr(MCMCglmm.BothUnrepBurst$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(log(MCMCglmm.BothUnrepBurst$VCV)) #plot of the posterior distribution of the variance (co)variance matrix
autocorr(MCMCglmm.BothUnrepBurst$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.BothUnrepBurst$Sol) #plot of the posterior distribution of the fixed effect
model.mcmcglmm.BothUnrepBurst.pvec #list of all p-values
hist(model.mcmcglmm.BothUnrepBurst.pvec) #histogram of posterior distribution of p-values
HPD.Sol.Both
HPD.VCV.Both

autocorr(MCMCglmm.UnrepBurst.UnrepBurst2$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(log(MCMCglmm.UnrepBurst.UnrepBurst2$VCV)) #plot of the posterior distribution of the variance (co)variance matrix
autocorr(MCMCglmm.UnrepBurst.UnrepBurst2$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.UnrepBurst.UnrepBurst2$Sol) #plot of the posterior distribution of the fixed effect
model.mcmcglmm.UnrepBurst.UnrepBurst2.pvec #list of all p-values
hist(model.mcmcglmm.UnrepBurst.UnrepBurst2.pvec) #histogram of posterior distribution of p-values
HPD.Sol.Unrep
HPD.VCV.Unrep
