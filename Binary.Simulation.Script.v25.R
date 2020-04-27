# Discrete correlation script
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 12/27/2016

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

# Threshold Model - Darwin Scenario
mcmcdarwin.r <- numeric(NSIMS) # vector for R values - Darwin Scenario Threshold Model
mcmcdarwin.r2 <- numeric(NSIMS) # vector for R^2 values - Darwin Scenario Threshold Model

# Threshold Model - Unreplicated Burst Scenario
mcmcunrepburst.r <- numeric(NSIMS) # vector for R values - UnrepBurst Scenario Threshold Model
mcmcunrepburst.r2 <- numeric(NSIMS) # vector for R^2 values - UnrepBurst Scenario Threshold Model

# MCMCglmm Model
model.mcmcglmm.darwin.pvec <- numeric(NSIMS) # vector for pMCMC values ("p-values") - MCMCglmm Model
model.mcmcglmm.unrepburst.pvec <- numeric(NSIMS) # vector for pMCMC values ("p-values") - MCMCglmm Model
HPD.Sol.Darwin <- numeric(NSIMS)
HPD.VCV.Darwin <- numeric(NSIMS)
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

#Unrep Burst 2 (divide tips by 8)
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

#Unrep Burst 3 (divide tips by 6)
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

#Unrep Burst 5 (divide by 2)
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

##Saves data in a table for each iteration. For explanations of data structure see Maddison and Fitzjohn. Syst. Biol. 64(1):127?136, 2015.
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
mcmc.darwin <- threshBayes(phy,darwindata.MatrixBind,types=c("discrete"),ngen=ngen)
burnin.darwin <- colMeans(mcmc.darwin$par[10:nrow(mcmc.darwin$par),2:6])

	#ThreshBayes Unreplicated Burst Scenario Analysis
mcmc.unrepburst <- threshBayes(phy,UnrepBurst.MatrixBind,types=c("discrete"),ngen=ngen)
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


# ##MCMCglmm Model in MCMCglmm
Ainv<-inverseA(phy, nodes="TIPS")$Ainv ###inverted tree matrix
prior <- list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))
prior1 <-list(R=list(V=10, fix =1), G=list(G1=list(V=1,nu=1, alpha.mu = 0, alpha.V = 1000))) #set prior
prior2<-list(B=list(mu=c(0,0), V=gelman.prior(~V3, data=darwindata.dataframe, scale=sqrt(pi^2/3+1))),R=list(V=1,fix=1),G=list(G1=list(V=1,nu=1, alpha.mu = 0, alpha.V = 1000)))
prior3 <- list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

	#MCMCglmm Darwin Scenario Data
darwindata.MatrixRownames <- cbind(animal = rownames(darwindata.MatrixBind), darwindata.MatrixBind)
darwindata.dataframe <- as.data.frame(darwindata.MatrixRownames)

	#MCMCglmm Unreplicated Burst Scenario Data
UnrepBurst.MatrixRownames <- cbind(animal = rownames(UnrepBurst.MatrixBind), UnrepBurst.MatrixBind)
UnrepBurst.dataframe <- as.data.frame(UnrepBurst.MatrixRownames)

	#MCMCglmm Darwin Scenario Analysis
MCMCglmm.darwin <-MCMCglmm(V2~V3, random =~animal, family="categorical", ginverse = list(animal = Ainv), prior = prior3, nitt=20000000, thin=7500, burnin=5000000, data=darwindata.dataframe)

	#MCMCglmm Unreplicated Burst Scenario Analysis
UnrepBurst.dataframe[,2]<- as.numeric(as.character(UnrepBurst.dataframe[,2])) #change second column from a factor to a numeric
UnrepBurst.dataframe[,3]<- as.numeric(as.character(UnrepBurst.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.unrepburst <-MCMCglmm(V2~V3, random =~animal, family="categorical", ginverse = list(animal = Ainv), prior = prior3, nitt=20000000, thin=7500, burnin=5000000, data=UnrepBurst.dataframe)

	#Saves data for the MCMCglmm models
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)
s <- summary(MCMCglmm.darwin)
model.mcmcglmm.darwin.pvec[i] <- s[["solutions"]][[10]]
HPD.Sol.Darwin [i] <- HPDinterval(MCMCglmm.darwin$Sol)
HPD.VCV.Darwin [i] <- HPDinterval(MCMCglmm.darwin$VCV)
s2 <- summary(MCMCglmm.unrepburst)
model.mcmcglmm.unrepburst.pvec[i] <- s2[["solutions"]][[10]]
HPD.Sol.Unrep [i] <- HPDinterval(MCMCglmm.unrepburst$Sol)
HPD.VCV.Unrep [i] <- HPDinterval(MCMCglmm.unrepburst$VCV)

}

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##General Results
length(trees)
length(cladesize)
cladesize
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

##Pagel Unreplicated Burst Scenario Results
model.Darwin.UnrepBurst.pvec
hist(model.Darwin.UnrepBurst.pvec)
min(model.Darwin.UnrepBurst.pvec)
max(model.Darwin.UnrepBurst.pvec)
median(model.Darwin.UnrepBurst.pvec)
plot(UnrepBurstCharacters.mean, model.Darwin.UnrepBurst.pvec)
plot(cladesize, model.Darwin.UnrepBurst.pvec)

##Pagel Progressive Departure Results
model.Darwin.UnrepBurst2.pvec
hist(model.Darwin.UnrepBurst2.pvec)
median(model.Darwin.UnrepBurst2.pvec)
model.Darwin.UnrepBurst3.pvec
hist(model.Darwin.UnrepBurst3.pvec)
median(model.Darwin.UnrepBurst3.pvec)
model.Darwin.UnrepBurst4.pvec
hist(model.Darwin.UnrepBurst4.pvec)
median(model.Darwin.UnrepBurst4.pvec)
model.Darwin.UnrepBurst5.pvec
hist(model.Darwin.UnrepBurst5.pvec)
median(model.Darwin.UnrepBurst5.pvec)

##Threshold Darwin Scenario Results
mcmcdarwin.r
hist(mcmcdarwin.r)
boxplot(mcmcdarwin.r)
min(mcmcdarwin.r)
max(mcmcdarwin.r)
median(mcmcdarwin.r)

##Threshold Unreplicated Burst Scenario Results
mcmcunrepburst.r
hist(mcmcunrepburst.r)
boxplot(mcmcunrepburst.r)
min(mcmcunrepburst.r)
max(mcmcunrepburst.r)
median(mcmcunrepburst.r)

##MCMCglmm Darwin Scenario Results
summary(MCMCglmm.darwin)
autocorr(MCMCglmm.darwin$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.darwin$VCV) #plot of the posterior distribution of the variance (co)variance matrix
autocorr(MCMCglmm.darwin$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.darwin$Sol) #plot of the posterior distribution of the fixed effect
model.mcmcglmm.darwin.pvec #list of p-values
hist(model.mcmcglmm.darwin.pvec) #histogram of posterior distribution of p-values
HPD.Sol.Darwin
HPD.VCV.Darwin

##MCMCglmm Unreplicated Burst Scenario Results
summary(MCMCglmm.unrepburst)
autocorr(MCMCglmm.unrepburst$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.unrepburst$VCV) #plot of the posterior distribution of the variance (co)variance matrix
autocorr(MCMCglmm.unrepburst$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.unrepburst$Sol) #plot of the posterior distribution of the fixed effect
model.mcmcglmm.unrepburst.pvec #list of p-values
hist(model.mcmcglmm.unrepburst.pvec) #histogram of posterior distribution of p-values
HPD.Sol.Unrep
HPD.VCV.Unrep
