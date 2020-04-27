# Discrete correlation script
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/24/2017

#library(doParallel)
start.time <- Sys.time()
#Parallel On Local Computer
#cl<-makeCluster(8)# for my desktop with 16 cores
#registerDoParallel(cl)

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

# MCMCglmm Model
Darwin.Interc.Post.Mean <- numeric(NSIMS) # Vector for Darwin scenario posterior distribution intercept mean
Darwin.Estim.Post.Mean <- numeric(NSIMS) # Vector for Darwin scenario posterior distribution estimate mean
Darwin.Interc.Lower.CI <- numeric(NSIMS) # Vector for Darwin scenario intercept 95% credible interval lower bound
Darwin.Estim.Lower.CI <- numeric(NSIMS) # Vector for Darwin scenario estimate 95% credible interval lower bound
Darwin.Interc.Upper.CI <- numeric(NSIMS) # Vector for Darwin scenario intercept 95% credible interval upper bound
Darwin.Estim.Upper.CI <- numeric(NSIMS) # Vector for Darwin scenario estimate 95% credible interval upper bound
Darwin.Interc.pMCMC <- numeric(NSIMS) # Vector for Darwin scenario intercept pMCMC ("p-value")
Darwin.Estim.pMCMC <- numeric(NSIMS) # Vector for Darwin scenario estimate pMCMC ("p-value")

Unrep.Interc.Post.Mean <- numeric(NSIMS) # Vector for Unreplicated Burst scenario posterior distribution intercept mean
Unrep.Estim.Post.Mean <- numeric(NSIMS) # Vector for Unreplicated Burst scenario posterior distribution estimate mean
Unrep.Interc.Lower.CI <- numeric(NSIMS) # Vector for Unreplicated Burst scenario intercept 95% credible interval lower bound
Unrep.Estim.Lower.CI <- numeric(NSIMS) # Vector for Unreplicated Burst scenario estimate 95% credible interval lower bound
Unrep.Interc.Upper.CI <- numeric(NSIMS) # Vector for Unreplicated Burst scenario intercept 95% credible interval upper bound
Unrep.Estim.Upper.CI <- numeric(NSIMS) # Vector for Unreplicated Burst scenario estimate 95% credible interval upper bound
Unrep.Interc.pMCMC <- numeric(NSIMS) # Vector for Unreplicated Burst scenario intercept pMCMC ("p-value")
Unrep.Estim.pMCMC <- numeric(NSIMS) # Vector for Unreplicated Burst scenario estimate pMCMC ("p-value")

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

	#ThreshBayes Darwin Scenario Data
darwindata.Matrix <- matrix(c(darwindata))
darwindata.MatrixBind <- cbind(darwindata.Matrix,darwindata.Matrix)
dimnames(darwindata.MatrixBind) <- list(as.list(phy$tip.label))

	#ThreshBayes Unreplicated Burst Scenario Data
UnrepBurst.Matrix <- matrix(c(unrepburstdata))
UnrepBurst.MatrixBind <- cbind(darwindata.Matrix,UnrepBurst.Matrix)
dimnames(UnrepBurst.MatrixBind) <- list(as.list(phy$tip.label))


# ##MCMCglmm Model in MCMCglmm
Ainv<-inverseA(phy, nodes="TIPS")$Ainv ###inverted tree matrix
#prior1 <-list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, fix = 1), G = list(G1 = list(V=2,nu=2)))
#prior2 <-list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, nu = 1))
prior3 <- list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

 #set prior, expected covariances: B = prior for the mean, G = random effects, R = prior for the variance, V = variance, nu = degree of belief parameter, mu = mean, alpha.V = variance structure component

	#MCMCglmm Darwin Scenario Data
darwindata.MatrixRownames <- cbind(animal = rownames(darwindata.MatrixBind), darwindata.MatrixBind)
darwindata.dataframe <- as.data.frame(darwindata.MatrixRownames)

	#MCMCglmm Unreplicated Burst Scenario Data
UnrepBurst.MatrixRownames <- cbind(animal = rownames(UnrepBurst.MatrixBind), UnrepBurst.MatrixBind)
UnrepBurst.dataframe <- as.data.frame(UnrepBurst.MatrixRownames)

	#MCMCglmm Darwin Scenario Analysis
darwindata.dataframe[,2]<- as.numeric(as.character(darwindata.dataframe[,2])) #change second column from a factor to a numeric
darwindata.dataframe[,3]<- as.numeric(as.character(darwindata.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.darwin <-MCMCglmm(V2~V3, random =~animal, family="categorical", ginverse = list(animal = Ainv), prior = prior3, nitt=5000000, thin=2500, burnin=1250000, data=darwindata.dataframe)
#MCMCglmm.data <-MCMCglmm(y~x, random =~animal, family="gaussian", ginverse = list(animal = Ainv), prior = prior, nitt=1000000, thin=100, burnin=250000, data=dataframe)

	#MCMCglmm Unreplicated Burst Scenario Analysis
UnrepBurst.dataframe[,2]<- as.numeric(as.character(UnrepBurst.dataframe[,2])) #change second column from a factor to a numeric
UnrepBurst.dataframe[,3]<- as.numeric(as.character(UnrepBurst.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.unrepburst <-MCMCglmm(V2~V3, random =~animal, family="categorical", ginverse = list(animal = Ainv), prior = prior3, nitt=5000000, thin=2500, burnin=1250000, data=UnrepBurst.dataframe)

	#Saves data for the MCMCglmm models
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips)
s <- summary(MCMCglmm.darwin)
Darwin.Interc.Post.Mean [i] <- s[["solutions"]][[1]]
Darwin.Estim.Post.Mean [i] <- s[["solutions"]][[2]]
Darwin.Interc.Lower.CI [i] <- s[["solutions"]][[3]]
Darwin.Estim.Lower.CI [i] <- s[["solutions"]][[4]]
Darwin.Interc.Upper.CI [i] <- s[["solutions"]][[5]]
Darwin.Estim.Upper.CI [i] <- s[["solutions"]][[6]]
Darwin.Interc.pMCMC [i] <- s[["solutions"]][[9]]
Darwin.Estim.pMCMC [i] <- s[["solutions"]][[10]]

s2 <- summary(MCMCglmm.unrepburst)
Unrep.Interc.Post.Mean [i] <- s2[["solutions"]][[1]]
Unrep.Estim.Post.Mean [i] <- s2[["solutions"]][[2]]
Unrep.Interc.Lower.CI [i] <- s2[["solutions"]][[3]]
Unrep.Estim.Lower.CI [i] <- s2[["solutions"]][[4]]
Unrep.Interc.Upper.CI [i] <- s2[["solutions"]][[5]]
Unrep.Estim.Upper.CI [i] <- s2[["solutions"]][[6]]
Unrep.Interc.pMCMC[i] <- s2[["solutions"]][[9]]
Unrep.Estim.pMCMC[i] <- s2[["solutions"]][[10]]
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
#stopCluster(cl)

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results
##General Results
length(trees)
length(cladesize)
length(DarwinCharacters)
length(UnrepBurstCharacters)
length(UnrepBurstCharacters.mean)
length(Darwin.Estim.pMCMC)
length(Unrep.Estim.pMCMC)
length(Darwin.Estim.pMCMC <= 0.05)
length(Unrep.Estim.pMCMC <= 0.05)

##MCMCglmm Darwin Scenario Results
autocorr(MCMCglmm.darwin$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(log(MCMCglmm.darwin$VCV)) #plot of the posterior distribution of the variance
autocorr(MCMCglmm.darwin$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.darwin$Sol) #plot of the posterior distribution of the coefficients
Darwin.Interc.Post.Mean
Darwin.Estim.Post.Mean
Darwin.Interc.Lower.CI
Darwin.Estim.Lower.CI
Darwin.Interc.Upper.CI
Darwin.Estim.Upper.CI
Darwin.Interc.pMCMC
Darwin.Estim.pMCMC
median(Darwin.Interc.Post.Mean)
median(Darwin.Estim.Post.Mean)
median(Darwin.Interc.Lower.CI)
median(Darwin.Estim.Lower.CI)
median(Darwin.Interc.Upper.CI)
median(Darwin.Estim.Upper.CI)
median(Darwin.Interc.pMCMC)
median(Darwin.Estim.pMCMC)
hist(Darwin.Estim.pMCMC) #histogram of posterior distribution of pMCMC-values
plot(cladesize, Darwin.Estim.pMCMC)
MCMCglmm.Darwin.lm <- lm(Darwin.Estim.pMCMC ~ cladesize)
summary(MCMCglmm.Darwin.lm)

##MCMCglmm Unreplicated Burst Scenario Results
autocorr(MCMCglmm.unrepburst$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(log(MCMCglmm.unrepburst$VCV)) #plot of the posterior distribution of the variance
autocorr(MCMCglmm.unrepburst$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.unrepburst$Sol) #plot of the posterior distribution of the coefficients
Unrep.Interc.Post.Mean
Unrep.Estim.Post.Mean
Unrep.Interc.Lower.CI
Unrep.Estim.Lower.CI
Unrep.Interc.Upper.CI
Unrep.Estim.Upper.CI
Unrep.Interc.pMCMC
Unrep.Estim.pMCMC
median(Unrep.Interc.Post.Mean)
median(Unrep.Estim.Post.Mean)
median(Unrep.Interc.Lower.CI)
median(Unrep.Estim.Lower.CI)
median(Unrep.Interc.Upper.CI)
median(Unrep.Estim.Upper.CI)
median(Unrep.Interc.pMCMC)
median(Unrep.Estim.pMCMC)
hist(Unrep.Estim.pMCMC) #histogram of posterior distribution of pMCMC-values
plot(cladesize, Unrep.Estim.pMCMC)
MCMCglmm.Unrep.lm <- lm(Unrep.Estim.pMCMC ~ cladesize)
summary(MCMCglmm.Unrep.lm)
