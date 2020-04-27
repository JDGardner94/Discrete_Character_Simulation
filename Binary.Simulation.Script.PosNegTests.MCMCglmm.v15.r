# Discrete correlation script - Control Tests
# Chris Organ, 12/30/2014
# edited by Jacob Gardner, 3/18/2017

start.time <- Sys.time()

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

# MCMCglmm Model
Positive.Interc.Post.Mean <- numeric(NSIMS) # Vector for Positive control scenario posterior distribution intercept mean
Positive.Estim.Post.Mean <- numeric(NSIMS) # Vector for Positive control scenario posterior distribution estimate mean
Positive.Interc.Lower.CI <- numeric(NSIMS) # Vector for Positive control scenario intercept 95% credible interval lower bound
Positive.Estim.Lower.CI <- numeric(NSIMS) # Vector for Positive control scenario estimate 95% credible interval lower bound
Positive.Interc.Upper.CI <- numeric(NSIMS) # Vector for Positive control scenario intercept 95% credible interval upper bound
Positive.Estim.Upper.CI <- numeric(NSIMS) # Vector for Positive control scenario estimate 95% credible interval upper bound
Positive.Interc.pMCMC <- numeric(NSIMS) # Vector for Positive control scenario intercept pMCMC ("p-value")
Positive.Estim.pMCMC <- numeric(NSIMS) # Vector for Positive control scenario estimate pMCMC ("p-value")
Positive.VCV.Stationarity <- numeric(NSIMS) # Vector for Heidelberg Welch stationarity convergence test for variance (0 = fail, 1 = pass)
Positive.VCV.Halfwidth <- numeric(NSIMS) # Vector for Heidelberg Welch halfwidth convergence test for variance (0 = fail, 1 = pass)
Positive.Sol.Pass <- list(NSIMS) # List of Heidelberg Welch convergence tests for fixed effects, is "passed" present? (True = Passed, False = Failed)
Positive.Sol.Fail <- list(NSIMS) # List of Heidelberg Welch convergence tests for fixed effects, is "failed" present? (True = Failed, False = Passed)

Negative.Interc.Post.Mean <- numeric(NSIMS) # Vector for Negative control scenario posterior distribution intercept mean
Negative.Estim.Post.Mean <- numeric(NSIMS) # Vector for Negative control scenario posterior distribution estimate mean
Negative.Interc.Lower.CI <- numeric(NSIMS) # Vector for Negative control scenario intercept 95% credible interval lower bound
Negative.Estim.Lower.CI <- numeric(NSIMS) # Vector for Negative control scenario estimate 95% credible interval lower bound
Negative.Interc.Upper.CI <- numeric(NSIMS) # Vector for Negative control scenario intercept 95% credible interval upper bound
Negative.Estim.Upper.CI <- numeric(NSIMS) # Vector for Negative control scenario estimate 95% credible interval upper bound
Negative.Interc.pMCMC <- numeric(NSIMS) # Vector for Negative control scenario intercept pMCMC ("p-value")
Negative.Estim.pMCMC <- numeric(NSIMS) # Vector for Negative control scenario estimate pMCMC ("p-value")
Negative.VCV.Stationarity <- numeric(NSIMS) # Vector for Heidelberg Welch stationarity convergence test for variance (0 = fail, 1 = pass)
Negative.VCV.Halfwidth <- numeric(NSIMS) # Vector for Heidelberg Welch halfwidth convergence test for variance (0 = fail, 1 = pass)
Negative.Sol.Pass <- list(NSIMS) # List of Heidelberg Welch convergence tests for fixed effects, is "passed" present? (True = Passed, False = Failed)
Negative.Sol.Fail <- list(NSIMS) # List of Heidelberg Welch convergence tests for fixed effects, is "failed" present? (True = Failed, False = Passed)


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


# ##MCMCglmm Model in MCMCglmm
Ainv<-inverseA(phy, nodes="TIPS")$Ainv ###inverted tree matrix
#prior1 <-list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, fix = 1)) #set prior, expected covariances: G = random (phylogenetic) effects, R = residuals, V = variance, nu = , alpha.mu = , alpha.V = variance structure component
prior3 <- list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

	#MCMCglmm Darwin Scenario Data
darwindata.MatrixRownames <- cbind(animal = rownames(darwindata.MatrixBind), darwindata.MatrixBind)
darwindata.dataframe <- as.data.frame(darwindata.MatrixRownames)

	#MCMCglmm Unreplicated Burst Scenario Data
UnrepBurst.MatrixRownames <- cbind(animal = rownames(UnrepBurst.MatrixBind), UnrepBurst.MatrixBind)
UnrepBurst.dataframe <- as.data.frame(UnrepBurst.MatrixRownames)

	#MCMCglmm Positive Control Data
BothUnrepBurst.MatrixRownames <- cbind(animal = rownames(BothUnrepBurst.MatrixBind), BothUnrepBurst.MatrixBind)
BothUnrepBurst.dataframe <- as.data.frame(BothUnrepBurst.MatrixRownames)

	#MCMCglmm Negative Control Data
UnrepBurst.UnrepBurst2.MatrixRownames <- cbind(animal = rownames(UnrepBurst.UnrepBurst2.MatrixBind), UnrepBurst.UnrepBurst2.MatrixBind)
UnrepBurst.UnrepBurst2.dataframe <- as.data.frame(UnrepBurst.UnrepBurst2.MatrixRownames)

	#MCMCglmm Positive Control Analysis
BothUnrepBurst.dataframe[,2]<- as.numeric(as.character(BothUnrepBurst.dataframe[,2])) #change second column from a factor to a numeric
BothUnrepBurst.dataframe[,3]<- as.numeric(as.character(BothUnrepBurst.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.BothUnrepBurst <-MCMCglmm(V2~V3, random=~animal, family="categorical", ginverse = list(animal = Ainv), prior = prior3, nitt=5000000, thin=2500, burnin=1250000, data=BothUnrepBurst.dataframe)

	#MCMCglmm Negative Control Analysis
UnrepBurst.UnrepBurst2.dataframe[,2]<- as.numeric(as.character(UnrepBurst.UnrepBurst2.dataframe[,2])) #change second column from a factor to a numeric
UnrepBurst.UnrepBurst2.dataframe[,3]<- as.numeric(as.character(UnrepBurst.UnrepBurst2.dataframe[,3])) #change third column from a factor to a numeric
MCMCglmm.UnrepBurst.UnrepBurst2 <-MCMCglmm(V2~V3, random=~animal, family="categorical", ginverse = list(animal = Ainv), prior = prior3, nitt=5000000, thin=2500, burnin=1250000, data=UnrepBurst.UnrepBurst2.dataframe)

	#Saves data for the MCMCglmm models
trees[i] <- write.tree(phy)
cladesize[i] <- length(tips)
DarwinCharacters <- rbind(DarwinCharacters,darwindata)
UnrepBurstCharacters <- rbind(UnrepBurstCharacters,unrepburstdata)
UnrepBurstCharacters.mean[i] <- length(unrepbursttips)/length(tips) 
s <- summary(MCMCglmm.BothUnrepBurst)
Positive.Interc.Post.Mean [i] <- s[["solutions"]][[1]]
Positive.Estim.Post.Mean [i] <- s[["solutions"]][[2]]
Positive.Interc.Lower.CI [i] <- s[["solutions"]][[3]]
Positive.Estim.Lower.CI [i] <- s[["solutions"]][[4]]
Positive.Interc.Upper.CI [i] <- s[["solutions"]][[5]]
Positive.Estim.Upper.CI [i] <- s[["solutions"]][[6]]
Positive.Interc.pMCMC [i] <- s[["solutions"]][[9]]
Positive.Estim.pMCMC [i] <- s[["solutions"]][[10]]
Positive.VCV.Heidel <- heidel.diag(MCMCglmm.BothUnrepBurst$VCV)
Positive.VCV.Stationarity [i] <- as.vector(Positive.VCV.Heidel)[[1]]
Positive.VCV.Halfwidth [i] <- as.vector(Positive.VCV.Heidel)[[7]]
Positive.Sol.Heidel <- heidel.diag(MCMCglmm.BothUnrepBurst$Sol)
Positive.Sol.Pass [i] <- '1' %in% as.vector(Positive.Sol.Heidel)
Positive.Sol.Fail [i] <- '0' %in% as.vector(Positive.Sol.Heidel)

s2 <- summary(MCMCglmm.UnrepBurst.UnrepBurst2)
Negative.Interc.Post.Mean [i] <- s2[["solutions"]][[1]]
Negative.Estim.Post.Mean [i] <- s2[["solutions"]][[2]]
Negative.Interc.Lower.CI [i] <- s2[["solutions"]][[3]]
Negative.Estim.Lower.CI [i] <- s2[["solutions"]][[4]]
Negative.Interc.Upper.CI [i] <- s2[["solutions"]][[5]]
Negative.Estim.Upper.CI [i] <- s2[["solutions"]][[6]]
Negative.Interc.pMCMC[i] <- s2[["solutions"]][[9]]
Negative.Estim.pMCMC[i] <- s2[["solutions"]][[10]]
Negative.VCV.Heidel <- heidel.diag(MCMCglmm.UnrepBurst.UnrepBurst2$VCV)
Negative.VCV.Stationarity [i] <- as.vector(Negative.VCV.Heidel)[[1]]
Negative.VCV.Halfwidth [i] <- as.vector(Negative.VCV.Heidel)[[7]]
Negative.Sol.Heidel <- heidel.diag(MCMCglmm.UnrepBurst.UnrepBurst2$Sol)
Negative.Sol.Pass [i] <- '1' %in% as.vector(Negative.Sol.Heidel)
Negative.Sol.Fail [i] <- '0' %in% as.vector(Negative.Sol.Heidel)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

DarwinCharacters <- DarwinCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok
UnrepBurstCharacters <- UnrepBurstCharacters[-1,] #remove the first row of data (all 1's) from the store table. The underlying data are ok

###Summarize results

##MCMCglmm Positive Control Scenario Results
autocorr(MCMCglmm.BothUnrepBurst$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(log(MCMCglmm.BothUnrepBurst$VCV)) #plot of the posterior distribution of the variance (co)variance matrix
Positive.VCV.Stationarity
Positive.VCV.Halfwidth
autocorr(MCMCglmm.BothUnrepBurst$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.BothUnrepBurst$Sol) #plot of the posterior distribution of the fixed effect
Positive.Sol.Pass
Positive.Sol.Fail
Positive.Interc.Post.Mean 
Positive.Estim.Post.Mean 
Positive.Interc.Lower.CI 
Positive.Estim.Lower.CI 
Positive.Interc.Upper.CI 
Positive.Estim.Upper.CI 
Positive.Interc.pMCMC 
Positive.Estim.pMCMC 
median(Positive.Interc.Post.Mean)
median(Positive.Estim.Post.Mean)
median(Positive.Interc.Lower.CI)
median(Positive.Estim.Lower.CI)
median(Positive.Interc.Upper.CI)
median(Positive.Estim.Upper.CI)
median(Positive.Interc.pMCMC)
median(Positive.Estim.pMCMC)
hist(Positive.Estim.pMCMC) #histogram of posterior distribution of pMCMC-values

##MCMCglmm Negative Control Scenario Results
autocorr(MCMCglmm.UnrepBurst.UnrepBurst2$VCV) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(log(MCMCglmm.UnrepBurst.UnrepBurst2$VCV)) #plot of the posterior distribution of the variance (co)variance matrix
Negative.VCV.Stationarity
Negative.VCV.Halfwidth
autocorr(MCMCglmm.UnrepBurst.UnrepBurst2$Sol) #measure of the autocorrelation effect - values of Lag 10 and onward should be close to 0
plot(MCMCglmm.UnrepBurst.UnrepBurst2$Sol) #plot of the posterior distribution of the fixed effect
Negative.Sol.Pass
Negative.Sol.Fail
Negative.Interc.Post.Mean 
Negative.Estim.Post.Mean 
Negative.Interc.Lower.CI 
Negative.Estim.Lower.CI 
Negative.Interc.Upper.CI 
Negative.Estim.Upper.CI 
Negative.Interc.pMCMC 
Negative.Estim.pMCMC 
median(Negative.Interc.Post.Mean)
median(Negative.Estim.Post.Mean)
median(Negative.Interc.Lower.CI)
median(Negative.Estim.Lower.CI)
median(Negative.Interc.Upper.CI)
median(Negative.Estim.Upper.CI)
median(Negative.Interc.pMCMC)
median(Negative.Estim.pMCMC)
hist(Negative.Estim.pMCMC) #histogram of posterior distribution of pMCMC-values



