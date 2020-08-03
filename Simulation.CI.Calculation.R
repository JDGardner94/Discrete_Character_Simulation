### Discrete simulation code
## Calculating consistency index
# Last edited: February 29, 2020, by Jacob Gardner

#install.packages("ape")
#install.packages("phangorn")
library("ape")
library("phangorn")

# Vectors for CI
NSIMS <- 1000
CI.Positive <- numeric(NSIMS)
CI.Negative <- numeric(NSIMS)
CI.Darwin <- numeric(NSIMS)
CI.Unrepburst <- numeric(NSIMS)
ClassImb.Positive <- numeric(NSIMS)
ClassImb.Negative <- numeric(NSIMS)
ClassImb.Darwin <- numeric(NSIMS)
ClassImb.Unrepburst <- numeric(NSIMS)
MultiImb.Positive <- numeric(NSIMS)
MultiImb.Negative <- numeric(NSIMS)
MultiImb.Darwin <- numeric(NSIMS)
MultiImb.Unrepburst <- numeric(NSIMS)
Positive.Class1 <- numeric(NSIMS)
Positive.Class2 <- numeric(NSIMS)
Positive.Class3 <- numeric(NSIMS)
Positive.Class4 <- numeric(NSIMS)
Negative.Class1 <- numeric(NSIMS)
Negative.Class2 <- numeric(NSIMS)
Negative.Class3 <- numeric(NSIMS)
Negative.Class4 <- numeric(NSIMS)
Darwin.Class1 <- numeric(NSIMS)
Darwin.Class2 <- numeric(NSIMS)
Darwin.Class3 <- numeric(NSIMS)
Darwin.Class4 <- numeric(NSIMS)
Unrepburst.Class1 <- numeric(NSIMS)
Unrepburst.Class2 <- numeric(NSIMS)
Unrepburst.Class3 <- numeric(NSIMS)
Unrepburst.Class4 <- numeric(NSIMS)


# Calculate CI - Positive Control
Positive.Files <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Control/RJMCMC_PositiveData/", pattern =".txt")
Positive.Trees <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Control/RJMCMC_PositiveData/", pattern =".nex")
setwd("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Control/RJMCMC_PositiveData/")
for(i in 1:length(Positive.Files)){
  PositiveData <- read.delim(Positive.Files[i], header=FALSE, row.names=1)
  tPositiveData <- t(PositiveData)
  PositiveFrame <- as.data.frame(tPositiveData, row.names = c("1","2"))
  PositiveFrame2 <- as.data.frame(PositiveData)
  Positive.Class1[i] <- nrow(PositiveFrame2[PositiveFrame2$V2==0 & PositiveFrame2$V3==0, ])
  Positive.Class2[i] <- nrow(PositiveFrame2[PositiveFrame2$V2==1 & PositiveFrame2$V3==0, ])
  Positive.Class3[i] <- nrow(PositiveFrame2[PositiveFrame2$V2==0 & PositiveFrame2$V3==1, ])
  Positive.Class4[i] <- nrow(PositiveFrame2[PositiveFrame2$V2==1 & PositiveFrame2$V3==1, ])
  ClassImb.Positive[i] <- rowMeans(PositiveFrame)
  PositiveDataPhy = phyDat(PositiveFrame, type="USER", levels=c("0","1"))
  PositiveTree <- read.nexus(Positive.Trees[i])
  CI.Positive[i] <- CI(PositiveTree, PositiveDataPhy)
}

# Calculate CI - Negative Control
Negative.Files <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Control/RJMCMC_NegativeData/", pattern =".txt")
Negative.Trees <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Control/RJMCMC_NegativeData/", pattern =".nex")
setwd("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Control/RJMCMC_NegativeData/")
for(i in 1:length(Negative.Files)){
  NegativeData <- read.delim(Negative.Files[i], header=FALSE, row.names=1)
  tNegativeData <- t(NegativeData)
  NegativeFrame <- as.data.frame(tNegativeData, row.names = c("1","2"))
  NegativeFrame2 <- as.data.frame(NegativeData)
  Negative.Class1[i] <- nrow(NegativeFrame2[NegativeFrame2$V2==0 & NegativeFrame2$V3==0, ])
  Negative.Class2[i] <- nrow(NegativeFrame2[NegativeFrame2$V2==1 & NegativeFrame2$V3==0, ])
  Negative.Class3[i] <- nrow(NegativeFrame2[NegativeFrame2$V2==0 & NegativeFrame2$V3==1, ])
  Negative.Class4[i] <- nrow(NegativeFrame2[NegativeFrame2$V2==1 & NegativeFrame2$V3==1, ])
  ClassImb.Negative[i] <- rowMeans(NegativeFrame)
  NegativeDataPhy = phyDat(NegativeFrame, type="USER", levels=c("0","1"))
  NegativeTree <- read.nexus(Negative.Trees[i])
  CI.Negative[i] <- CI(NegativeTree, NegativeDataPhy)
}

# Calculate CI - Darwin's Scenario
Darwin.Files <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Experimental/RJMCMC_DarwinData/", pattern =".txt")
Darwin.Trees <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Experimental/RJMCMC_DarwinData/", pattern =".nex")
setwd("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Experimental/RJMCMC_DarwinData/")
for(i in 1:length(Darwin.Files)){
  DarwinData <- read.delim(Darwin.Files[i], header=FALSE, row.names=1)
  tDarwinData <- t(DarwinData)
  DarwinFrame <- as.data.frame(tDarwinData, row.names = c("1","2"))
  DarwinFrame2 <- as.data.frame(DarwinData)
  Darwin.Class1[i] <- nrow(DarwinFrame2[DarwinFrame2$V2==0 & DarwinFrame2$V3==0, ])
  Darwin.Class2[i] <- nrow(DarwinFrame2[DarwinFrame2$V2==1 & DarwinFrame2$V3==0, ])
  Darwin.Class3[i] <- nrow(DarwinFrame2[DarwinFrame2$V2==0 & DarwinFrame2$V3==1, ])
  Darwin.Class4[i] <- nrow(DarwinFrame2[DarwinFrame2$V2==1 & DarwinFrame2$V3==1, ])
  ClassImb.Darwin[i] <- rowMeans(DarwinFrame)
  DarwinDataPhy = phyDat(DarwinFrame, type="USER", levels=c("0","1"))
  DarwinTree <- read.nexus(Darwin.Trees[i])
  CI.Darwin[i] <- CI(DarwinTree, DarwinDataPhy)
}

# Calculate CI - Unreplicated Burst Scenario
Unrepburst.Files <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Experimental/RJMCMC_UnrepburstData/", pattern =".txt")
Unrepburst.Trees <- dir("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Experimental/RJMCMC_UnrepburstData/", pattern =".nex")
setwd("~/Desktop/Simulation Supplementary Material/Code/RJMCMC/RJMCMC Experimental/RJMCMC_UnrepburstData/")
for(i in 1:length(Unrepburst.Files)){
  UnrepburstData <- read.delim(Unrepburst.Files[i], header=FALSE, row.names=1)
  tUnrepburstData <- t(UnrepburstData)
  UnrepburstFrame <- as.data.frame(tUnrepburstData, row.names = c("1","2"))
  UnrepburstFrame2 <- as.data.frame(UnrepburstData)
  Unrepburst.Class1[i] <- nrow(UnrepburstFrame2[UnrepburstFrame2$V2==0 & UnrepburstFrame2$V3==0, ])
  Unrepburst.Class2[i] <- nrow(UnrepburstFrame2[UnrepburstFrame2$V2==1 & UnrepburstFrame2$V3==0, ])
  Unrepburst.Class3[i] <- nrow(UnrepburstFrame2[UnrepburstFrame2$V2==0 & UnrepburstFrame2$V3==1, ])
  Unrepburst.Class4[i] <- nrow(UnrepburstFrame2[UnrepburstFrame2$V2==1 & UnrepburstFrame2$V3==1, ])
  ClassImb.Unrepburst[i] <- rowMeans(UnrepburstFrame)
  UnrepburstDataPhy = phyDat(UnrepburstFrame, type="USER", levels=c("0","1"))
  UnrepburstTree <- read.nexus(Unrepburst.Trees[i])
  CI.Unrepburst[i] <- CI(UnrepburstTree, UnrepburstDataPhy)
}

# Results
median(CI.Positive) # low CI
median(CI.Negative) # low CI
median(CI.Darwin) # CI = 1
median(CI.Unrepburst) # CI < 1

median(ClassImb.Positive) # low within-trait imbalance
median(ClassImb.Negative) # low imbalance
median(ClassImb.Darwin) # low imbalance
median(ClassImb.Unrepburst) # low imbalance

# median multi-class imablance (c1=53,c2=0,c3=0,c4=47)
median(Positive.Class1) 
median(Positive.Class2)
median(Positive.Class3)
median(Positive.Class4)

# median multi-class imablance (c1=28,c2=25,c3=25,c4=22)
median(Negative.Class1) 
median(Negative.Class2)
median(Negative.Class3)
median(Negative.Class4)

# median multi-class imablance (c1=53,c2=0,c3=0,c4=47)
median(Darwin.Class1) 
median(Darwin.Class2)
median(Darwin.Class3)
median(Darwin.Class4)

# median multi-class imablance (c1=52,c2=24,c3=0,c4=24)
median(Unrepburst.Class1)
median(Unrepburst.Class2)
median(Unrepburst.Class3)
median(Unrepburst.Class4)

