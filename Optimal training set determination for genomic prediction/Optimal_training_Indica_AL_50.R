rm(list=ls())
library(TSDFGS)
library(SKM)
library(BGLR) #load BLR
library(caret)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)

################Example 1##########################################################
################Example 1##########################################################
load("Indica_AL.RData", verbose = TRUE)
#load the data that comes with the package
dataset <- "Indica_AL"
ls()
GRM=Geno
Pheno=Pheno
rownames(Geno)=1:ncol(Geno)
colnames(Geno)=1:ncol(Geno)
dim(GRM)
head(Pheno)
head(Geno[1:5,1:5])
head(GRM[1:5,1:5])
dim(Pheno)
Name_traits=colnames(Pheno)[2:5]
Name_traits
Sample_Size <- 50
iterations_number <- 3000
burn_in <- 2000

set.seed(123)

Cand <- 1:nrow(Geno)
GVoverall_result <- optTrain(
  geno = Geno,
  cand = Cand,
  n.train = 50,
  method = "rScore"
)

Cand <- 1:nrow(Geno)
GVoverall_result <- optTrain(
  geno = Geno,
  cand = Cand,
  n.train = 50,
  method = "PEV"
)

Cand <- 1:nrow(Geno)
GVoverall_result <- optTrain(
  geno = Geno,
  cand = Cand,
  n.train = 50,
  method = "CD"
)

TRS_indices <- GVoverall_result$train
TRS_indices
Sample_Size=50
Results=data.frame()
for (t  in 1: length(Name_traits)){
#  t=1
Trait_t=Name_traits[t]
Y=Pheno[,Trait_t]
#colnames(GRM)=1:ncol(GRM)
#rownames(GRM)=1:ncol(GRM)
dim(GRM)
varP <- c(var(Y)) # lsmean is a vector of phenotypes
h2 <- 0.30	# Trait heritility
Candidates=colnames(GRM)
length(Candidates)
ETA <- list(Line=list(model='RKHS', K=GRM))
y <- Y

##########Parameters##############
##DesN denotes the required multiple initial experiments to use this is a scalar
##Candidates denotes the names of the candidate set (this is a vector)
##Amat denotes the pedigree or genomic relationship matrix
##criteria denotes the optimization criteria "A" or "D" A for A-optimality and D for D-optimality
##n denotes the size of the required optimal training, this is a scalar
##varP denotes the phenotypic variance of the trait under evaluation
##h2 denotes the heritability of the trait under evaluation
##contrast should be TRUE or FALSE TRUE when PEV is computed as the
##difference between Testing and Average training.
##Testing denotes the testing set you need to provide otherwise provide NULL
##n_iter denotes the number of iterations 

ETA=list(Line=list(model='RKHS',K=GRM))
y=Y



# rownames(Geno)=1:ncol(Geno)
# colnames(Geno)=1:ncol(Geno)
OT_RScore=optTrain(
  Geno,
  cand=1:ncol(Geno),
  n.train=Sample_Size,
  subpop = NULL,
  test = NULL,
  method = "rScore",
  min.iter = NULL
)

SS_RScore=OT_RScore$OPTtrain
Pos_optimal_trn=SS_RScore
y_ff=y
y_ff[-Pos_optimal_trn]=NA
##############Training the regression model with a random sample#########
model_f<-BGLR::BGLR(
  y = y_ff,
  ETA = ETA,
  response_type = "gaussian",
  nIter = iterations_number,
  burnIn = burn_in,
  verbose = FALSE
)

Observed2=y[-Pos_optimal_trn]
Predicted2=model_f$yHat[-Pos_optimal_trn]
COR_Conv2=cor(Observed2,Predicted2)
MSE_Conv2=mse(Observed2,Predicted2)
NRMSE_Conv2=nrmse(Observed2,Predicted2, type="mean")

#####rScore, PEV and CD.
OT_PEV=optTrain(
  Geno,
  cand=1:ncol(Geno),
  n.train=Sample_Size,
  subpop = NULL,
  test = NULL,
  method = "PEV",
  min.iter = NULL
)

SS_PEV=OT_PEV$OPTtrain

Pos_optimal_trn3=SS_PEV
y_ff=y
y_ff[-Pos_optimal_trn3]=NA
##############Training the regression model with a random sample#########
model_f<-BGLR::BGLR(
  y = y_ff,
  ETA = ETA,
  response_type = "gaussian",
  nIter = iterations_number,
  burnIn = burn_in,
  verbose = FALSE
)

Observed3=y[-Pos_optimal_trn3]
Predicted3=model_f$yHat[-Pos_optimal_trn3]
COR_Conv3=cor(Observed3,Predicted3)
MSE_Conv3=mse(Observed3,Predicted3)
NRMSE_Conv3=nrmse(Observed3,Predicted3, type="mean")


#####rScore, PEV and CD.
OT_CD=optTrain(
  Geno,
  cand=1:ncol(Geno),
  n.train=Sample_Size,
  subpop = NULL,
  test = NULL,
  method = "CD",
  min.iter = NULL
)

SS_CD=OT_CD$OPTtrain

Pos_optimal_trn4=SS_CD
y_ff=y
y_ff[-Pos_optimal_trn4]=NA
##############Training the regression model with a random sample#########
model_f<-BGLR::BGLR(
  y = y_ff,
  ETA = ETA,
  response_type = "gaussian",
  nIter = iterations_number,
  burnIn = burn_in,
  verbose = FALSE
)

Observed4=y[-Pos_optimal_trn4]
Predicted4=model_f$yHat[-Pos_optimal_trn4]
COR_Conv4=cor(Observed4,Predicted4)
MSE_Conv4=mse(Observed4,Predicted4)
NRMSE_Conv4=nrmse(Observed4,Predicted4, type="mean")

# Results=rbind(Results, data.frame(Dataset=dataset,Trait=Trait_t,COR_Conv=COR_Conv, COR_RScore=COR_Conv2,COR_PEV=COR_Conv3, COR_CD=COR_Conv4,MSE_Conv=MSE_Conv, MSE_RScore=MSE_Conv2,MSE_PEV=MSE_Conv3,MSE_CD=MSE_Conv4, NRMSE_Conv=NRMSE_Conv, NRMSE_RScore=NRMSE_Conv2, NRMSE_PEV=NRMSE_Conv3, NRMSE_CD=NRMSE_Conv4,Sample_Size=Sample_Size))
# Results
for(r in 1:10){
  
  set.seed(r)
  
  pos_trn=sample(1:length(Candidates),Sample_Size)
  pos_trn
  iterations_number=3000
  burn_in=2000
  y_f=y
  y_f[-pos_trn]=NA
  ##############Training the regression model with a random sample#########
  model_f<-BGLR::BGLR(
    y = y_f,
    ETA = ETA,
    response_type = "gaussian",
    nIter = iterations_number,
    burnIn = burn_in,
    verbose = FALSE
  )
  
  Observed=y[-pos_trn]
  Predicted=model_f$yHat[-pos_trn]
  COR_Conv_val <- cor(Observed, Predicted)
  MSE_Conv_val <- mse(Observed, Predicted)
  NRMSE_Conv_val <- nrmse(Observed, Predicted, type="mean")
  Results_trait <- data.frame(
    Dataset = dataset,
    Trait = Trait_t,
    Rep = r,
    Sample_Size = Sample_Size,
    
    COR_Conv = COR_Conv_val,
    MSE_Conv = MSE_Conv_val,
    NRMSE_Conv = NRMSE_Conv_val,
    
    COR_RScore = COR_Conv2,
    COR_PEV = COR_Conv3,
    COR_CD = COR_Conv4,
    
    MSE_RScore = MSE_Conv2,
    MSE_PEV = MSE_Conv3,
    MSE_CD = MSE_Conv4,
    
    NRMSE_RScore = NRMSE_Conv2,
    NRMSE_PEV = NRMSE_Conv3,
    NRMSE_CD = NRMSE_Conv4
  )
  
  Results <- rbind(Results, Results_trait)
}
}

Results
write.csv(Results,file=paste("Results_",dataset,"_Sample_Size_",Sample_Size,".csv",sep=""))
