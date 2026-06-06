rm(list=ls())

library(TSDFGS)
library(SKM)
library(BGLR)
library(caret)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)

################ Cargar datos ################
load("2021_Winter_Wheat_Whashington_USA_Carter_2024.RData", verbose = TRUE)

dataset <- "2021_Winter_Wheat_Whashington_USA_Carter_2024"
Geno <- Geno[, -1]
colnames(Geno) <- gsub("\\.", "-", colnames(Geno))
Pheno$GID <- Pheno$Line
Pheno$Loc <- Pheno$Env
common_ids = intersect(Pheno$GID, colnames(Geno))
Pheno = Pheno[Pheno$GID %in% common_ids,]
Geno = Geno[common_ids]
Geno_mat <- as.matrix(t(Geno))

#GRM <- Geno
#GRM <- as.matrix(GRM)
#rownames(GRM) = colnames(GRM)

Name_traits <- colnames(Pheno)[5:18]
set.seed(123)
Sample_Size <- 50
Results <- data.frame()


################ LOOP ################
for (Env_Name in unique(Pheno$Loc)) {
  

  Pheno_env <- Pheno[Pheno$Loc == Env_Name, ]
  

  if (nrow(Pheno_env) <= Sample_Size) {
    print(paste("Skipping", Env_Name, " - Not enough plants for Sample_Size"))
    next
  }
    
 
  ids_env <- as.character(Pheno_env$GID)
  Geno_env <- Geno_mat[ids_env, ]
  GRM_env <- Geno_env %*% t(Geno_env) / ncol(Geno_env) 
  
  for (t in 1:length(Name_traits)) {
    Trait_name <- Name_traits[t]
    Y <- Pheno_env[[Trait_name]]
    Candidates <- 1:nrow(Geno_env)
  
  
    ################ MODELO BASE (aleatorio) ################
    ETA <- list(Line = list(model = 'RKHS', K = GRM_env))
    
    pos_trn <- sample(Candidates, Sample_Size)
    
    y_f <- Y
    y_f[-pos_trn] <- NA
    
    model_f <- BGLR::BGLR(
      y = y_f,
      ETA = ETA,
      response_type = "gaussian",
      nIter = 3000,
      burnIn = 2000,
      verbose = FALSE
    )
    
    Observed <- Y[-pos_trn]
    Predicted <- model_f$yHat[-pos_trn]
    
    COR_Conv <- cor(Observed, Predicted, use = "complete.obs")
    MSE_Conv <- mse(Observed, Predicted)
    NRMSE_Conv <- nrmse(Observed, Predicted, type = "mean")
    
    ################ rScore ################
    OT_RScore <- optTrain(
      geno = Geno_env,
      cand = Candidates,
      n.train = Sample_Size,
      method = "rScore"
    )
    
    pos_opt <-  OT_RScore$OPTtrain
    pos_opt
    
    y_ff <- Y
    y_ff[-pos_opt] <- NA
    
    model_f <- BGLR::BGLR(
      y = y_ff,
      ETA = ETA,
      response_type = "gaussian",
      nIter = 3000,
      burnIn = 2000,
      verbose = FALSE
    )
    
    Observed2 <- Y[-pos_opt]
    Predicted2 <- model_f$yHat[-pos_opt]
    
    COR_RScore <- cor(Observed2, Predicted2, use = "complete.obs")
    MSE_RScore <- mse(Observed2, Predicted2)
    NRMSE_RScore <- nrmse(Observed2, Predicted2, type = "mean")
    
    ################ PEV ################
    OT_PEV <- optTrain(
      geno = Geno_env,
      cand = Candidates,
      n.train = Sample_Size,
      method = "PEV"
    )
    
    pos_opt <- OT_PEV$OPTtrain
    
    y_ff <- Y
    y_ff[-pos_opt] <- NA
    
    model_f <- BGLR::BGLR(
      y = y_ff,
      ETA = ETA,
      response_type = "gaussian",
      nIter = 3000,
      burnIn = 2000,
      verbose = FALSE
    )
    
    Observed3 <- Y[-pos_opt]
    Predicted3 <- model_f$yHat[-pos_opt]
    
    COR_PEV <- cor(Observed3, Predicted3, use = "complete.obs")
    MSE_PEV <- mse(Observed3, Predicted3)
    NRMSE_PEV <- nrmse(Observed3, Predicted3, type = "mean")
    
    ################ CD ################
    OT_CD <- optTrain(
      geno = Geno_env,
      cand = Candidates,
      n.train = Sample_Size,
      method = "CD"
    )
    
    pos_opt <- OT_CD$OPTtrain
    
    y_ff <- Y
    y_ff[-pos_opt] <- NA
    
    model_f <- BGLR::BGLR(
      y = y_ff,
      ETA = ETA,
      response_type = "gaussian",
      nIter = 3000,
      burnIn = 2000,
      verbose = FALSE
    )
    
    Observed4 <- Y[-pos_opt]
    Predicted4 <- model_f$yHat[-pos_opt]
    

    
    if(length(Observed4) > 1 && sd(Observed4, na.rm=T) > 0) {
      COR_CD <- cor(Observed4, Predicted4, use = "complete.obs")
    } else {
      COR_CD <- NA
    }
    
    COR_CD <- cor(Observed4, Predicted4, use = "complete.obs")
    MSE_CD <- mse(Observed4, Predicted4)
    NRMSE_CD <- nrmse(Observed4, Predicted4, type = "mean")
    
    
    ################ GUARDAR RESULTADOS ################
    Results <- rbind(Results, data.frame(
      Dataset = dataset,
      Environment = Env_Name,
      Trait = Trait_name,
      COR_Conv = COR_Conv,
      COR_RScore = COR_RScore,
      COR_PEV = COR_PEV,
      COR_CD = COR_CD,
      MSE_Conv = MSE_Conv,
      MSE_RScore = MSE_RScore,
      MSE_PEV = MSE_PEV,
      MSE_CD = MSE_CD,
      NRMSE_Conv = NRMSE_Conv,
      NRMSE_RScore = NRMSE_RScore,
      NRMSE_PEV = NRMSE_PEV,
      NRMSE_CD = NRMSE_CD,
      Sample_Size = Sample_Size
    ))
    
    print(paste("Done:", Env_Name, "-", Trait_name))
  }
}


################ RESULTADO FINAL ################
Results
write.csv(Results,file=paste("Results_",dataset,"_Sample_Size_",Sample_Size,".csv",sep=""))

