rm(list = ls())
library(SKM)
library(plyr)
library(tidyr)
library(dplyr)
library(BGLR)

#####Loading data set#############
dataset_file="DMario_Family_14"
##load(paste0(dataset_file,".RData"))

Geno <- read.csv("C:/Users/Lenovo/Desktop/Osval/Unimodal/DMARIO/Geno_DMario_Family_14.csv")
Pheno <- read.csv("C:/Users/Lenovo/Desktop/Osval/Unimodal/DMARIO/Pheno_DMario_Family_14.csv")

Geno <- Geno[, -1]
ls()
Pheno=Pheno
Geno=Geno
head(Pheno)
#Pheno$GID=Pheno$Line
Pheno$Line = Pheno$GID

######Selecting the traits to be evaluated###############
Traits_to_evaluate=colnames(Pheno)[c(5:7)]

iterations_number <- 10000
burn_in <- 5000

##############Metric to compute percentage of matching###
best_lines_match <- function(Data, proportion = 0.1) {
  best_lines_number <- floor(nrow(Data) * proportion)
  
  best_lines <- Data %>%
    arrange(desc(Observed)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  predicted_lines <- Data %>%
    arrange(desc(Predicted)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  percentage_matching <- sum(predicted_lines %in% best_lines) /
    best_lines_number *
    100
  
  return(percentage_matching)
}

Pheno <- Pheno %>% arrange(Loc, GID)
Pheno <- droplevels(Pheno)
# Select the unique lines both in pheno and geno and sort them
rownames(Geno) <- colnames(Geno)
final_geno_lines <- intersect(Pheno$GID, rownames(Geno)) %>% sort()    ##change rownames(Geno) for colnames(Geno) and it works ####
Geno <- Geno[final_geno_lines, final_geno_lines]

###########Directory to save the results#####
results_dir <- file.path(
  "No_priors_in_G_and_E-mod ETA Lines",
  dataset_file)
mkdir(results_dir)

# Data preparation -------------------------------------------------------------
Predictions_Envs=data.frame()
Sumary_Envs=data.frame()

folds <- cv_random_line(Pheno$Line,folds_number=10, testing_proportion =0.50)  
folds
##########Sorting lines in Geno
Predictions_Traits=data.frame()
Sumary_Traits=data.frame()

ZL=model.matrix(~0+Line,data=Pheno)
ZE=model.matrix(~0+Loc,data=Pheno)
Geno=data.matrix(Geno)
K_E=ZE%*%t(ZE)
###Geno <- Geno+ 5 * diag(nrow(Geno))
L=t(chol(Geno))
K_G=ZL%*%Geno%*%t(ZL)
K_GE=K_E*K_G
XL=ZL%*%L
Predictions_Final_traits=data.frame()
Sumary_all_traits=data.frame()

for (t in 1:length(Traits_to_evaluate)){
  #  t=1
  Trait_t=Traits_to_evaluate[t]
  SKM::echo("\t*** Trait %s / %s ***", t, length(Traits_to_evaluate))
  #########Response variable in Obregon
  y <- Pheno[,Trait_t]
  y_f <- y
  
  y_f[testing_indices]=NA
  y_f
  Predictions_Final=data.frame()
  Sumary_all=data.frame()
  for(i in seq_along(folds )) {
    #i=1
    fold_i=folds[[i]]
    
    testing_indices <-fold_i$testing
    
    ##########ETA1 con Env y Lines ==P1############## uni-inputs
    ETA1=list(Lines=list(model="FIXED",X=XL))
    
    ##############Training with the whole oregon data set#############################
    model_f<-BGLR::BGLR(
      y = y_f,
      ETA = ETA1,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedA=model_f$yHat[testing_indices]
    Observed=y[testing_indices]
    ##########ETA2 con Env y Lines ==P1##############
    
    
    ETA2=list(Lines=list(model="BRR",X=XL))
    
    
    ##############Training with the whole oregon data set#############################
    model_f2<-BGLR::BGLR(
      y = y_f,
      ETA = ETA2,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedB=model_f2$yHat[testing_indices]
    
    
    ETA3=list(Lines=list(model="BayesA",X=XL))
    
    
    ##############Training with the whole oregon data set#############################
    model_f3<-BGLR::BGLR(
      y = y_f,
      ETA = ETA3,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedC=model_f3$yHat[testing_indices]
    
    
    ETA4=list(Lines=list(model="BayesB",X=XL))
    
    ##############Training with the whole oregon data set#############################
    model_f4<-BGLR::BGLR(
      y = y_f,
      ETA = ETA4,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedD=model_f4$yHat[testing_indices]
    
    
    ETA5=list(Lines=list(model="BayesC",X=XL))
    
    
    ##############Training with the whole oregon data set#############################
    model_f5<-BGLR::BGLR(
      y = y_f,
      ETA = ETA5,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedE=model_f5$yHat[testing_indices]
    
    ETA6=list(Lines=list(model="BL",X=XL))
    
    
    ##############Training with the whole oregon data set#############################
    model_f6<-BGLR::BGLR(
      y = y_f,
      ETA = ETA6,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedF=model_f6$yHat[testing_indices]
    
    
    ETA7=list(Lines=list(model="RKHS",K=K_G))
    
    ##############Training with the whole oregon data set#############################
    model_f7<-BGLR::BGLR(
      y = y_f,
      ETA = ETA7,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    PredictedG=model_f7$yHat[testing_indices]
    
    #U_Pred=data.frame(Line=Pheno$Line[fold$testing],Pred=g_Pred_testing1, Observed=g_True_testing1)
    Data_A=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedA)
    
    #####Metrics hole testing
    COR_A=cor(Observed,PredictedA)
    MSE_A=mse(Observed,PredictedA)
    NRMSE_A=nrmse(Observed,PredictedA, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_A_10=best_lines_match(Data=Data_A,proportion = 0.1) 
    PM_A_10
    PM_A_20=best_lines_match(Data=Data_A,proportion = 0.2) 
    PM_A_20
    PM_A_30=best_lines_match(Data=Data_A,proportion = 0.3) 
    PM_A_30
    
    #U_Pred=data.frame(Line=Pheno$Line[fold$testing],Pred=g_Pred_testing1, Observed=g_True_testing1)
    Data_B=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedB)
    
    #####Metrics hole testing
    COR_B=cor(Observed,PredictedB)
    MSE_B=mse(Observed,PredictedB)
    NRMSE_B=nrmse(Observed,PredictedB, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_B_10=best_lines_match(Data=Data_B,proportion = 0.1) 
    PM_B_10
    PM_B_20=best_lines_match(Data=Data_B,proportion = 0.2) 
    PM_B_20
    PM_B_30=best_lines_match(Data=Data_B,proportion = 0.3) 
    PM_B_30
    
    #U_Pred=data.frame(Line=Pheno$Line[fold$testing],Pred=g_Pred_testing1, Observed=g_True_testing1)
    Data_C=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedC)
    
    #####Metrics hole testing
    COR_C=cor(Observed,PredictedC)
    MSE_C=mse(Observed,PredictedC)
    NRMSE_C=nrmse(Observed,PredictedC, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_C_10=best_lines_match(Data=Data_C,proportion = 0.1) 
    PM_C_10
    PM_C_20=best_lines_match(Data=Data_C,proportion = 0.2) 
    PM_C_20
    PM_C_30=best_lines_match(Data=Data_C,proportion = 0.3) 
    PM_C_30
    
    Data_D=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedD)
    
    #####Metrics hole testing
    COR_D=cor(Observed,PredictedD)
    MSE_D=mse(Observed,PredictedD)
    NRMSE_D=nrmse(Observed,PredictedD, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_D_10=best_lines_match(Data=Data_D,proportion = 0.1) 
    PM_D_10
    PM_D_20=best_lines_match(Data=Data_D,proportion = 0.2) 
    PM_D_20
    PM_D_30=best_lines_match(Data=Data_D,proportion = 0.3) 
    PM_D_30
    
    Data_E=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedE)
    
    #####Metrics hole testing
    COR_E=cor(Observed,PredictedE)
    MSE_E=mse(Observed,PredictedE)
    NRMSE_E=nrmse(Observed,PredictedE, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_E_10=best_lines_match(Data=Data_E,proportion = 0.1) 
    PM_E_10
    PM_E_20=best_lines_match(Data=Data_E,proportion = 0.2) 
    PM_E_20
    PM_E_30=best_lines_match(Data=Data_E,proportion = 0.3) 
    PM_E_30
    
    Data_F=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedF)
    
    #####Metrics hole testing
    COR_F=cor(Observed,PredictedF)
    MSE_F=mse(Observed,PredictedF)
    NRMSE_F=nrmse(Observed,PredictedF, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_F_10=best_lines_match(Data=Data_F,proportion = 0.1) 
    PM_F_10
    PM_F_20=best_lines_match(Data=Data_F,proportion = 0.2) 
    PM_F_20
    PM_F_30=best_lines_match(Data=Data_F,proportion = 0.3) 
    PM_F_30
    
    Data_G=data.frame(Line=Pheno$Line[testing_indices],Observed=Observed, Predicted=PredictedG)
    
    #####Metrics hole testing
    COR_G=cor(Observed,PredictedG)
    MSE_G=mse(Observed,PredictedG)
    NRMSE_G=nrmse(Observed,PredictedG, type="mean")
    
    
    ###########Percentage of metrics####
    
    PM_G_10=best_lines_match(Data=Data_G,proportion = 0.1) 
    PM_G_10
    PM_G_20=best_lines_match(Data=Data_G,proportion = 0.2) 
    PM_G_20
    PM_G_30=best_lines_match(Data=Data_G,proportion = 0.3) 
    PM_G_30
    
    Sumary=data.frame(
      Dataset =dataset_file,
      Trait = Trait_t,
      Fold=i,
      COR_A=COR_A,
      COR_B=COR_B,
      COR_C=COR_C,
      COR_D=COR_D,
      COR_E=COR_E,
      COR_F=COR_F,
      COR_G=COR_G,
      MSE_A=MSE_A,
      MSE_B=MSE_B,
      MSE_C=MSE_C,
      MSE_D=MSE_D,
      MSE_E=MSE_E,
      MSE_F=MSE_F,
      MSE_G=MSE_G,
      NRMSE_A=NRMSE_A,
      NRMSE_B=NRMSE_B,
      NRMSE_C=NRMSE_C,
      NRMSE_D=NRMSE_D,
      NRMSE_E=NRMSE_E,
      NRMSE_F=NRMSE_F,
      NRMSE_G=NRMSE_G,
      PM_A_10=PM_A_10,
      PM_B_10=PM_B_10,
      PM_C_10=PM_C_10,
      PM_D_10=PM_D_10,
      PM_E_10=PM_E_10,
      PM_F_10=PM_F_10,
      PM_G_10=PM_G_10,
      PM_A_20=PM_A_20,
      PM_B_20=PM_B_20,
      PM_C_20=PM_C_20,
      PM_D_20=PM_D_20,
      PM_E_20=PM_E_20,
      PM_F_20=PM_F_20,
      PM_G_20=PM_G_20,
      PM_A_30=PM_A_30,
      PM_B_30=PM_B_30,
      PM_C_30=PM_C_30,
      PM_D_30=PM_D_30,
      PM_E_30=PM_E_30,
      PM_F_30=PM_F_30,
      PM_G_30=PM_G_30
    )
    Sumary
    Sumary_all=rbind(Sumary_all,Sumary)
    
    Predictions_i=data.frame(
      Dataset =dataset_file,
      Trait =Trait_t,
      Env =Pheno$Loc[testing_indices],
      Fold=i,
      Line = Pheno$Line[testing_indices],
      Observed =Observed,
      Predicted_A =PredictedA,
      Predicted_B =PredictedB,
      Predicted_C =PredictedC,
      Predicted_D =PredictedD,
      Predicted_E =PredictedE,
      Predicted_F =PredictedF,
      Predicted_G =PredictedG
    )
    Predictions_Final <-rbind(Predictions_Final,
                              Predictions_i
    )
  }
  Predictions_Final_traits=rbind(Predictions_Final_traits,Predictions_Final)
  Sumary_all_traits=rbind(Sumary_all_traits,Sumary_all)
}

###write.csv(Predictions_Final_traits, Predictions_Final_Nopriors.csv)
###write.csv(Sumary_all_traits , Sumary_all_traits_Nopriors.csv)
write_csv(Predictions_Final_traits, file.path(results_dir, "Predictions_Final_DMario_Family_14_Improved.csv"))
write_csv(Sumary_all_traits, file.path(results_dir, "Summary_all_traits_DMario_Family_14_Improved.csv"))
