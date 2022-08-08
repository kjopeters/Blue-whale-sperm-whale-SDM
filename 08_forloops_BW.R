#SDM blue and sperm whale
#running the models in for-loops (blue whales)

#04.10.2021

#..............................................................................
#rm(list = ls())
library(biomod2)
library(marmap)
library(RColorBrewer)
library(raster)
library(dplyr)
library(usdm)
library(PresenceAbsence)
library(caTools)
library(ggplot2)
library(MASS)
library(nnet)
library(rpart)
library(earth)
library(mda)
library(randomForest)
library(gbm)
library(cowplot)
library(maxnet)
library(rJava)
library(ModelMetrics)
library(ROCR)
library(caret)
library(dismo)
if (is.element('package:gam', search())) detach('package:gam') ## make sure the gam package is not loaded to avoid conflicts
library(mgcv)
library(Hmisc)

#..............................................................................
#..............................................................................
setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM")

source("scripts/TSSCalculation.R")

#read in data
BW = read.csv("data/bw_data_incl_PsA.csv")
BW$depth = as.numeric(BW$depth)

env_data = readRDS(file="data/env_data_s.rds")

#import prepared files for results
Test_results = readRDS("results/models/BW/Test_results.rds")
Test_results_TSS = readRDS("results/models/BW/Test_results_TSS.rds")

#these are the subsets created by using a different  70% (calib) and 30% (eval) of the original data each time
calib1 = read.csv("results/models/BW/calib1.csv")
eval1 = read.csv("results/models/BW/eval1.csv")
calib2 = read.csv("results/models/BW/calib2.csv")
eval2 = read.csv("results/models/BW/eval2.csv")
calib3 = read.csv("results/models/BW/calib3.csv")
eval3 = read.csv("results/models/BW/eval3.csv")
calib4 = read.csv("results/models/BW/calib4.csv")
eval4 = read.csv("results/models/BW/eval4.csv")
calib5 = read.csv("results/models/BW/calib5.csv")
eval5 = read.csv("results/models/BW/eval5.csv")
calib6 = read.csv("results/models/BW/calib6.csv")
eval6 = read.csv("results/models/BW/eval6.csv")
calib7 = read.csv("results/models/BW/calib7.csv")
eval7 = read.csv("results/models/BW/eval7.csv")
calib8 = read.csv("results/models/BW/calib8.csv")
eval8 = read.csv("results/models/BW/eval8.csv")
calib9 = read.csv("results/models/BW/calib9.csv")
eval9 = read.csv("results/models/BW/eval9.csv")
calib10 = read.csv("results/models/BW/calib10.csv")
eval10 = read.csv("results/models/BW/eval10.csv")
calib11 = read.csv("results/models/BW/calib11.csv")
eval11 = read.csv("results/models/BW/eval11.csv")
calib12 = read.csv("results/models/BW/calib12.csv")
eval12 = read.csv("results/models/BW/eval12.csv")
calib13 = read.csv("results/models/BW/calib13.csv")
eval13 = read.csv("results/models/BW/eval13.csv")
calib14 = read.csv("results/models/BW/calib14.csv")
eval14 = read.csv("results/models/BW/eval14.csv")
calib15 = read.csv("results/models/BW/calib15.csv")
eval15 = read.csv("results/models/BW/eval15.csv")
calib16 = read.csv("results/models/BW/calib16.csv")
eval16 = read.csv("results/models/BW/eval16.csv")
calib17 = read.csv("results/models/BW/calib17.csv")
eval17 = read.csv("results/models/BW/eval17.csv")
calib18 = read.csv("results/models/BW/calib18.csv")
eval18 = read.csv("results/models/BW/eval18.csv")
calib19 = read.csv("results/models/BW/calib19.csv")
eval19 = read.csv("results/models/BW/eval19.csv")
calib20 = read.csv("results/models/BW/calib20.csv")
eval20 = read.csv("results/models/BW/eval20.csv")

env_data = readRDS(file="data/env_data_s.rds")
FP_had_26 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP26.rds")
FP_had_45 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP45.rds")
FP_had_85 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP85.rds")
FP_IPSL_26 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP26.rds")
FP_IPSL_45 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP45.rds")
FP_IPSL_85 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP85.rds")
FP_MIROC_26 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP26.rds")
FP_MIROC_45 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP45.rds")
FP_MIROC_85 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP85.rds")
FP_MPI_26 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP26.rds")
FP_MPI_45 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP45.rds")
FP_MPI_85 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP85.rds")



calib <- list(calib1, calib2, calib3, calib4, calib5, calib6, calib7, calib8, calib9, calib10, calib11, calib12, calib13, calib14, calib15, calib16, calib17, calib18, calib19, calib20)

eval <- list(eval1, eval2, eval3, eval4, eval5, eval6, eval7, eval8, eval9, eval10, eval11, eval12, eval13, eval14, eval15, eval16, eval17, eval18, eval19, eval20)

env <- list(env_data, FP_had_26, FP_had_45, FP_had_85, FP_IPSL_26, FP_IPSL_45, FP_IPSL_85, FP_MIROC_26, FP_MIROC_45, FP_MIROC_85, FP_MPI_26, FP_MPI_45, FP_MPI_85)

pred_results <- list("Pred_results_env_data", "Pred_results_had_26", "Pred_results_had_45", "Pred_results_had_85", "Pred_results_IPSL_26", "Pred_results_IPSL_45","Pred_results_IPSL_85", "Pred_results_MIROC_26","Pred_results_MIROC_45", "Pred_results_MIROC_85", "Pred_results_MPI_26", "Pred_results_MPI_45", "Pred_results_MPI_85")

#*********************************************************************
#*********************************************************************
#
#           Run models
#
#*#*********************************************************************
#*#*********************************************************************

#########################################################
#   SRE 
#########################################################

ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  Pred_test  = sre(Response = calib[[i]]$PA, Explanatory = calib[[i]][,c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], NewData = eval[[i]][,c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], Quant = 0.025); Test_results['SRE', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_SRE_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]])
  Test_results_TSS['SRE', i] <- TSS_SRE_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'SRE', i]  = sre(Response = calib[[i]]$PA, Explanatory = calib[[i]][,c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], NewData = env[[j]][,c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], Quant = 0)
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
  
}
saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")



#########################################################
#   GLM 
#########################################################

ni = 20
nj = 13

glm.formula = formula("PA ~ 1 + poly(depth, 2) + poly(slope, 2) + poly(aspect, 2) + poly(distance_shore, 2) + poly(distance_200m, 2) + poly(distance_1000m, 2) + poly(chll, 2) + poly(sst, 2) + depth:slope + depth:aspect + depth:distance_shore + depth:distance_200m + depth:distance_1000m + depth:chll + depth:sst + slope:aspect + slope:distance_shore + slope:distance_200m + slope:distance_1000m + slope:chll + slope:sst + aspect:distance_shore + aspect:distance_200m + aspect:distance_1000m + aspect:chll + aspect:sst + distance_shore:distance_200m + distance_shore:distance_1000m + distance_shore:chll + distance_shore:sst + distance_200m:distance_1000m + distance_200m:chll + distance_200m:sst + distance_1000m:chll + distance_1000m:sst + chll:sst")

for (i in 1:ni){
  print((i/ni)*100)
  glmStart = glm(PA ~ 1, data = calib[[i]], family = "binomial")
 # glmModBIC = stepAIC(glmStart, glm.formula, data = calib[[i]], direction = "both", trace = F, k=2, control=glm.control(maxit = 100))
  glmModBIC = stepAIC(glmStart, glm.formula, direction = "both", trace = F, k = log(nrow(calib[[i]])), control = glm.control(maxit = 100))
  # prediction on the evaluation data and evaluation using the # AUC approach
  Pred_test <- predict.glm(glmModBIC, eval[[i]], type = 'response'); Test_results['GLM', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_GLM_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]])
  Test_results_TSS['GLM', i] <- TSS_GLM_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'GLM', i] <- predict(glmModBIC, env[[j]], type = 'response')
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}

saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")


#########################################################
#   GAM 
#########################################################
ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  gam_mgcv <- gam(PA ~ s(depth) + s(slope) + s(aspect) + s(distance_shore)+ s(distance_200m) + s(distance_1000m) +s(chll) + s(sst), data = calib[[i]], family = 'binomial')
  # prediction on the evaluation data and evaluation using the # AUC approach
  Pred_test <- predict.gam(gam_mgcv, newdata = eval[[i]], type = "response"); Test_results['GAM', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_GAM_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]])
  Test_results_TSS['GAM', i] <- TSS_GAM_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'GAM', i] <- predict.gam(gam_mgcv, newdata = env[[j]], type = 'response')
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}

saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")




#########################################################
#    MARS 
#########################################################

ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  Mars_int1 <- earth(PA ~ 1 + depth + slope + aspect + aspect + distance_shore + distance_200m + distance_1000m + chll +sst, data = calib[[i]], degree = 1, glm = list(family = binomial))
  # prediction on the evaluation data and evaluation using the # AUC approach
  Pred_test <- predict(Mars_int1, eval[[i]], type = 'response'); Test_results['MARS', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_MARS_res = getTSS(Pred_test, eval[[i]]$PA, eval[[1]])
  Test_results_TSS['MARS', i] <- TSS_MARS_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'MARS', i] <- predict(Mars_int1, env[[j]], type = 'response')
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}

saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")


#########################################################
#    FDA 
#########################################################

ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  fda_mod <- fda(PA ~ 1 + depth + slope + aspect + aspect + distance_shore + distance_200m + distance_1000m + chll + sst, data = calib[[i]], method = mars)
  # prediction on the evaluation data and evaluation using the # AUC approach
  Pred_test <- predict(fda_mod, eval[[i]], type = 'posterior')[, 2]; Test_results['FDA', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_FDA_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]])
  Test_results_TSS['FDA', i] <- TSS_FDA_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'FDA', i] <- predict(fda_mod, env[[j]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")],type = 'posterior')[,2]
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}  



saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")

#########################################################
#    ANN 
#########################################################

ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  set.seed(555)
  CV_nnet <- biomod2:::.CV.nnet(Input = calib[[i]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")],Target = calib[[i]]$PA)
  nnet.Final <- nnet(calib[[i]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], calib[[i]]$PA, size = CV_nnet[1, 1],rang = 0.1, decay = CV_nnet[1, 2], maxit = 200, trace = F)
  # prediction on the evaluation data and evaluation using the # AUC approach
  pred.test <- predict(nnet.Final, eval[[i]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")]); Test_results['ANN', i] <- somers2(pred.test, eval[[i]]$PA)['C']
  #get TSS
  TSS_ANN_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]])
  Test_results_TSS['ANN', i] <- TSS_ANN_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'ANN', i] <- predict(nnet.Final, env[[j]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")])
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}  

saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")



#########################################################
#    RF 
#########################################################

ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  RF_mod <- randomForest(x = calib[[i]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], y = as.factor(calib[[i]]$PA),ntree = 1000, importance = TRUE)
  # prediction on the evaluation data and evaluation using the # AUC approach
  Pred_test <- predict(RF_mod, eval[[i]], type = 'prob')[, 2]; Test_results['RF', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_RF_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]]); Test_results_TSS['RF', i] <- TSS_RF_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'RF', i] <- predict(RF_mod, env[[j]], type = 'prob')[, 2]
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}  

saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")



#########################################################
#    BRT 
#########################################################
library(biomod2)
library(marmap)
library(RColorBrewer)
library(raster)
library(dplyr)
library(usdm)
library(PresenceAbsence)
library(caTools)
library(ggplot2)
library(MASS)
library(nnet)
library(rpart)
library(earth)
library(mda)
library(randomForest)
library(gbm)
library(cowplot)
library(maxnet)
library(rJava)
library(ModelMetrics)
library(ROCR)
library(caret)
library(dismo)
library(gam)
library(Hmisc)

#..............................................................................
#..............................................................................
setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM")

source("scripts/TSSCalculation.R")

#read in data
BW = read.csv("data/bw_data_incl_PsA.csv")
BW$depth = as.numeric(BW$depth)

env_data = readRDS(file="data/env_data_s.rds")

#import prepared files for results
Test_results = readRDS("results/models/BW/Test_results.rds")
Test_results_TSS = readRDS("results/models/BW/Test_results_TSS.rds")


calib1 = read.csv("results/models/BW/calib1.csv")
eval1 = read.csv("results/models/BW/eval1.csv")
calib2 = read.csv("results/models/BW/calib2.csv")
eval2 = read.csv("results/models/BW/eval2.csv")
calib3 = read.csv("results/models/BW/calib3.csv")
eval3 = read.csv("results/models/BW/eval3.csv")
calib4 = read.csv("results/models/BW/calib4.csv")
eval4 = read.csv("results/models/BW/eval4.csv")
calib5 = read.csv("results/models/BW/calib5.csv")
eval5 = read.csv("results/models/BW/eval5.csv")
calib6 = read.csv("results/models/BW/calib6.csv")
eval6 = read.csv("results/models/BW/eval6.csv")
calib7 = read.csv("results/models/BW/calib7.csv")
eval7 = read.csv("results/models/BW/eval7.csv")
calib8 = read.csv("results/models/BW/calib8.csv")
eval8 = read.csv("results/models/BW/eval8.csv")
calib9 = read.csv("results/models/BW/calib9.csv")
eval9 = read.csv("results/models/BW/eval9.csv")
calib10 = read.csv("results/models/BW/calib10.csv")
eval10 = read.csv("results/models/BW/eval10.csv")
calib11 = read.csv("results/models/BW/calib11.csv")
eval11 = read.csv("results/models/BW/eval11.csv")
calib12 = read.csv("results/models/BW/calib12.csv")
eval12 = read.csv("results/models/BW/eval12.csv")
calib13 = read.csv("results/models/BW/calib13.csv")
eval13 = read.csv("results/models/BW/eval13.csv")
calib14 = read.csv("results/models/BW/calib14.csv")
eval14 = read.csv("results/models/BW/eval14.csv")
calib15 = read.csv("results/models/BW/calib15.csv")
eval15 = read.csv("results/models/BW/eval15.csv")
calib16 = read.csv("results/models/BW/calib16.csv")
eval16 = read.csv("results/models/BW/eval16.csv")
calib17 = read.csv("results/models/BW/calib17.csv")
eval17 = read.csv("results/models/BW/eval17.csv")
calib18 = read.csv("results/models/BW/calib18.csv")
eval18 = read.csv("results/models/BW/eval18.csv")
calib19 = read.csv("results/models/BW/calib19.csv")
eval19 = read.csv("results/models/BW/eval19.csv")
calib20 = read.csv("results/models/BW/calib20.csv")
eval20 = read.csv("results/models/BW/eval20.csv")

env_data = readRDS(file="data/env_data_s.rds")
FP_had_26 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP26.rds")
FP_had_45 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP45.rds")
FP_had_85 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP85.rds")
FP_IPSL_26 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP26.rds")
FP_IPSL_45 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP45.rds")
FP_IPSL_85 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP85.rds")
FP_MIROC_26 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP26.rds")
FP_MIROC_45 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP45.rds")
FP_MIROC_85 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP85.rds")
FP_MPI_26 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP26.rds")
FP_MPI_45 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP45.rds")
FP_MPI_85 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP85.rds")


calib <- list(calib1, calib2, calib3, calib4, calib5, calib6, calib7, calib8, calib9, calib10, calib11, calib12, calib13, calib14, calib15, calib16, calib17, calib18, calib19, calib20)

eval <- list(eval1, eval2, eval3, eval4, eval5, eval6, eval7, eval8, eval9, eval10, eval11, eval12, eval13, eval14, eval15, eval16, eval17, eval18, eval19, eval20)

env <- list(env_data, FP_had_26, FP_had_45, FP_had_85, FP_IPSL_26, FP_IPSL_45, FP_IPSL_85, FP_MIROC_26, FP_MIROC_45, FP_MIROC_85, FP_MPI_26, FP_MPI_45, FP_MPI_85)


pred_results <- list("Pred_results_env_data", "Pred_results_had_26", "Pred_results_had_45", "Pred_results_had_85", "Pred_results_IPSL_26", "Pred_results_IPSL_45","Pred_results_IPSL_85", "Pred_results_MIROC_26","Pred_results_MIROC_45", "Pred_results_MIROC_85", "Pred_results_MPI_26", "Pred_results_MPI_45", "Pred_results_MPI_85")


ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  GBM.mod <- gbm(PA ~ depth + slope + aspect + aspect + distance_shore + distance_200m + distance_1000m + chll + sst, data = calib[[i]], distribution = 'bernoulli', n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,bag.fraction = 0.5, cv.folds = 10)
  gbm.mod.perf <- gbm.perf(GBM.mod, method = 'cv', plot.it = F)
  # prediction on the evaluation data and evaluation using the # AUC approach
  Pred_test <- predict(GBM.mod, newdata = eval[[i]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], type = 'response', n.trees = gbm.mod.perf); Test_results['BRT', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_BRT_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]]); Test_results_TSS['BRT', i] <- TSS_BRT_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    Pred_results[, 'BRT', i] <-  predict(GBM.mod, newdata = env[[j]][, c("depth", "slope", "aspect", "distance_shore", "distance_200m", "distance_1000m", "chll", "sst")], type = 'response', n.trees = gbm.mod.perf)
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}  

saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")



#########################################################
#    MaxEnt 
#########################################################
setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM")
options(java.parameters = "-Xmx1g" )

library(biomod2)
library(rJava)
library(dismo)
library(Hmisc)
library(raster)


source("scripts/TSSCalculation.R")

#read in data
BW = read.csv("data/bw_data_incl_PsA.csv")
BW$depth = as.numeric(BW$depth)



#import prepared files for results
Test_results = readRDS("results/models/BW/Test_results.rds")
Test_results_TSS = readRDS("results/models/BW/Test_results_TSS.rds")


calib1 = read.csv("results/models/BW/calib1.csv")
eval1 = read.csv("results/models/BW/eval1.csv")
calib2 = read.csv("results/models/BW/calib2.csv")
eval2 = read.csv("results/models/BW/eval2.csv")
calib3 = read.csv("results/models/BW/calib3.csv")
eval3 = read.csv("results/models/BW/eval3.csv")
calib4 = read.csv("results/models/BW/calib4.csv")
eval4 = read.csv("results/models/BW/eval4.csv")
calib5 = read.csv("results/models/BW/calib5.csv")
eval5 = read.csv("results/models/BW/eval5.csv")
calib6 = read.csv("results/models/BW/calib6.csv")
eval6 = read.csv("results/models/BW/eval6.csv")
calib7 = read.csv("results/models/BW/calib7.csv")
eval7 = read.csv("results/models/BW/eval7.csv")
calib8 = read.csv("results/models/BW/calib8.csv")
eval8 = read.csv("results/models/BW/eval8.csv")
calib9 = read.csv("results/models/BW/calib9.csv")
eval9 = read.csv("results/models/BW/eval9.csv")
calib10 = read.csv("results/models/BW/calib10.csv")
eval10 = read.csv("results/models/BW/eval10.csv")
calib11 = read.csv("results/models/BW/calib11.csv")
eval11 = read.csv("results/models/BW/eval11.csv")
calib12 = read.csv("results/models/BW/calib12.csv")
eval12 = read.csv("results/models/BW/eval12.csv")
calib13 = read.csv("results/models/BW/calib13.csv")
eval13 = read.csv("results/models/BW/eval13.csv")
calib14 = read.csv("results/models/BW/calib14.csv")
eval14 = read.csv("results/models/BW/eval14.csv")
calib15 = read.csv("results/models/BW/calib15.csv")
eval15 = read.csv("results/models/BW/eval15.csv")
calib16 = read.csv("results/models/BW/calib16.csv")
eval16 = read.csv("results/models/BW/eval16.csv")
calib17 = read.csv("results/models/BW/calib17.csv")
eval17 = read.csv("results/models/BW/eval17.csv")
calib18 = read.csv("results/models/BW/calib18.csv")
eval18 = read.csv("results/models/BW/eval18.csv")
calib19 = read.csv("results/models/BW/calib19.csv")
eval19 = read.csv("results/models/BW/eval19.csv")
calib20 = read.csv("results/models/BW/calib20.csv")
eval20 = read.csv("results/models/BW/eval20.csv")



env_data = readRDS(file="data/env_data_s.rds")
env_depth = env_data[,c(1,2,3)]
env_slope = env_data[,c(1,2,4)]
env_aspect = env_data[,c(1,2,5)]
env_distance_shore = env_data[,c(1,2,6)]
env_distance_200m = env_data[,c(1,2,7)]
env_distance_1000m = env_data[,c(1,2,8)]
env_sst = env_data[,c(1,2,9)]
env_chll = env_data[,c(1,2,10)]

env_depth.r = rasterFromXYZ(env_depth); names(env_depth.r) = "depth"
env_slope.r = rasterFromXYZ(env_slope); names(env_slope.r) = "slope"
env_aspect.r = rasterFromXYZ(env_aspect); names(env_aspect.r) = "aspect"
env_distance_shore.r = rasterFromXYZ(env_distance_shore); names(env_distance_shore.r) = "distance_shore"
env_distance_200m.r = rasterFromXYZ(env_distance_200m); names(env_distance_200m.r) = "distance_200m"
env_distance_1000m.r = rasterFromXYZ(env_distance_1000m); names(env_distance_1000m.r) = "distance_1000m"
env_sst.r = rasterFromXYZ(env_sst); names(env_sst.r) = "sst"
env_chll.r = rasterFromXYZ(env_chll); names(env_chll.r) = "chll"

env_data_stack = raster::stack(env_depth.r, env_slope.r, env_aspect.r, env_distance_shore.r, env_distance_200m.r, env_distance_1000m.r, env_sst.r, env_chll.r)




had_26 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP26.rds")
had_26_depth = had_26[,c(1,2,3)]
had_26_slope = had_26[,c(1,2,4)]
had_26_aspect = had_26[,c(1,2,5)]
had_26_distance_shore = had_26[,c(1,2,6)]
had_26_distance_200m = had_26[,c(1,2,7)]
had_26_distance_1000m = had_26[,c(1,2,8)]
had_26_sst = had_26[,c(1,2,9)]
had_26_chll = had_26[,c(1,2,10)]

had_26_depth.r = rasterFromXYZ(had_26_depth); names(had_26_depth.r) = "depth"
had_26_slope.r = rasterFromXYZ(had_26_slope); names(had_26_slope.r) = "slope"
had_26_aspect.r = rasterFromXYZ(had_26_aspect); names(had_26_aspect.r) = "aspect"
had_26_distance_shore.r = rasterFromXYZ(had_26_distance_shore); names(had_26_distance_shore.r) = "distance_shore"
had_26_distance_200m.r = rasterFromXYZ(had_26_distance_200m); names(had_26_distance_200m.r) = "distance_200m"
had_26_distance_1000m.r = rasterFromXYZ(had_26_distance_1000m); names(had_26_distance_1000m.r) = "distance_1000m"
had_26_sst.r = rasterFromXYZ(had_26_sst); names(had_26_sst.r) = "sst"
had_26_chll.r = rasterFromXYZ(had_26_chll); names(had_26_chll.r) = "chll"

had_26_stack = raster::stack(had_26_depth.r, had_26_slope.r, had_26_aspect.r, had_26_distance_shore.r, had_26_distance_200m.r, had_26_distance_1000m.r, had_26_sst.r, had_26_chll.r)


had_45 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP45.rds")
had_45_depth = had_45[,c(1,2,3)]
had_45_slope = had_45[,c(1,2,4)]
had_45_aspect = had_45[,c(1,2,5)]
had_45_distance_shore = had_45[,c(1,2,6)]
had_45_distance_200m = had_45[,c(1,2,7)]
had_45_distance_1000m = had_45[,c(1,2,8)]
had_45_sst = had_45[,c(1,2,9)]
had_45_chll = had_45[,c(1,2,10)]

had_45_depth.r = rasterFromXYZ(had_45_depth); names(had_45_depth.r) = "depth"
had_45_slope.r = rasterFromXYZ(had_45_slope); names(had_45_slope.r) = "slope"
had_45_aspect.r = rasterFromXYZ(had_45_aspect); names(had_45_aspect.r) = "aspect"
had_45_distance_shore.r = rasterFromXYZ(had_45_distance_shore); names(had_45_distance_shore.r) = "distance_shore"
had_45_distance_200m.r = rasterFromXYZ(had_45_distance_200m); names(had_45_distance_200m.r) = "distance_200m"
had_45_distance_1000m.r = rasterFromXYZ(had_45_distance_1000m); names(had_45_distance_1000m.r) = "distance_1000m"
had_45_sst.r = rasterFromXYZ(had_45_sst); names(had_45_sst.r) = "sst"
had_45_chll.r = rasterFromXYZ(had_45_chll); names(had_45_chll.r) = "chll"

had_45_stack = raster::stack(had_45_depth.r, had_45_slope.r, had_45_aspect.r, had_45_distance_shore.r, had_45_distance_200m.r, had_45_distance_1000m.r, had_45_sst.r, had_45_chll.r)

had_85 = readRDS(file="data/climate_data_future/HadGEM2-ES_2100_RCP85.rds")
had_85_depth = had_85[,c(1,2,3)]
had_85_slope = had_85[,c(1,2,4)]
had_85_aspect = had_85[,c(1,2,5)]
had_85_distance_shore = had_85[,c(1,2,6)]
had_85_distance_200m = had_85[,c(1,2,7)]
had_85_distance_1000m = had_85[,c(1,2,8)]
had_85_sst = had_85[,c(1,2,9)]
had_85_chll = had_85[,c(1,2,10)]

had_85_depth.r = rasterFromXYZ(had_85_depth); names(had_85_depth.r) = "depth"
had_85_slope.r = rasterFromXYZ(had_85_slope); names(had_85_slope.r) = "slope"
had_85_aspect.r = rasterFromXYZ(had_85_aspect); names(had_85_aspect.r) = "aspect"
had_85_distance_shore.r = rasterFromXYZ(had_85_distance_shore); names(had_85_distance_shore.r) = "distance_shore"
had_85_distance_200m.r = rasterFromXYZ(had_85_distance_200m); names(had_85_distance_200m.r) = "distance_200m"
had_85_distance_1000m.r = rasterFromXYZ(had_85_distance_1000m); names(had_85_distance_1000m.r) = "distance_1000m"
had_85_sst.r = rasterFromXYZ(had_85_sst); names(had_85_sst.r) = "sst"
had_85_chll.r = rasterFromXYZ(had_85_chll); names(had_85_chll.r) = "chll"

had_85_stack = raster::stack(had_85_depth.r, had_85_slope.r, had_85_aspect.r, had_85_distance_shore.r, had_85_distance_200m.r, had_85_distance_1000m.r, had_85_sst.r, had_85_chll.r)

IPSL_26 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP26.rds")
IPSL_26_depth = IPSL_26[,c(1,2,3)]
IPSL_26_slope = IPSL_26[,c(1,2,4)]
IPSL_26_aspect = IPSL_26[,c(1,2,5)]
IPSL_26_distance_shore = IPSL_26[,c(1,2,6)]
IPSL_26_distance_200m = IPSL_26[,c(1,2,7)]
IPSL_26_distance_1000m = IPSL_26[,c(1,2,8)]
IPSL_26_sst = IPSL_26[,c(1,2,9)]
IPSL_26_chll = IPSL_26[,c(1,2,10)]

IPSL_26_depth.r = rasterFromXYZ(IPSL_26_depth); names(IPSL_26_depth.r) = "depth"
IPSL_26_slope.r = rasterFromXYZ(IPSL_26_slope); names(IPSL_26_slope.r) = "slope"
IPSL_26_aspect.r = rasterFromXYZ(IPSL_26_aspect); names(IPSL_26_aspect.r) = "aspect"
IPSL_26_distance_shore.r = rasterFromXYZ(IPSL_26_distance_shore); names(IPSL_26_distance_shore.r) = "distance_shore"
IPSL_26_distance_200m.r = rasterFromXYZ(IPSL_26_distance_200m); names(IPSL_26_distance_200m.r) = "distance_200m"
IPSL_26_distance_1000m.r = rasterFromXYZ(IPSL_26_distance_1000m); names(IPSL_26_distance_1000m.r) = "distance_1000m"
IPSL_26_sst.r = rasterFromXYZ(IPSL_26_sst); names(IPSL_26_sst.r) = "sst"
IPSL_26_chll.r = rasterFromXYZ(IPSL_26_chll); names(IPSL_26_chll.r) = "chll"

IPSL_26_stack = raster::stack(IPSL_26_depth.r, IPSL_26_slope.r, IPSL_26_aspect.r, IPSL_26_distance_shore.r, IPSL_26_distance_200m.r, IPSL_26_distance_1000m.r, IPSL_26_sst.r, IPSL_26_chll.r)

IPSL_45 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP45.rds")
IPSL_45_depth = IPSL_45[,c(1,2,3)]
IPSL_45_slope = IPSL_45[,c(1,2,4)]
IPSL_45_aspect = IPSL_45[,c(1,2,5)]
IPSL_45_distance_shore = IPSL_45[,c(1,2,6)]
IPSL_45_distance_200m = IPSL_45[,c(1,2,7)]
IPSL_45_distance_1000m = IPSL_45[,c(1,2,8)]
IPSL_45_sst = IPSL_45[,c(1,2,9)]
IPSL_45_chll = IPSL_45[,c(1,2,10)]

IPSL_45_depth.r = rasterFromXYZ(IPSL_45_depth); names(IPSL_45_depth.r) = "depth"
IPSL_45_slope.r = rasterFromXYZ(IPSL_45_slope); names(IPSL_45_slope.r) = "slope"
IPSL_45_aspect.r = rasterFromXYZ(IPSL_45_aspect); names(IPSL_45_aspect.r) = "aspect"
IPSL_45_distance_shore.r = rasterFromXYZ(IPSL_45_distance_shore); names(IPSL_45_distance_shore.r) = "distance_shore"
IPSL_45_distance_200m.r = rasterFromXYZ(IPSL_45_distance_200m); names(IPSL_45_distance_200m.r) = "distance_200m"
IPSL_45_distance_1000m.r = rasterFromXYZ(IPSL_45_distance_1000m); names(IPSL_45_distance_1000m.r) = "distance_1000m"
IPSL_45_sst.r = rasterFromXYZ(IPSL_45_sst); names(IPSL_45_sst.r) = "sst"
IPSL_45_chll.r = rasterFromXYZ(IPSL_45_chll); names(IPSL_45_chll.r) = "chll"

IPSL_45_stack = raster::stack(IPSL_45_depth.r, IPSL_45_slope.r, IPSL_45_aspect.r, IPSL_45_distance_shore.r, IPSL_45_distance_200m.r, IPSL_45_distance_1000m.r, IPSL_45_sst.r, IPSL_45_chll.r)

IPSL_85 = readRDS(file="data/climate_data_future/IPSL-CM5A-LR_2100_RCP85.rds")
IPSL_85_depth = IPSL_85[,c(1,2,3)]
IPSL_85_slope = IPSL_85[,c(1,2,4)]
IPSL_85_aspect = IPSL_85[,c(1,2,5)]
IPSL_85_distance_shore = IPSL_85[,c(1,2,6)]
IPSL_85_distance_200m = IPSL_85[,c(1,2,7)]
IPSL_85_distance_1000m = IPSL_85[,c(1,2,8)]
IPSL_85_sst = IPSL_85[,c(1,2,9)]
IPSL_85_chll = IPSL_85[,c(1,2,10)]

IPSL_85_depth.r = rasterFromXYZ(IPSL_85_depth); names(IPSL_85_depth.r) = "depth"
IPSL_85_slope.r = rasterFromXYZ(IPSL_85_slope); names(IPSL_85_slope.r) = "slope"
IPSL_85_aspect.r = rasterFromXYZ(IPSL_85_aspect); names(IPSL_85_aspect.r) = "aspect"
IPSL_85_distance_shore.r = rasterFromXYZ(IPSL_85_distance_shore); names(IPSL_85_distance_shore.r) = "distance_shore"
IPSL_85_distance_200m.r = rasterFromXYZ(IPSL_85_distance_200m); names(IPSL_85_distance_200m.r) = "distance_200m"
IPSL_85_distance_1000m.r = rasterFromXYZ(IPSL_85_distance_1000m); names(IPSL_85_distance_1000m.r) = "distance_1000m"
IPSL_85_sst.r = rasterFromXYZ(IPSL_85_sst); names(IPSL_85_sst.r) = "sst"
IPSL_85_chll.r = rasterFromXYZ(IPSL_85_chll); names(IPSL_85_chll.r) = "chll"

IPSL_85_stack = raster::stack(IPSL_85_depth.r, IPSL_85_slope.r, IPSL_85_aspect.r, IPSL_85_distance_shore.r, IPSL_85_distance_200m.r, IPSL_85_distance_1000m.r, IPSL_85_sst.r, IPSL_85_chll.r)

MIROC_26 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP26.rds")
MIROC_26_depth = MIROC_26[,c(1,2,3)]
MIROC_26_slope = MIROC_26[,c(1,2,4)]
MIROC_26_aspect = MIROC_26[,c(1,2,5)]
MIROC_26_distance_shore = MIROC_26[,c(1,2,6)]
MIROC_26_distance_200m = MIROC_26[,c(1,2,7)]
MIROC_26_distance_1000m = MIROC_26[,c(1,2,8)]
MIROC_26_sst = MIROC_26[,c(1,2,9)]
MIROC_26_chll = MIROC_26[,c(1,2,10)]

MIROC_26_depth.r = rasterFromXYZ(MIROC_26_depth); names(MIROC_26_depth.r) = "depth"
MIROC_26_slope.r = rasterFromXYZ(MIROC_26_slope); names(MIROC_26_slope.r) = "slope"
MIROC_26_aspect.r = rasterFromXYZ(MIROC_26_aspect); names(MIROC_26_aspect.r) = "aspect"
MIROC_26_distance_shore.r = rasterFromXYZ(MIROC_26_distance_shore); names(MIROC_26_distance_shore.r) = "distance_shore"
MIROC_26_distance_200m.r = rasterFromXYZ(MIROC_26_distance_200m); names(MIROC_26_distance_200m.r) = "distance_200m"
MIROC_26_distance_1000m.r = rasterFromXYZ(MIROC_26_distance_1000m); names(MIROC_26_distance_1000m.r) = "distance_1000m"
MIROC_26_sst.r = rasterFromXYZ(MIROC_26_sst); names(MIROC_26_sst.r) = "sst"
MIROC_26_chll.r = rasterFromXYZ(MIROC_26_chll); names(MIROC_26_chll.r) = "chll"

MIROC_26_stack = raster::stack(MIROC_26_depth.r, MIROC_26_slope.r, MIROC_26_aspect.r, MIROC_26_distance_shore.r, MIROC_26_distance_200m.r, MIROC_26_distance_1000m.r, MIROC_26_sst.r, MIROC_26_chll.r)

MIROC_45 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP45.rds")
MIROC_45_depth = MIROC_45[,c(1,2,3)]
MIROC_45_slope = MIROC_45[,c(1,2,4)]
MIROC_45_aspect = MIROC_45[,c(1,2,5)]
MIROC_45_distance_shore = MIROC_45[,c(1,2,6)]
MIROC_45_distance_200m = MIROC_45[,c(1,2,7)]
MIROC_45_distance_1000m = MIROC_45[,c(1,2,8)]
MIROC_45_sst = MIROC_45[,c(1,2,9)]
MIROC_45_chll = MIROC_45[,c(1,2,10)]

MIROC_45_depth.r = rasterFromXYZ(MIROC_45_depth); names(MIROC_45_depth.r) = "depth"
MIROC_45_slope.r = rasterFromXYZ(MIROC_45_slope); names(MIROC_45_slope.r) = "slope"
MIROC_45_aspect.r = rasterFromXYZ(MIROC_45_aspect); names(MIROC_45_aspect.r) = "aspect"
MIROC_45_distance_shore.r = rasterFromXYZ(MIROC_45_distance_shore); names(MIROC_45_distance_shore.r) = "distance_shore"
MIROC_45_distance_200m.r = rasterFromXYZ(MIROC_45_distance_200m); names(MIROC_45_distance_200m.r) = "distance_200m"
MIROC_45_distance_1000m.r = rasterFromXYZ(MIROC_45_distance_1000m); names(MIROC_45_distance_1000m.r) = "distance_1000m"
MIROC_45_sst.r = rasterFromXYZ(MIROC_45_sst); names(MIROC_45_sst.r) = "sst"
MIROC_45_chll.r = rasterFromXYZ(MIROC_45_chll); names(MIROC_45_chll.r) = "chll"

MIROC_45_stack = raster::stack(MIROC_45_depth.r, MIROC_45_slope.r, MIROC_45_aspect.r, MIROC_45_distance_shore.r, MIROC_45_distance_200m.r, MIROC_45_distance_1000m.r, MIROC_45_sst.r, MIROC_45_chll.r)

MIROC_85 = readRDS(file="data/climate_data_future/MIROC-ESM-CHEM_2100_RCP85.rds")
MIROC_85_depth = MIROC_85[,c(1,2,3)]
MIROC_85_slope = MIROC_85[,c(1,2,4)]
MIROC_85_aspect = MIROC_85[,c(1,2,5)]
MIROC_85_distance_shore = MIROC_85[,c(1,2,6)]
MIROC_85_distance_200m = MIROC_85[,c(1,2,7)]
MIROC_85_distance_1000m = MIROC_85[,c(1,2,8)]
MIROC_85_sst = MIROC_85[,c(1,2,9)]
MIROC_85_chll = MIROC_85[,c(1,2,10)]

MIROC_85_depth.r = rasterFromXYZ(MIROC_85_depth); names(MIROC_85_depth.r) = "depth"
MIROC_85_slope.r = rasterFromXYZ(MIROC_85_slope); names(MIROC_85_slope.r) = "slope"
MIROC_85_aspect.r = rasterFromXYZ(MIROC_85_aspect); names(MIROC_85_aspect.r) = "aspect"
MIROC_85_distance_shore.r = rasterFromXYZ(MIROC_85_distance_shore); names(MIROC_85_distance_shore.r) = "distance_shore"
MIROC_85_distance_200m.r = rasterFromXYZ(MIROC_85_distance_200m); names(MIROC_85_distance_200m.r) = "distance_200m"
MIROC_85_distance_1000m.r = rasterFromXYZ(MIROC_85_distance_1000m); names(MIROC_85_distance_1000m.r) = "distance_1000m"
MIROC_85_sst.r = rasterFromXYZ(MIROC_85_sst); names(MIROC_85_sst.r) = "sst"
MIROC_85_chll.r = rasterFromXYZ(MIROC_85_chll); names(MIROC_85_chll.r) = "chll"

MIROC_85_stack = raster::stack(MIROC_85_depth.r, MIROC_85_slope.r, MIROC_85_aspect.r, MIROC_85_distance_shore.r, MIROC_85_distance_200m.r, MIROC_85_distance_1000m.r, MIROC_85_sst.r, MIROC_85_chll.r)

MPI_26 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP26.rds")
MPI_26_depth = MPI_26[,c(1,2,3)]
MPI_26_slope = MPI_26[,c(1,2,4)]
MPI_26_aspect = MPI_26[,c(1,2,5)]
MPI_26_distance_shore = MPI_26[,c(1,2,6)]
MPI_26_distance_200m = MPI_26[,c(1,2,7)]
MPI_26_distance_1000m = MPI_26[,c(1,2,8)]
MPI_26_sst = MPI_26[,c(1,2,9)]
MPI_26_chll = MPI_26[,c(1,2,10)]

MPI_26_depth.r = rasterFromXYZ(MPI_26_depth); names(MPI_26_depth.r) = "depth"
MPI_26_slope.r = rasterFromXYZ(MPI_26_slope); names(MPI_26_slope.r) = "slope"
MPI_26_aspect.r = rasterFromXYZ(MPI_26_aspect); names(MPI_26_aspect.r) = "aspect"
MPI_26_distance_shore.r = rasterFromXYZ(MPI_26_distance_shore); names(MPI_26_distance_shore.r) = "distance_shore"
MPI_26_distance_200m.r = rasterFromXYZ(MPI_26_distance_200m); names(MPI_26_distance_200m.r) = "distance_200m"
MPI_26_distance_1000m.r = rasterFromXYZ(MPI_26_distance_1000m); names(MPI_26_distance_1000m.r) = "distance_1000m"
MPI_26_sst.r = rasterFromXYZ(MPI_26_sst); names(MPI_26_sst.r) = "sst"
MPI_26_chll.r = rasterFromXYZ(MPI_26_chll); names(MPI_26_chll.r) = "chll"

MPI_26_stack = raster::stack(MPI_26_depth.r, MPI_26_slope.r, MPI_26_aspect.r, MPI_26_distance_shore.r, MPI_26_distance_200m.r, MPI_26_distance_1000m.r, MPI_26_sst.r, MPI_26_chll.r)

MPI_45 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP45.rds")
MPI_45_depth = MPI_45[,c(1,2,3)]
MPI_45_slope = MPI_45[,c(1,2,4)]
MPI_45_aspect = MPI_45[,c(1,2,5)]
MPI_45_distance_shore = MPI_45[,c(1,2,6)]
MPI_45_distance_200m = MPI_45[,c(1,2,7)]
MPI_45_distance_1000m = MPI_45[,c(1,2,8)]
MPI_45_sst = MPI_45[,c(1,2,9)]
MPI_45_chll = MPI_45[,c(1,2,10)]

MPI_45_depth.r = rasterFromXYZ(MPI_45_depth); names(MPI_45_depth.r) = "depth"
MPI_45_slope.r = rasterFromXYZ(MPI_45_slope); names(MPI_45_slope.r) = "slope"
MPI_45_aspect.r = rasterFromXYZ(MPI_45_aspect); names(MPI_45_aspect.r) = "aspect"
MPI_45_distance_shore.r = rasterFromXYZ(MPI_45_distance_shore); names(MPI_45_distance_shore.r) = "distance_shore"
MPI_45_distance_200m.r = rasterFromXYZ(MPI_45_distance_200m); names(MPI_45_distance_200m.r) = "distance_200m"
MPI_45_distance_1000m.r = rasterFromXYZ(MPI_45_distance_1000m); names(MPI_45_distance_1000m.r) = "distance_1000m"
MPI_45_sst.r = rasterFromXYZ(MPI_45_sst); names(MPI_45_sst.r) = "sst"
MPI_45_chll.r = rasterFromXYZ(MPI_45_chll); names(MPI_45_chll.r) = "chll"

MPI_45_stack = raster::stack(MPI_45_depth.r, MPI_45_slope.r, MPI_45_aspect.r, MPI_45_distance_shore.r, MPI_45_distance_200m.r, MPI_45_distance_1000m.r, MPI_45_sst.r, MPI_45_chll.r)

MPI_85 = readRDS(file="data/climate_data_future/MPI-ESM-LR_2100_RCP85.rds")
MPI_85_depth = MPI_85[,c(1,2,3)]
MPI_85_slope = MPI_85[,c(1,2,4)]
MPI_85_aspect = MPI_85[,c(1,2,5)]
MPI_85_distance_shore = MPI_85[,c(1,2,6)]
MPI_85_distance_200m = MPI_85[,c(1,2,7)]
MPI_85_distance_1000m = MPI_85[,c(1,2,8)]
MPI_85_sst = MPI_85[,c(1,2,9)]
MPI_85_chll = MPI_85[,c(1,2,10)]

MPI_85_depth.r = rasterFromXYZ(MPI_85_depth); names(MPI_85_depth.r) = "depth"
MPI_85_slope.r = rasterFromXYZ(MPI_85_slope); names(MPI_85_slope.r) = "slope"
MPI_85_aspect.r = rasterFromXYZ(MPI_85_aspect); names(MPI_85_aspect.r) = "aspect"
MPI_85_distance_shore.r = rasterFromXYZ(MPI_85_distance_shore); names(MPI_85_distance_shore.r) = "distance_shore"
MPI_85_distance_200m.r = rasterFromXYZ(MPI_85_distance_200m); names(MPI_85_distance_200m.r) = "distance_200m"
MPI_85_distance_1000m.r = rasterFromXYZ(MPI_85_distance_1000m); names(MPI_85_distance_1000m.r) = "distance_1000m"
MPI_85_sst.r = rasterFromXYZ(MPI_85_sst); names(MPI_85_sst.r) = "sst"
MPI_85_chll.r = rasterFromXYZ(MPI_85_chll); names(MPI_85_chll.r) = "chll"

MPI_85_stack = raster::stack(MPI_85_depth.r, MPI_85_slope.r, MPI_85_aspect.r, MPI_85_distance_shore.r, MPI_85_distance_200m.r, MPI_85_distance_1000m.r, MPI_85_sst.r, MPI_85_chll.r)



calib <- list(calib1, calib2, calib3, calib4, calib5, calib6, calib7, calib8, calib9, calib10, calib11, calib12, calib13, calib14, calib15, calib16, calib17, calib18, calib19, calib20)

eval <- list(eval1, eval2, eval3, eval4, eval5, eval6, eval7, eval8, eval9, eval10, eval11, eval12, eval13, eval14, eval15, eval16, eval17, eval18, eval19, eval20)

env <- list(env_data_stack, had_26_stack, had_45_stack, had_85_stack, IPSL_26_stack, IPSL_45_stack, IPSL_85_stack, MIROC_26_stack, MIROC_45_stack, MIROC_85_stack, MPI_26_stack, MPI_45_stack, MPI_85_stack)


pred_results <- list("Pred_results_env_data", "Pred_results_had_26", "Pred_results_had_45", "Pred_results_had_85", "Pred_results_IPSL_26", "Pred_results_IPSL_45","Pred_results_IPSL_85", "Pred_results_MIROC_26","Pred_results_MIROC_45", "Pred_results_MIROC_85", "Pred_results_MPI_26", "Pred_results_MPI_45", "Pred_results_MPI_85")





Sys.setenv(NOAWT=TRUE) 


ni = 20
nj = 13

for (i in 1:ni){
  print((i/ni)*100)
  train_data_ME = calib[[i]][,c(4,5,6,7,8,9,10,11)] #get only predictor variables
  train_PA_ME = calib[[i]][,c(12)] # get only PAs
  test_ME = eval[[i]][,c(4,5,6,7,8,9,10,11)] #get only predictor variables
  maxent.model.train = maxent(train_data_ME, train_PA_ME)
  #predict on test dataset
  Pred_test = dismo::predict(maxent.model.train, test_ME, progress = "text"); Test_results['MAXENT', i] <- somers2(Pred_test, eval[[i]]$PA)['C']
  #get TSS
  TSS_ME_res = getTSS(Pred_test, eval[[i]]$PA, eval[[i]]); Test_results_TSS['MAXENT', i] <- TSS_ME_res[1,1]
  
  # prediction on the total dataset
  for (j in 1: nj){
    print(pred_results[[j]])
    setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM/results/models/BW")
    Pred_results <- readRDS(paste0(pred_results[[j]],".rds"))
    
    pred.maxent = dismo::predict(maxent.model.train, env[[j]], progress = "text")
    pred.maxent.df = rasterToPoints(pred.maxent)
    Pred_results[, 'MAXENT', i] = pred.maxent.df[,3]
    
    saveRDS(Pred_results, file= paste0(pred_results[[j]],".rds"))
    rm(Pred_results)
  }
}  


saveRDS(Test_results, file = "Test_results.rds")
saveRDS(Test_results_TSS, file = "Test_results_TSS.rds")

