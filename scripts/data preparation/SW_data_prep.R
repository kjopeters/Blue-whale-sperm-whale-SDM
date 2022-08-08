#SDM preparation
#24.09.2021
#
#creating the data packages to run the model loops
#..............................................................................
#running thel loops manually
library(biomod2)

#..............................................................................
setwd("~/Documents/Work/Massey/Projects/NCGIS/SDM3/R_project/SB_SDM")

source("scripts/TSSCalculation.R")

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

#read in data
SW = read.csv("data/sw_data_incl_PsA.csv")
SW$depth = as.numeric(SW$depth)

#make frames for results (only need to do once)
nRow <- nrow(env_data)

Test_results <- as.data.frame(matrix(0, ncol = 20, nrow = 9, dimnames = list(c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT" ), NULL)))
Test_results_TSS <- as.data.frame(matrix(0, ncol = 20, nrow = 9, dimnames = list(c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT" ), NULL)))

Pred_results <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_had_26 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_had_45 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_had_85 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_IPSL_26 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_IPSL_45 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_IPSL_85 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_MIROC_26 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_MIROC_45 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_MIROC_85 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_MPI_26 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_MPI_45 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

Pred_results_MPI_85 <- array(0, c(nRow, 9, 20), dimnames = list(seq(1:nRow), c("SRE", "GLM", "GAM", "MARS", "FDA", "ANN", "RF", "BRT", "MAXENT"), seq(1:20)))

saveRDS(Test_results, file = "results/models/SW/Test_results.rds")
saveRDS(Test_results_TSS, file = "results/models/SW/Test_results_TSS.rds")
saveRDS(Pred_results, file = "results/models/SW/Pred_results.rds")

saveRDS(Pred_results_had_26, file = "results/models/SW/Pred_results_had_26.rds")
saveRDS(Pred_results_had_45, file = "results/models/SW/Pred_results_had_45.rds")
saveRDS(Pred_results_had_85, file = "results/models/SW/Pred_results_had_85.rds")
saveRDS(Pred_results_IPSL_26, file = "results/models/SW/Pred_results_IPSL_26.rds")
saveRDS(Pred_results_IPSL_45, file = "results/models/SW/Pred_results_IPSL_45.rds")
saveRDS(Pred_results_IPSL_85, file = "results/models/SW/Pred_results_IPSL_85.rds")
saveRDS(Pred_results_MIROC_26, file = "results/models/SW/Pred_results_MIROC_26.rds")
saveRDS(Pred_results_MIROC_45, file = "results/models/SW/Pred_results_MIROC_45.rds")
saveRDS(Pred_results_MIROC_85, file = "results/models/SW/Pred_results_MIROC_85.rds")
saveRDS(Pred_results_MPI_26, file = "results/models/SW/Pred_results_MPI_26.rds")
saveRDS(Pred_results_MPI_45, file = "results/models/SW/Pred_results_MPI_45.rds")
saveRDS(Pred_results_MPI_85, file = "results/models/SW/Pred_results_MPI_85.rds")


#-----------------make datasets------------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib1.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval1.csv", row.names = F)

#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib2.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval2.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib3.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval3.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib4.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval4.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib5.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval5.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib6.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval6.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib7.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval7.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib8.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval8.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib9.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval9.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib10.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval10.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib11.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval11.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib12.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval12.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib13.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval13.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib14.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval14.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib15.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval15.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib16.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval16.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib17.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval17.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib18.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval18.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib19.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval19.csv", row.names = F)
#--------------
a <- SampleMat2(ref = SW$PA,
                ratio = 0.7) # function from the biomod2 package
calib <- SW[a$calibration, ]
eval <- SW[a$evaluation, ]

write.csv(calib, file = "results/models/SW/calib20.csv", row.names = F)
write.csv(eval, file = "results/models/SW/eval20.csv", row.names = F)

