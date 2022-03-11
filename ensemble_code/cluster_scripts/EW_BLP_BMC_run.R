library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)

# only train on server for now
# info
args <- commandArgs(TRUE)
test_s <- as.numeric(args[1])
if(test_s == 17){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016")
  lname <- "EW_BLP17"
} else if (test_s == 18){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016","2016/2017")
  lname <- "EW_BLP18"
} else if (test_s == 19){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016","2016/2017","2017/2018")
  lname <- "EW_BLP19"
}
targets<-c("1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead")

##--------------------------------------------------------------------------------------------------------------------##
# set path
#path <- "./ch1dat/"
path <- "/Users/ahctun_woon/git/end-of-2019-2020/cdc-flusight-ensemble/calibration_work/BLPdata_app/"
##--------------------------------------------------------------------------------------------------------------------##
# get in scores and pits
pitByModel<-read.csv(paste0(path,"pit_modelcol.csv")) %>%
  dplyr::filter(calendar_week %in% c(43:53,1:18),
                season %in% train_season)
prepitByModel<-read.csv(paste0(path,"prepit_modelcol.csv")) %>%
  dplyr::filter(calendar_week %in% c(43:53,1:18),
                season %in% train_season)
# pitByModel<-read.csv(paste0(path,"calibration_work/BLPdata_app/pit_modelcol.csv")) %>%
#   dplyr::filter(calendar_week %in% c(43:53,1:18),
#                 season %in% train_season)
# prepitByModel<-read.csv(paste0(path,"calibration_work/BLPdata_app/prepit_modelcol.csv")) %>%
#   dplyr::filter(calendar_week %in% c(43:53,1:18),
#                 season %in% train_season)
##--------------------------------------------------------------------------------------------------------------------##
# info
componentModel_list <- colnames(pitByModel)[5:ncol(pitByModel)]
mod_num <- length(componentModel_list)

##--------------------------------------------------Functions to get params-------------------------------------------##
source("calibration_work/cluster_scripts/ch1_funcs.R")

##----------------------------------------------------- BLP -----------------------------------------------------------##

#k=1
EW_BLP_params <- lapply(1:4, function(x) cv_func(targets[x],
                                              NULL,
                                              "EW-BLP",
                                              K=1,
                                              M=mod_num,
                                              component_pits=pitByModel,
                                              component_prepits=prepitByModel))


# cv (later)
for (i in 2:5){
  for (j in 1:4){
    assign(paste0("EW_BMC",i,"_params_cv_target",j),
           lapply(train_season, 
                  function(x) cv_func(targets[j],
                                      x,
                                      "EW-BLP",
                                      K=i,
                                      M=mod_num,
                                      component_pits=pitByModel,
                                      component_prepits=prepitByModel)
           )
    )
  }
}
##---------------------------------------------- save outputs ---------------------------------------------------------##
assign(lname,list(EW_BLP_params=EW_BLP_params,
                  EW_BMC2_params_cv_target1=EW_BMC2_params_cv_target1,
                  EW_BMC3_params_cv_target1=EW_BMC3_params_cv_target1,
                  EW_BMC4_params_cv_target1=EW_BMC4_params_cv_target1,
                  EW_BMC5_params_cv_target1=EW_BMC5_params_cv_target1,
                  EW_BMC2_params_cv_target2=EW_BMC2_params_cv_target2,
                  EW_BMC3_params_cv_target2=EW_BMC3_params_cv_target2,
                  EW_BMC4_params_cv_target2=EW_BMC4_params_cv_target2,
                  EW_BMC5_params_cv_target2=EW_BMC5_params_cv_target2,
                  EW_BMC2_params_cv_target3=EW_BMC2_params_cv_target3,
                  EW_BMC3_params_cv_target3=EW_BMC3_params_cv_target3,
                  EW_BMC4_params_cv_target3=EW_BMC4_params_cv_target3,
                  EW_BMC5_params_cv_target3=EW_BMC5_params_cv_target3,
                  EW_BMC2_params_cv_target4=EW_BMC2_params_cv_target4,
                  EW_BMC3_params_cv_target4=EW_BMC3_params_cv_target4,
                  EW_BMC4_params_cv_target4=EW_BMC4_params_cv_target4,
                  EW_BMC5_params_cv_target4=EW_BMC5_params_cv_target4))
save(list=lname,file=paste0(lname,".rda"))