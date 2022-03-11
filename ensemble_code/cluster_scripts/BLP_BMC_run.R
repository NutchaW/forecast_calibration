library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)

# only train on server for now
# info
args <- commandArgs(TRUE)
test_s <- as.numeric(args[1])
tar_num <- as.numeric(args[2])
b_num <- as.numeric(args[3])

if(test_s == 17){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016")
} else if (test_s == 18){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016","2016/2017")
} else if (test_s == 19){
  train_season<-c("2010/2011","2011/2012","2012/2013","2013/2014","2014/2015","2015/2016","2016/2017","2017/2018")
}
targets<-c("1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead")

##--------------------------------------------------------------------------------------------------------------------##
# set path
path <- "./ch1dat/"
#path <- "/Users/ahctun_woon/git/end-of-2019-2020/cdc-flusight-ensemble/"
##--------------------------------------------------------------------------------------------------------------------##
# get in scores and pits
pitByModel<-read.csv(paste0(path,"pit_modelcol.csv")) %>%
dplyr::filter(calendar_week %in% c(43:53,1:18),
              season %in% train_season,
              target == targets[tar_num])
prepitByModel<-read.csv(paste0(path,"prepit_modelcol.csv")) %>%
dplyr::filter(calendar_week %in% c(43:53,1:18),
              season %in% train_season,
              target == targets[tar_num])
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
source("ch1_funcs.R")
##----------------------------------------------------- BLP -----------------------------------------------------------##

if(b_num==1){
  #k=1
  lname <- paste0("BLP",test_s,"_target",tar_num)
  assign(lname,
         cv_func(targets[tar_num],
                 NULL,
                 "BLP",
                 K=1,
                 M=mod_num,
                 component_pits=pitByModel,
                 component_prepits=prepitByModel))
} else {
  lname <- paste0("BMC",b_num,test_s,"_target",tar_num,"_cv")
  assign(lname,
           lapply(train_season, 
                  function(x) cv_func(targets[tar_num],
                                      x,
                                      "BLP",
                                      K=b_num,
                                      M=mod_num,
                                      component_pits=pitByModel,
                                      component_prepits=prepitByModel)
           )
    )
}
##---------------------------------------------- save outputs ---------------------------------------------------------##

save(list=lname,file=paste0(lname,".rda"))