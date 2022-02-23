library(tidyverse);library(cdcfluview);library(Matrix)
library(data.table);library(reshape2)
library(ranger);library(xtable)
library(here);library(stats);library(mixtools)
library(ggforce);library(grid);library(FluSight)
library(rmutil);library(R.utils)

##------------------------------------------- What does this do? -------------------------------------------------##

# get existing results for FSNetwork-TW and EW and calculate PITs

##--------------------------------------------- Specifications ----------------------------------------------------##
source("/Users/ahctun_woon/git/end-of-2019-2020/cdc-flusight-ensemble/calibration_work/BLPwork_functions_app/functions_app.R")

targets<-c("1 wk ahead","2 wk ahead","3 wk ahead","4 wk ahead")
model_list<-c("CU_EAKFC_SEIRS","CU_EAKFC_SIRS","CU_EKF_SEIRS","CU_EKF_SIRS","CU_RHF_SEIRS","CU_RHF_SIRS","CUBMA",
              "Delphi_BasisRegression","Delphi_EmpiricalFutures","Delphi_EmpiricalTrajectories", "Delphi_ExtendedDeltaDensity",
              "Delphi_MarkovianDeltaDensity", "Delphi_Uniform", "LANL_DBMplus", "Protea_Cheetah", "Protea_Kudu",
              "Protea_Springbok", "ReichLab_kcde", "ReichLab_kde", "ReichLab_sarima_seasonal_difference_FALSE",
              "ReichLab_sarima_seasonal_difference_TRUE", "equal-weights","target-based-weights")
componentModel_list<-c("CUBMA","CU_EAKFC_SEIRS","CU_EAKFC_SIRS","CU_EKF_SEIRS","CU_EKF_SIRS","CU_RHF_SEIRS",
                       "CU_RHF_SIRS","Delphi_BasisRegression","Delphi_EmpiricalFutures","Delphi_EmpiricalTrajectories",
                       "Delphi_ExtendedDeltaDensity","Delphi_MarkovianDeltaDensity","Delphi_Uniform","LANL_DBMplus",
                       "Protea_Cheetah","Protea_Kudu","Protea_Springbok",
                       "ReichLab_kcde_backfill_none","ReichLab_kcde_backfill_post_hoc","ReichLab_kde",
                       "ReichLab_sarima_seasonal_difference_FALSE","ReichLab_sarima_seasonal_difference_TRUE",
                       "FluOutlook_Mech","FluOutlook_MechAug",
                       "FluX_ARLR","FluX_LSTM","UA_EpiCos")
componentModel_list<-model_list[-c(22:23)]
model_list_sc<-c("CU-BMA","CU-EAKFC_SEIRS","CU-EAKFC_SIRS", "CU-EKF_SEIRS", "CU-EKF_SIRS", "CU-RHF_SEIRS", "CU-RHF_SIRS",
                 "Delphi-BasisRegression", "Delphi-EmpiricalFuture", "Delphi-EmpiricalTraj", "Delphi-DeltaDensity1",
                 "Delphi-DeltaDensity2", "Delphi-Uniform", "LANL-DBMplus", "Protea Analytics-Cheetah", "Protea Analytics-Kudu",
                 "Protea Analytics-Springbok", "ReichLab-KCDE", "ReichLab-KDE", "ReichLab-SARIMA1", "ReichLab-SARIMA2",
                 "FSNetwork-CW", "FSNetwork-EW", "FSNetwork-TRW", "FSNetwork-TW", "FSNetwork-TTW")
test_seasons<-c("2016/2017","2017/2018","2018/2019")
# get truths
truths<-read.csv('/Users/ahctun_woon/git/end-of-2018-2019/cdc-flusight-ensemble/scores/target-multivals_1819.csv') %>%
  dplyr::filter(Target %in% targets, Calendar.Week %in% c(43:53,1:18))
truths$Valid.Bin_start_incl<-as.numeric(as.character(truths$Valid.Bin_start_incl))
# set non-BMC ensemble names
nonBMC_ensembles <- c("TLP","EW","BLP","EW-BLP")

##--------------------------------------------- RUN ----------------------------------------------------##
# loop through each ensemble
for (ensemble in nonBMC_ensembles){
  # loop through each test season
  # set initial values (for now - calibration params is 1,1 (uniform transformation), weights are equal)
  inits <- c(1,1,rep(1/length(componentModel_list),length(componentModel_list)))
  if (ensemble=="TLP"){
    inits <- inits[-c(1:2)]
    test_seasons1 <- test_seasons
  } else if (ensemble=="BLP"){
    inits <- inits
    test_seasons1 <- test_seasons[1:2]
  } else if (ensemble=="EW"){
    inits <- NULL
    test_seasons1 <- test_seasons
  } else if (ensemble=="EW-BLP"){
    inits <- inits[1:2]
    test_seasons1 <- test_seasons[1:2]
  }
  # loop through each season
  for (ts in test_seasons1) {
    short_ts <- substr(ts, nchar(ts)-1, nchar(ts))
    # save results dynamically
    assign(paste0(ensemble,"_",short_ts),
           make_freq_ensemble(ensemble_name=ensemble, 
                              target_list=targets, 
                              single_test_season=ts, 
                              truths=truths, 
                              initial_vals=inits))
  }
}
# start here!!!!!!
# separate run for BLP and EW-BLP for 18/19 with different initial parameters
# still doesn't work, use stan is an option here
for (ensemble in nonBMC_ensembles[4]){
  # loop through each test season
  # set initial values (use tlp weights and change thecalibration)
  w <- TLP_18[[1]]$params
  inits <- c(1,1,w)
  if (ensemble=="BLP"){
    inits <- inits
  }else if (ensemble=="EW-BLP"){
    inits <- inits[1:2]
  }
  # loop through each season
  for (ts in test_seasons[3]) {
    short_ts <- substr(ts, nchar(ts)-1, nchar(ts))
    # save results dynamically
    assign(paste0(ensemble,"_",short_ts),
           make_freq_ensemble(ensemble_name=ensemble, 
                              target_list=targets, 
                              single_test_season=ts, 
                              truths=truths, 
                              initial_vals=inits))
  }
}


# optim TLP and EM flusight are mean log score are close to  the first/second decimal point (with in 0.05 margin)
# fw_sc <- read.csv("/Users/ahctun_woon/git/end-of-2018-2019/cdc-flusight-ensemble/BLPdata_app/scores_rename.csv") %>%
#          dplyr::filter(Model=="target-based-weights")
# fw_sc2 <- fw_sc %>%
#   dplyr::group_by(Target,Season) %>%
#   dplyr::mutate(mls=mean(Score)) %>%
#   ungroup()
# tlp_optimw <- data.frame()
# for (ts in test_seasons) {
#   for (tar in 1:length(targets)){
#     tlp_opw <- data.frame(weights=matrix(get(paste0("TLP_",substr(ts, nchar(ts)-1, nchar(ts))))[[tar]][[1]],ncol=1))
#     tlp_opw$component_model_id <- model_list_sc[1:21]
#     tlp_opw$season <- ts
#     tlp_opw$target <- targets[tar]
#     rbind(tlp_optimw,tlp_opw) -> tlp_optimw
#   } 
# }
# 
# tw_weights <- target_based_weights %>%
#   dplyr::filter(season %in% test_seasons,target %in% targets) %>%
#   dplyr::left_join(tlp_optimw, by=c("season","component_model_id","target")) %>%
#   dplyr::mutate(diff = round(weight-weights,4))

# check ew scores
# ew_sc <- read.csv("/Users/ahctun_woon/git/end-of-2018-2019/cdc-flusight-ensemble/BLPdata_app/scores_rename.csv") %>%
#   dplyr::filter(Model=="equal-weights", Season %in% test_seasons)
# 
# ew_sc2 <- ew_sc %>%
#   dplyr::group_by(Target,Season) %>%
#   dplyr::mutate(mls=mean(Score)) %>%
#   ungroup()
##------------------------------------------- Save Results  --------------------------------------------------##
# set main path
main_path <- "/Users/ahctun_woon/git/end-of-2018-2019/cdc-flusight-ensemble/BLPdata_app/ensemble_results/"
# build a function to collapse results??

for (ensemble in nonBMC_ensembles){
  if(ensemble=="EW-BLP") {
    ens_path <- paste0(main_path,"EW_BLP/")
  }
  else {
    ens_path <- paste0(main_path,ensemble,"/")
  }
  combine_ens <- NULL
  for (ts in test_seasons[1:2]) {
    short_ts <- substr(ts, nchar(ts)-1, nchar(ts))
    combine_ts <- NULL
    for (tar in targets){
      tar_num <- substr(tar,1,1)
      sp_path <- paste0(ens_path,short_ts,"/t",tar_num)
      # write files dynamically
      tar_results <- get(paste0(ensemble,"_",short_ts))[[as.numeric(tar_num)]]
      # combine params, mean train ls, and mean test ls
      combine_res <- data.frame(cbind(matrix(tar_results$params,nrow=1),
                                      tar_results$mean_train_ls,
                                      tar_results$mean_test_ls))
      names(combine_res)[ncol(combine_res)-1] <- "mean_train_ls"
      names(combine_res)[ncol(combine_res)] <- "mean_test_ls"
      combine_res$target <- tar
      combine_res$test_season <- ts
      write.csv(tar_results$test_epdf, file=paste0(sp_path,"/test_epdf.csv"), quote = FALSE, row.names = FALSE)
      write.csv(tar_results$pit_ls_tarloc, file=paste0(sp_path,"/test_pit_ls_info.csv"), quote = FALSE, row.names = FALSE)
      rbind(combine_ts,combine_res) -> combine_ts
    }
    rbind(combine_ens,combine_ts) -> combine_ens
    write.csv(combine_ens, file=paste0(ens_path,"params_mean_ls.csv"), quote = FALSE, row.names = FALSE)
  }
}
