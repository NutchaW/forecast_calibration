# Presentation functions for markdown
cd_calc <- function(values){
  library(covidHubUtils)
  # get quantile
  # tau_F_vector <- tau_G_vector <- seq(0.025, 0.975, 0.025)
  # qf <- quantile(values, tau_F_vector)
  # qg <- qunif(tau_F_vector, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
  # approx_cd <- covidHubUtils::calc_cramers_dist_one_model_pair(qf, tau_F_vector, qg, tau_G_vector, "approximation1")
  # n <- length(values)
  # y <- runif(n)
  # a <- sum(c(sapply(values, function(x) abs(x-y))))/(n^2) 
  # b <- sum(c(dist(values)))/(n^2) 
  # c <- sum(c(dist(y)))/(n^2)  
  # approx_cd <- (a-(b/2)-(c/2))
  # integration
  c <- ecdf(values) 
  int_res <- integrate(function(x) ((c(x) - x)^2), 0, 1, subdivisions=2000)
  return(list(2,int_res[[1]]))
}
# cm <- cd_calc(pnorm(runif(10000)))
#cm2 <- cd_calc(punif(seq(0, 1, 0.01)))

# check <- lapply(some_files1, FUN=read.csv)
# check1 <- do.call(rbind.data.frame, check)
reli_cd <- function(file1,file2,test_seas,test=TRUE,cd_type=1){
  tmp <- lapply(file1, FUN=read.csv)
  frame1 <- do.call(rbind.data.frame, tmp)
  tmp2 <- lapply(file2, FUN=read.csv)
  frame2 <- do.call(rbind.data.frame, tmp2) 
  # %>%
  #   dplyr::mutate(model_name=recode_factor(model_name,`BLP`="LP"))
  frame <- rbind(frame1,frame2)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  if(test){
    frame <- frame %>%
      dplyr::filter(season==paste0("20",as.numeric(test_seas)-1,"/20",test_seas)) %>%
      dplyr::group_by(target,location,season,bin_start_incl,Method) %>%
      dplyr::mutate(rand_pit=sample(cdf_vals,1)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c("calendar_week","cdf_vals","value","ls")) %>%
      distinct() %>%
      dplyr::group_by(target,Method) %>%
      dplyr::mutate(cd= cd_calc(rand_pit)[[cd_type]]) %>%
      dplyr::ungroup()
    dis <- frame %>%
      dplyr::select(target,Method,cd) %>%
      distinct()
  } else {
    frame <- frame %>%
      dplyr::filter(season!=paste0("20",as.numeric(test_seas)-1,"/20",test_seas)) %>%
      dplyr::group_by(target,location,season,bin_start_incl,Method) %>%
      dplyr::mutate(rand_pit=sample(cdf_vals,1)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c("calendar_week","cdf_vals","value","ls")) %>%
      distinct() %>%
      dplyr::group_by(target,Method) %>%
      dplyr::mutate(cd= cd_calc(rand_pit)[[cd_type]]) %>%
      dplyr::ungroup()
    dis <- frame %>%
      dplyr::select(target,Method,cd) %>%
      distinct()
  }
  return(dis)
}

reli_cd2 <- function(file1,file2,file3,file4,file5,file6,test=TRUE,cd_type=1){
  tmp <- lapply(file1, FUN=read.csv)
  frame1 <- do.call(rbind.data.frame, tmp)
  tmp2 <- lapply(file2, FUN=read.csv)
  frame2 <- do.call(rbind.data.frame, tmp2) 
  tmp3 <- lapply(file3, FUN=read.csv)
  frame3 <- do.call(rbind.data.frame, tmp3)
  tmp4 <- lapply(file4, FUN=read.csv)
  frame4 <- do.call(rbind.data.frame, tmp4) 
  tmp5 <- lapply(file5, FUN=read.csv)
  frame5 <- do.call(rbind.data.frame, tmp5)
  tmp6 <- lapply(file6, FUN=read.csv)
  frame6 <- do.call(rbind.data.frame, tmp6) 
  # %>%
  #   dplyr::mutate(model_name=recode_factor(model_name,`BLP`="LP"))
  frame01 <- rbind(frame1,frame2)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  frame02 <- rbind(frame3,frame4)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  frame03 <- rbind(frame5,frame6)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  # have to fix these
  if(test){
    frame01 <- frame01 %>%
      dplyr::filter(season==paste0("20",as.numeric(17)-1,"/20",17))
    frame02 <- frame02 %>%
      dplyr::filter(season==paste0("20",as.numeric(18)-1,"/20",18)) 
    frame03 <- frame03 %>%
      dplyr::filter(season==paste0("20",as.numeric(19)-1,"/20",19)) 
  } else {
    frame01 <- frame01 %>%
      dplyr::filter(season!=paste0("20",as.numeric(17)-1,"/20",17))
    frame02 <- frame02 %>%
      dplyr::filter(season!=paste0("20",as.numeric(18)-1,"/20",18))
    frame03 <- frame03 %>%
      dplyr::filter(season!=paste0("20",as.numeric(19)-1,"/20",19))
  }
  dis <- rbind(frame01,frame02,frame03) %>%
    dplyr::group_by(target,location,season,bin_start_incl,Method) %>%
    dplyr::mutate(rand_pit=sample(cdf_vals,1)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(target,Method) %>%
    dplyr::mutate(cd= cd_calc(rand_pit)[[cd_type]]) %>%
    dplyr::ungroup() %>%
    dplyr::select(target,Method,cd) %>%
    distinct()
  return(dis)
}
####------------------------------------ Plot PIT function -------------------------------------------------#####
# work on making a diagonal plot
pitplot <- function(combined_dat,target_name,seasons, ens_name){
  dat_target<-combined_dat %>% 
    dplyr::filter(target==target_name,
                  season %in% seasons,
                  model_name==ens_name)
    ens_n <- ifelse(ens_name=="LP","FSNetwork-TW",
                    ifelse(ens_name=="EW","EW-LP",
                           ifelse(ens_name=="EW_BLP","EW-BLP",
                                  ifelse(ens_name=="BMC2",expression("BMC"[2]),
                                         ifelse(ens_name=="EW_BMC2",expression("EW-BMC"[2]),
                                                ifelse(ens_name=="EW_BMC5",expression("EW-BMC"[5]),
                                                       expression("EW-BMC"[3])))))))
    hist(dat_target$cdf_vals, breaks = 0:10/10, col = grey(0.5), freq = FALSE,
         ylim = c(0, 2), main = ens_n,xlab="",cex.main=0.8)
    abline( h = 1, lty = 2)
}

pitplot2 <- function(combined_dat,target_name,seasons, ens_name, ylimp){
  dat_target<-combined_dat %>% 
    dplyr::filter(target==target_name,
                  season %in% seasons,
                  model_name==ens_name)
  ens_n <- ifelse(ens_name=="LP",expression("FSNetwork-TW"),
                  ifelse(ens_name=="EW",expression("EW-LP"),
                         ifelse(ens_name=="EW_BLP",expression("EW-BLP"),
                                ifelse(ens_name=="BMC2",expression("BMC"[2]),
                                       ifelse(ens_name=="EW_BMC2",expression("EW-BMC"[2]),
                                              ifelse(ens_name=="EW_BMC5",expression("EW-BMC"[5]),
                                                     expression("EW-BMC"[3])))))))
  hist(dat_target$cdf_vals, breaks = 0:10/10, col = grey(0.5), freq = FALSE,
       ylim = ylimp, main = ens_n,xlab="",cex.main=0.8)
  abline( h = 1, lty = 2)
}
####------------------------------------ Create ls and params table-------------------------------------------------#####
# LStable_target<-function(combined_dat,target){
#   library(knitr)
#   avgDat<-combined_dat %>%
#     dplyr::filter(target==target) %>%
#     group_by(model_name) %>%
#     summarize(Average = mean(Score))
#   tab<-combined_dat %>%
#     dplyr::filter(Target==target) %>%
#     group_by(Model, Season) %>%
#     summarize(AvgLogScore = mean(Score)) %>%
#     pivot_wider(id_cols = Model, 
#                 names_from = Season, 
#                 values_from = AvgLogScore)
#   table<-tab %>% left_join(avgDat,by="Model")
#   return(kable(table,digits=3))
# }
# 
# LStable<-function(combined_dat){
#   library(knitr)
#   avgDat<-combined_dat %>%
#     group_by(Model) %>%
#     summarize(Average = mean(Score))
#   tab<-combined_dat %>%
#     group_by(Model, Season) %>%
#     summarize(AvgLogScore = mean(Score))%>%
#     pivot_wider(id_cols = Model,
#                 names_from = Season,
#                 values_from = AvgLogScore)
#   table<-tab %>% left_join(avgDat,by="Model") 
#   return(kable(table,digits=3))
# }

make_sim_table <- function(forecast_names,tar,tar_names,ts,test_seasons,path){
  param_file <- paste0("target",tar,"_sea",16+ts,"_params")
  load(paste0(path,"/train_params/final_train/",param_file,".rda"))
  params_bmc <- get(param_file)
  # arrange
  load(paste0(path,"/BLPdata_app/model_list.rda"))
  params_LP2 <- read.csv("~/git/end-of-2019-2020/cdc-flusight-ensemble/weights/target-based-weights.csv") %>%
    dplyr::filter(season==test_seasons[ts],
                  target==tar_names[tar]) 
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-BMA"] <- "CUBMA"
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-EAKFC_SEIRS"] <- "CU_EAKFC_SEIRS"
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-EAKFC_SIRS"] <- "CU_EAKFC_SIRS"
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-EKF_SEIRS"] <- "CU_EKF_SEIRS"
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-EKF_SIRS"] <- "CU_EKF_SIRS"
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-RHF_SEIRS"] <- "CU_RHF_SEIRS"
  params_LP2$component_model_id[params_LP2$component_model_id=="CU-RHF_SIRS"] <- "CU_RHF_SIRS"
  params_LP2$component_model_id[params_LP2$component_model_id=="Delphi-BasisRegression"] <- "Delphi_BasisRegression"
  params_LP2$component_model_id[params_LP2$component_model_id=="Delphi-EmpiricalFuture"] <- "Delphi_EmpiricalFutures"
  params_LP2$component_model_id[params_LP2$component_model_id=="Delphi-EmpiricalTraj"] <- "Delphi_EmpiricalTrajectories"
  params_LP2$component_model_id[params_LP2$component_model_id=="Delphi-DeltaDensity1"] <- "Delphi_ExtendedDeltaDensity"
  params_LP2$component_model_id[params_LP2$component_model_id=="Delphi-DeltaDensity2"] <- "Delphi_MarkovianDeltaDensity"
  params_LP2$component_model_id[params_LP2$component_model_id=="Delphi-Uniform"] <- "Delphi_Uniform"
  params_LP2$component_model_id[params_LP2$component_model_id=="LANL-DBMplus"] <- "LANL_DBMplus"
  params_LP2$component_model_id[params_LP2$component_model_id=="Protea-Cheetah"] <- "Protea_Cheetah"
  params_LP2$component_model_id[params_LP2$component_model_id=="Protea-Kudu"] <- "Protea_Kudu"
  params_LP2$component_model_id[params_LP2$component_model_id=="Protea-Springbok"] <- "Protea_Springbok"
  params_LP2$component_model_id[params_LP2$component_model_id=="ReichLab-KCDE_PHB"] <- "ReichLab_kcde_backfill_post_hoc"
  params_LP2$component_model_id[params_LP2$component_model_id=="ReichLab-KCDE_NB"] <- "ReichLab_kcde_backfill_none"
  params_LP2$component_model_id[params_LP2$component_model_id=="ReichLab-KDE"] <- "ReichLab_kde"
  params_LP2$component_model_id[params_LP2$component_model_id=="ReichLab-SARIMA1"] <- "ReichLab_sarima_seasonal_difference_FALSE"
  params_LP2$component_model_id[params_LP2$component_model_id=="ReichLab-SARIMA2"] <- "ReichLab_sarima_seasonal_difference_TRUE"
  params_LP2$component_model_id[params_LP2$component_model_id=="FSNetwork-TW"] <- "target-based-weights"
  params_LP2$component_model_id[params_LP2$component_model_id=="FSNetwork-EW"] <- "equal-weights"
  params_LP2$component_model_id[params_LP2$component_model_id=="FluOutlook-Mech"] <- "FluOutlook_Mech"
  params_LP2$component_model_id[params_LP2$component_model_id=="FluOutlook-MechAug"] <- "FluOutlook_MechAug"
  params_LP2$component_model_id[params_LP2$component_model_id=="FluX-ARLR"] <- "FluX_ARLR"
  params_LP2$component_model_id[params_LP2$component_model_id=="FluX-FluX_LSTM"] <- "FluX_LSTM"
  params_LP2$component_model_id[params_LP2$component_model_id=="UA-EpiCos"] <- "UA_EpiCos"
  params_LP <- params_LP2 %>%
    dplyr::arrange(factor(component_model_id, levels = model_list))
  # make table and kable it
  omega_list <- sapply(1:5, 
                       function(x) sapply(1:27, 
                                          function(y) paste0("$\\omega_{",x,y,"}$")))
  omega <- sapply(1:5, 
                       function(x) sapply(1:27, 
                                          function(y) paste0("omega_k[",x,y,"]")))
  param_names <- c("$w_1$", "$w_2$", "$w_3$","$w_4$","$w_5$", 
                   "$\\alpha_1$", "$\\beta_1$", "$\\alpha_2$", "$\\beta_2$", "$\\alpha_3$", "$\\beta_3$", 
                   "$\\alpha_4$", "$\\beta_4$","$\\alpha_5$", "$\\beta_5$",omega_list)
  param_cols <- c("w[1]", "w[2]", "w[3]","w[4]","w[5]", 
                  "alpha[1]", "beta[1]", "alpha[2]", "beta[2]", "alpha[3]", "beta[3]", 
                  "alpha[4]", "beta[4]","alpha[5]", "beta[5]",omega)
  table_params <- data.frame(matrix(NA, nrow = 6, ncol = length(param_cols)))
  colnames(table_params) <- param_cols
  table_params$Method <- forecast_names
  table_params <- table_params[c(ncol(table_params),1:(ncol(table_params)-1))]
  # params
  params1 <- round(1/27,3)
  params2 <- round(params_LP$weight,3)
  params3 <- round(params_bmc[[1]],3)
  params4 <- round(params_bmc[[3]],3)
  params5 <- round(params_bmc[[2]],3)
  params6 <- round(params_bmc[[4]],3)
    
  # fill tables
  non_mix <- forecast_names[forecast_names %in% c("BLP","LP","EW","EW_BLP")]
  for(ens in non_mix){
    for(i in 1:27){
      if(ens=="LP"){
        table_params[table_params$Method==ens, paste0("omega_k[",1,i,"]")] <- params2[i]
        } else if(ens=="EW"){
          table_params[table_params$Method==ens, paste0("omega_k[",1,i,"]")] <- params1
        } else if(ens=="BLP"){
          table_params[table_params$Method==ens, paste0("omega_k[",1,i,"]")] <- params3[i+2]
          table_params[table_params$Method==ens, paste0("alpha[",1,"]")] <- params3[1]
          table_params[table_params$Method==ens, paste0("beta[",1,"]")] <- params3[2]
        } else if(ens=="EW_BLP"){
          table_params[table_params$Method==ens, paste0("omega_k[",1,i,"]")] <- params1
          table_params[table_params$Method==ens, paste0("alpha[",1,"]")] <- params4[1]
          table_params[table_params$Method==ens, paste0("beta[",1,"]")] <- params4[2]
        }
      }
    }
    # fill again
  mix <- forecast_names[!(forecast_names %in% non_mix)]
  for(ens in mix){
    K <- as.numeric(substr(ens, nchar(ens), nchar(ens)))
    if(grepl("EW_BMC",ens,fixed = TRUE)){
      for(k in 1:K){
        table_params[table_params$Method==ens, paste0("alpha[",k,"]")] <- params6[k]
        table_params[table_params$Method==ens, paste0("beta[",k,"]")] <- params6[K+k]
        table_params[table_params$Method==ens, paste0("w[",k,"]")] <- params6[(2*K)+k]
        for(m in 1:27){
          table_params[table_params$Method==ens, paste0("omega_k[",k,m,"]")] <- params1
          }
      }
    } else {
      for(k in 1:K){
        table_params[table_params$Method==ens, paste0("alpha[",k,"]")] <- params5[k]
        table_params[table_params$Method==ens, paste0("beta[",k,"]")] <- params5[K+k]
        table_params[table_params$Method==ens, paste0("w[",k,"]")] <- params5[(2*K)+(K*27)+k]
        for(m in 1:27){
          table_params[table_params$Method==ens, paste0("omega_k[",k,m,"]")] <- params5[(2*K)+((m*K)-1)+(k-1)]
        }
      }
      }
    }
    
  # table_params[,2][table_params$Method=="EW_BMC1"|table_params$Method=="BMC1"] <- NA
  # table_params$Method[table_params$Method=="EW_BMC1"] <- "EW_BLP"
  # table_params$Method[table_params$Method=="BMC1"] <- "BLP"
  names(table_params)[2:ncol(table_params)] <- param_names
  table_params <- table_params[rowSums(is.na(table_params))<ncol(table_params),
                               colSums(is.na(table_params))<nrow(table_params)]
  table_params$Method <- gsub("_", "-", table_params$Method)
  return(table_params)
} 

#------------------------------------------------------new plot-------------------------------------------------------#

distance_heatmap <- function(cv_score,pname){
   name <- sapply(as.matrix(cv_score[,2]),
                   function(x)
                     ifelse(x=="BMC2",expression("BMC"[2]),
                            ifelse(x=="BMC3",expression("BMC"[3]),
                                   ifelse(x=="BMC4",expression("BMC"[4]),
                                          ifelse(x=="BMC5",expression("BMC"[5]),
                                                 ifelse(x=="EW_BMC2",expression("EW-BMC"[2]),
                                                        ifelse(x=="EW_BMC3",expression("EW-BMC"[3]),
                                                               ifelse(x=="EW_BMC4",expression("EW-BMC"[4]),
                                                                      expression("EW-BMC"[5])
                                                                      )
                                                               )
                                                        )
                                                 )
                                          )
                                   )
                            )
                   )
    cv_score %>%
      dplyr::mutate(group_ew=grepl("EW",`Model Name`)) %>%
      dplyr::group_by(Target,group_ew) %>%
      dplyr::mutate(max = max(round(`Mean validation log score`,2))
                    #,sds= sd(`Mean validation log score`),
                    ) %>%
      ungroup() %>%
      dplyr::mutate(
        #wthin= ifelse(`Mean validation log score` >= (min+sds)|`Mean validation log score` <= (min-sds),1,0)
        wthin= ifelse(grepl("2", `Model Name`),1,0)
                    ) %>%
    ggplot(aes(Target,`Model Name`)) +
      geom_tile(aes(fill= factor(rank)),show.legend = FALSE) +
      labs(title=pname)+
      scale_y_discrete(labels=name,limits = rev)+
      scale_x_discrete(expand = c(0, 0))+
      xlab("") +
      ylab("") +
      #scale_fill_distiller(palette = "PuBu",direction=-1) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(9,"PuBu")[1:4],
                        name="rank")+
      geom_text(aes(label = sprintf("%0.2f", round(`Mean validation log score`, 2)),
                    colour=factor(wthin)),size=3,show.legend = FALSE) +
      scale_colour_manual(values=c("black","red"))+
      theme_bw()+
      theme(plot.title = element_text(size=8),
            axis.text.x=element_text(size=rel(0.8)),
            axis.text.y=element_text(size=rel(0.8)),
            # legend.title = element_text(size=5),
            # legend.key.size = unit(0.3, 'cm'),
            # legend.position = "bottom",
            # legend.text = element_text(size=4),
            plot.margin=unit(c(0,0,0,0),"cm"),
            panel.background = element_blank(),
            plot.background = element_blank())
}


loc_heatmap <- function(loc_frame){
  lorder <- c("US National",rev(paste0("HHS Region ",1:10)))
  mod_order <- c("EW","TLP","EW_BLP","BLP","EW_BMC2","BMC2")
  name <- c(expression("EW-LP"),expression("LP"),expression("EW-BLP"),expression("BLP"),
            expression("EW-BMC"[2]),expression("BMC"[2]))
  loc <- loc_frame %>%
    dplyr::arrange(target,location,desc(loc_mean)) %>%
    dplyr::group_by(target,location) %>%
    dplyr::mutate(rank=rank(-round(loc_mean,2), ties.method="min")) %>%
    ungroup() 
  ggplot(loc,aes(factor(model_name,levels=mod_order),factor(location,levels=lorder))) +
    geom_tile(aes(fill= factor(rank))) +
    facet_wrap(~target) +
    xlab("") +
    ylab("") +
    #scale_fill_distiller(palette = "PuBu",direction = -1,trans = 'reverse') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9,"PuBu")[c(2:7)],
                      name="rank")+
    geom_text(aes(label = sprintf("%0.2f", round(loc_mean, 2))),size=3) +
    theme_bw()+
    theme(plot.title = element_text(size=8),
          legend.title = element_text(size=5),
          legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size=4),
          plot.margin=unit(c(0,0,0,0),"cm"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
          axis.text.x=element_text(size=rel(0.8)),
          axis.text.y=element_text(size=rel(0.8)))+
    scale_x_discrete(labels=name,expand = c(0, 0))+
    scale_y_discrete(expand = c(0, 0))
}

mean_heat <- function(yscore,pname){
  mod_order <- c("EW","TLP","EW_BLP","BLP","EW_BMC2","BMC2")
  name <- c(expression("EW-LP"),expression("LP"),expression("EW-BLP"),expression("BLP"),
            expression("EW-BMC"[2]),expression("BMC"[2]))
  yscore %>%
    ggplot(aes(factor(`Model Name`,levels=mod_order),test_train)) +
    geom_tile(aes(fill=factor(rank))) +
    theme(axis.text.x=element_text(size=rel(0.8)),
          axis.text.y=element_text(size=rel(0.8)))+
    facet_wrap(~Target,ncol=2) +
    #labs(title=pname)+
    xlab("") +
    ylab("") +
    ggtitle(pname) +
    # scale_fill_distiller(palette = "PuBu",direction=-1,trans = 'reverse') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9,"PuBu")[c(2:7)],
                      name="rank")+
    geom_text(aes(label = sprintf("%0.2f", round(value, 2))),size=3.5) +
    guides(fill=guide_legend(ncol=6)) +
    theme(plot.title = element_text(size=8),
          legend.title = element_text(size=5),
          legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size=4),
          legend.position = "bottom",
          plot.margin=unit(c(0.2,0.2,-0.4,0.2),"cm"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))+
    scale_x_discrete(labels=name,expand = c(0, 0))+
    scale_y_discrete(expand = c(0, 0))
}

# add sampling for each target, do within tar loc sea
reli_plot <- function(file1,file2,test_seas,test=TRUE){
  tmp <- lapply(file1, FUN=read.csv)
  frame1 <- do.call(rbind.data.frame, tmp)
  tmp2 <- lapply(file2, FUN=read.csv)
  frame2 <- do.call(rbind.data.frame, tmp2) 
  # %>%
  #   dplyr::mutate(model_name=recode_factor(model_name,`BLP`="LP"))
  frame <- rbind(frame1,frame2)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
                  )
  if(test){
    frame <- frame %>%
      dplyr::filter(season==paste0("20",as.numeric(test_seas)-1,"/20",test_seas)) %>%
      dplyr::group_by(target,location,season,bin_start_incl,Method) %>%
      dplyr::mutate(rand_pit=sample(cdf_vals,1)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c("calendar_week","cdf_vals","value","ls")) %>%
      distinct() 
  } else {
    frame <- frame %>%
      dplyr::filter(season!=paste0("20",as.numeric(test_seas)-1,"/20",test_seas)) %>%
      dplyr::group_by(target,location,season,bin_start_incl,Method) %>%
      dplyr::mutate(rand_pit=sample(cdf_vals,1)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c("calendar_week","cdf_vals","value","ls")) %>%
      distinct() 
  }
  ## create relative observed frequency
  ggplot(frame, aes(rand_pit,color=Method)) +
    stat_ecdf(geom = "step")+
    geom_abline(intercept = 0, slope = 1, linetype =3)+
    facet_wrap(~target,ncol=2)+
    labs(title = paste0("20",as.numeric(test_seas)-1,"/20",test_seas),
         x = "PIT values",
         y = "CDF")+ 
    scale_colour_manual(
      # values = c(RColorBrewer::brewer.pal(11,"RdBu")[c(4,3,8,10)],
      #            RColorBrewer::brewer.pal(11,"PuOr")[c(8,10)]),
      values = c("#fdc086","#d95f02","#beaed4","#7570b3","#7fc97f","#1b9e77"),
      labels = parse_format())+
    guides(fill=guide_legend(ncol=6)) +
    # scale_colour_discrete(values = RColorBrewer::brewer.pal(6,"Dark2"),
    #                       labels = parse_format())+
    theme_bw()+
    theme(legend.text.align = 0,
          legend.position = "bottom",
          plot.title = element_text(size=9),
          axis.text.x=element_text(size=rel(0.8)),
          axis.text.y=element_text(size=rel(0.8)),
          axis.title=element_text(size=8),
          plot.margin=unit(c(0.05,0.15,0.05,0.2), "cm"),
          strip.text.x = element_text(size=8,margin = margin(.1, 0, .1, 0, "cm")))+
    scale_x_continuous(expand = c(0, 0),labels=seq(0,1,0.2), breaks=seq(0,1,0.2), limits=c(0,1)) + 
    scale_y_continuous(expand = c(0, 0),labels=seq(0,1,0.2), breaks=seq(0,1,0.2), limits=c(0,1)) 
}

reli_plot2 <- function(file1,file2,file3,file4,file5,file6,test=TRUE){
  tmp <- lapply(file1, FUN=read.csv)
  frame1 <- do.call(rbind.data.frame, tmp)
  tmp2 <- lapply(file2, FUN=read.csv)
  frame2 <- do.call(rbind.data.frame, tmp2) 
  tmp3 <- lapply(file3, FUN=read.csv)
  frame3 <- do.call(rbind.data.frame, tmp3)
  tmp4 <- lapply(file4, FUN=read.csv)
  frame4 <- do.call(rbind.data.frame, tmp4) 
  tmp5 <- lapply(file5, FUN=read.csv)
  frame5 <- do.call(rbind.data.frame, tmp5)
  tmp6 <- lapply(file5, FUN=read.csv)
  frame6 <- do.call(rbind.data.frame, tmp6) 
  # %>%
  #   dplyr::mutate(model_name=recode_factor(model_name,`BLP`="LP"))
  frame01 <- rbind(frame1,frame2)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  frame02 <- rbind(frame3,frame4)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  frame03 <- rbind(frame5,frame6)  %>%
    dplyr::mutate(Method=recode_factor(model_name, 
                                       `EW` = "EW-LP",
                                       `TLP` = "LP",
                                       `EW_BLP` = "EW-BLP",
                                       `BLP` = "BLP",
                                       `EW_BMC2` = "EW-BMC[2]",
                                       `BMC2` = "BMC[2]")
    )
  # have to fix these
  if(test){
    frame01 <- frame01 %>%
      dplyr::filter(season==paste0("20",as.numeric(17)-1,"/20",17))
    frame02 <- frame02 %>%
      dplyr::filter(season==paste0("20",as.numeric(18)-1,"/20",18)) 
    frame03 <- frame03 %>%
      dplyr::filter(season==paste0("20",as.numeric(19)-1,"/20",19)) 
  } else {
    frame01 <- frame01 %>%
      dplyr::filter(season!=paste0("20",as.numeric(17)-1,"/20",17))
    frame02 <- frame02 %>%
      dplyr::filter(season!=paste0("20",as.numeric(18)-1,"/20",18))
    frame03 <- frame03 %>%
      dplyr::filter(season!=paste0("20",as.numeric(19)-1,"/20",19))
  }
  frame <- rbind(frame01,frame02,frame03) %>%
    dplyr::group_by(target,location,season,bin_start_incl,Method) %>%
    dplyr::mutate(rand_pit=sample(cdf_vals,1)) %>%
    dplyr::ungroup() 
  ## create relative observed frequency
  ggplot(frame, aes(rand_pit,color=Method)) +
    stat_ecdf(geom = "step")+
    geom_abline(intercept = 0, slope = 1, linetype =3)+
    facet_wrap(~target,ncol=2)+
    labs(x = "PIT values",
         y = "CDF")+ 
    ggtitle(ifelse(test,"(b) Probability plots of ensemble forecasts in the test period",
                   "(a) Probability plots of ensemble forecasts in the training period"))+
    scale_colour_manual(
      # values = c(RColorBrewer::brewer.pal(11,"RdBu")[c(4,3,8,10)],
      #            RColorBrewer::brewer.pal(11,"PuOr")[c(8,10)]),
      values = c("#fdc086","#d95f02","#beaed4","#7570b3","#7fc97f","#1b9e77"),
      labels = parse_format())+
    guides(fill=guide_legend(ncol=6)) +
    # scale_colour_discrete(values = RColorBrewer::brewer.pal(6,"Dark2"),
    #                       labels = parse_format())+
    theme_bw()+
    theme(legend.text.align = 0,
          legend.position = "bottom",
          plot.title = element_text(size=12),
          axis.text=element_text(size=rel(0.8)),
          plot.margin=unit(c(0.1,0.25,0.2,0.2), "cm"))+
    scale_x_continuous(expand = c(0, 0),labels=seq(0,1,0.2), breaks=seq(0,1,0.2), limits=c(0,1)) + 
    scale_y_continuous(expand = c(0, 0),labels=seq(0,1,0.2), breaks=seq(0,1,0.2), limits=c(0,1))
  
}

mean_heat_all <- function(yscore){
  mod_order <- c("EW","TLP","EW_BLP","BLP","EW_BMC2","BMC2")
  name <- c(expression("EW-LP"),expression("LP"),expression("EW-BLP"),expression("BLP"),
            expression("EW-BMC"[2]),expression("BMC"[2]))
  yscore %>%
    ggplot(aes(factor(`Model Name`,levels=mod_order),test_train)) +
    geom_tile(aes(fill=factor(rank))) +
    theme(axis.text.x=element_text(size=rel(0.8)),
          axis.text.y=element_text(size=rel(0.8)))+
    xlab("") +
    ylab("") +
    labs(title="Average across all targets and seasons")+
    scale_fill_manual(values = RColorBrewer::brewer.pal(9,"PuBu")[c(2:7)],
                      name="rank")+
    geom_text(aes(label = sprintf("%0.2f", round(mean_across, 2))),size=3) +
    guides(fill=guide_legend(ncol=6)) +
    theme(plot.title = element_text(size=10),
          legend.title = element_text(size=7),
          legend.key.size = unit(0.3, 'cm'),
          legend.position = "bottom",
          legend.text = element_text(size=6),
          plot.margin=unit(c(0,2,1,2),"cm"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))+
    scale_x_discrete(labels=name,expand = c(0, 0))+ 
    scale_y_discrete(expand = c(0, 0))
}

mean_heat_all2 <- function(yscore,titlen){
  mod_order <- c("EW","TLP","EW_BLP","BLP","EW_BMC2","BMC2")
  name <- c(expression("EW-LP"),expression("LP"),expression("EW-BLP"),expression("BLP"),
            expression("EW-BMC"[2]),expression("BMC"[2]))
  yscore %>%
    ggplot(aes(factor(`Model Name`,levels=mod_order),test_train)) +
    geom_tile(aes(fill=factor(rank))) +
    facet_wrap(~`Test Season`,ncol=2)+
    xlab("") +
    ylab("") +
   # labs(title=titlen)+
    ggtitle(titlen) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9,"PuBu")[c(2:7)],
                      name="rank")+
    geom_text(aes(label = sprintf("%0.2f", round(mean_across, 2))),size=2.7) +
    guides(fill=guide_legend(ncol=6)) +
    theme(plot.title = element_text(size=10),
          axis.text.x=element_text(size=rel(0.73)),
          axis.text.y=element_text(size=rel(0.8)),
          legend.title = element_text(size=7),
          legend.key.size = unit(0.3, 'cm'),
          legend.position = "bottom",
          legend.text = element_text(size=7),
          plot.margin=unit(c(0.2,0.2,-0.4,0.2),"cm"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text.x = element_text(margin = margin(.13, 0, .13, 0, "cm"),size=7))+
    scale_x_discrete(labels=name,expand = c(0, 0))+ 
    scale_y_discrete(expand = c(0, 0))
}

mean_heat_all3 <- function(yscore,titlen){
  mod_order <- c("EW","TLP","EW_BLP","BLP","EW_BMC2","BMC2")
  name <- c(expression("EW-LP"),expression("LP"),expression("EW-BLP"),expression("BLP"),
            expression("EW-BMC"[2]),expression("BMC"[2]))
  yscore %>%
    ggplot(aes(factor(`Model Name`,levels=mod_order),test_train)) +
    geom_tile(aes(fill=factor(rank))) +
    facet_wrap(~Target,ncol=2)+
    xlab("") +
    ylab("") +
   # labs(title=titlen)+
    ggtitle(titlen) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9,"PuBu")[c(2:7)],
                      name="rank")+
    geom_text(aes(label = sprintf("%0.2f", round(mean_across, 2))),size=2.7) +
    guides(fill=guide_legend(ncol=6)) +
    theme(plot.title = element_text(size=10),
          axis.text.x=element_text(size=rel(0.73)),
          axis.text.y=element_text(size=rel(0.8)),
          legend.title = element_text(size=7),
          legend.key.size = unit(0.3, 'cm'),
          legend.position = "bottom",
          legend.text = element_text(size=7),
          plot.margin=unit(c(0.2,0.2,-0.4,0.2),"cm"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text.x = element_text(margin = margin(.13, 0, .13, 0, "cm"),size=7))+
    scale_x_discrete(labels=name,expand = c(0, 0))+ 
    scale_y_discrete(expand = c(0, 0))
}


