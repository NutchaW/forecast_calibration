library(tidyverse)
library(scales)
reli_plot <- function(test_seas){
  some_files <- list.files(paste0("calibration_work/ensembles/pit_ls_frame/train_test/",test_seas,"/"),
                           full.names=TRUE, recursive = TRUE)
  tmp <- lapply(some_files, FUN=read.csv)
  frame1 <- do.call(rbind.data.frame, tmp)
  some_files1 <- list.files(paste0("calibration_work/ensembles/pit_ls_frame/train_test/EW_TLP/",test_seas,"/"),
                            full.names=TRUE, recursive = TRUE)
  tmp2 <- lapply(some_files1, FUN=read.csv)
  frame2 <- do.call(rbind.data.frame, tmp2) %>%
    dplyr::mutate(model_name=ifelse(model_name=="EW","EW","TLP")) 
  frame <- rbind(frame1,frame2)  %>%
    dplyr::filter(season==paste0("20",as.numeric(test_seas)-1,"/20",test_seas)) 
  ## create relative observed frequency
  ggplot(frame, aes(cdf_vals,color=model_name)) +
    stat_ecdf(geom = "step")+
    facet_wrap(~target,ncol=2)
}



reli_plot(17)

# check <- EW %>%
#   dplyr::left_join(truth,
#                    by=c("calendar_week"="Calendar.Week","season"="Season","target"="Target","location"="Location")) %>%
#   dplyr::select(-c("Model.Week","Year")) 




ggplot(check,aes(bin_prob_val,m_rel_obs)) +
  geom_line(aes(color=model_name)) +
  facet_wrap(~target,ncol=2)


