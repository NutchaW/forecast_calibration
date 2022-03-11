library(tidyverse)
library(cdcfluview)
library(Matrix)
library(data.table)
library(reshape2)
library(ranger)
library(xtable)
library(here)
library(stats)
library(mixtools)
library(grid)
library(FluSight)
library(rmutil)
library(R.utils)

# combine forecast and truth files per season
# separate all dis data into different target/season
all_season <-
  c(
    "2010/2011",
    "2011/2012",
    "2012/2013",
    "2013/2014",
    "2014/2015",
    "2015/2016",
    "2016/2017",
    "2017/2018",
    "2018/2019"
  )
targets <- c("1 wk ahead", "2 wk ahead", "3 wk ahead", "4 wk ahead")
load("./data_transformed/model_list.rda")
allDis <- read.csv("./data_transformed/rawdata/allDis.csv")
model_list <- unique(as.character(allDis$model_name))[c(1:27, 29, 31)]
componentModel_list <- model_list[1:27]
# loop to write all files by model/target/season
for (i in 1:length(all_season)) {
  for (j in 1:length(targets)) {
    for (k in 1:length(model_list)) {
      componentMod <- allDis %>%
        dplyr::filter(
          target == targets[j],
          season == all_season[i],
          model_name == model_list[k],
          calendar_week %in% c(43:53, 1:18)
        )
      # make sub folders for season-target
      write.csv(
        componentMod,
        file = paste0(
          "./data_transformed/rawdata/",
          substr(all_season[i], 8, 9),
          "/t",
          j,
          "/",
          model_list[k],
          ".csv"
        ),
        quote = FALSE,
        row.names = FALSE
      )
    }
  }
}

## cdf
for (i in 11:19) {
  for (j in 1:4) {
    for (k in 1:length(model_list)) {
      model_data <-
        read.csv(paste0(
          "./data_transformed/rawdata/",
          i,
          "/t",
          j,
          "/",
          model_list[k],
          ".csv"
        )) %>%
        dplyr::mutate(
          calendar_week = factor(calendar_week, levels = c(43:53, 1:18)),
          location = factor(location, levels = unique(location)),
          bin_start_incl = as.numeric(as.character(bin_start_incl))
        ) %>%
        dplyr::group_by(location, calendar_week) %>%
        dplyr::arrange(bin_start_incl) %>%
        dplyr::mutate(cdf = cumsum(value)) %>%
        ungroup() %>%
        data.frame()
      # cdf<-get_cdf(model_data)
      # cdfdata<-data.frame(cbind(model_data,cdf))
      write.csv(
        model_data,
        file = paste0(
          "./data_transformed/cdf/",
          i,
          "/t",
          j,
          "/",
          model_list[k],
          ".csv"
        ),
        quote = FALSE,
        row.names = FALSE
      )
    }
  }
}

# calculate PIT by target/season/model/week
truths <- read.csv("scores/target-multivals.csv") %>%
  dplyr::filter(Target %in% targets, Calendar.Week %in% c(43:53, 1:18))
truths$Season <- as.character(truths$Season)
truths$Location <- as.character(truths$Location)
truths$Target <- as.character(truths$Target)
names(truths) <- tolower(colnames(truths))

temp <- truths %>%
  dplyr::mutate(
    calendar_week = calendar.week,
    valid.bin_start_incl = round(as.numeric(valid.bin_start_incl), 1),
    valid.bin_lag = ifelse((as.numeric(
      valid.bin_start_incl
    ) - 0.1) < 0,
    round(0, 1),
    round((
      as.numeric(valid.bin_start_incl) - 0.1
    ), 1))
  ) %>%
  dplyr::select(-c("year", "model.week", "calendar.week"))
temp <- temp[, c(2:3, 1, 5, 4, 6)]

for (i in 11:19) {
  for (j in 1:4) {
    for (k in 1:length(model_list)) {
      temp1 <- temp %>%
        dplyr::filter(target == targets[j],
                      season == all_season[i - 10]) %>%
        dplyr::arrange(factor(location, levels = unique(location)),
                       factor(calendar_week, levels = c(43:53, 1:18)))
      model_data <-
        read.csv(paste0(
          "./data_transformed/cdf/",
          i,
          "/t",
          j,
          "/",
          model_list[k],
          ".csv"
        )) %>%
        dplyr::select(
          "location",
          "target",
          "season",
          "calendar_week",
          "bin_start_incl",
          "cdf",
          "model_name"
        )
      pitdata <- temp1 %>%
        dplyr::select(-"valid.bin_lag") %>%
        dplyr::left_join(
          model_data,
          by = c(
            "valid.bin_start_incl" = "bin_start_incl",
            "location" = "location",
            "target" = "target",
            "season" = "season",
            "calendar_week" = "calendar_week"
          )
        )
      pitdata_lag <- temp1 %>%
        dplyr::select(-"valid.bin_start_incl") %>%
        dplyr::left_join(
          model_data,
          by = c(
            "valid.bin_lag" = "bin_start_incl",
            "location" = "location",
            "target" = "target",
            "season" = "season",
            "calendar_week" = "calendar_week"
          )
        )
      write.csv(
        pitdata[, c(7, 1:6)],
        file = paste0(
          "./data_transformed/pit/",
          i,
          "/t",
          j,
          "/",
          model_list[k],
          ".csv"
        ),
        quote = FALSE,
        row.names = FALSE
      )
      write.csv(
        pitdata_lag[, c(7, 1:6)],
        file = paste0(
          "./data_transformed/pit_lag/",
          i,
          "/t",
          j,
          "/",
          model_list[k],
          ".csv"
        ),
        quote = FALSE,
        row.names = FALSE
      )
    }
  }
}

# get combined pits

fix_files <- function(filename) {
  require(dplyr)
  dat <- read.csv(filename)
  dis <- dat %>%
    dplyr::select(-"valid.bin_start_incl")
  return(dis)
}
fix_files2 <- function(filename) {
  require(dplyr)
  dat <- read.csv(filename)
  dis <- dat %>%
    dplyr::select(-"valid.bin_lag")
  return(dis)
}

## get all model files without the metadata
some_files <-
  list.files("./data_transformed/pit/",
             full.names = TRUE,
             recursive = TRUE)

## extract point estimates and put into one dataframe
tmp <- lapply(some_files, FUN = fix_files)
tmp.df <- do.call(rbind.data.frame, tmp)
write.csv(tmp.df,
          file = "./data_transformed/pit/combined_pit.csv",
          quote = FALSE,
          row.names = FALSE)

sc <- read.csv("./data_transformed/scores_rename.csv") %>%
  dplyr::filter(Epiweek %in% c(43:53, 1:18), Model %in% model_list)
pt <- read.csv("./data_transformed/pit/combined_pit.csv")

sc$Model <- as.character(sc$Model)
pt$model_name <- as.character(pt$model_name)

combined_pitscore <- pt %>%
  dplyr::left_join(
    sc,
    by = c(
      model_name = "Model",
      location = "Location",
      target = "Target",
      season = "Season",
      calendar_week = "Epiweek"
    )
  ) %>%
  dplyr::select(-"Multi.bin.score")

# transform
library(data.table)
library(reshape2)

comb2 <- combined_pitscore %>%
  dplyr::select(-"Model.Week", -"Year", -"Model.Week", -"Score") %>%
  dcast(target + location + season + calendar_week ~ model_name, value.var =
          "cdf") %>%
  dplyr::select(-c("target-based-weights", "equal-weights"))

write.csv(comb2,
          file = "./data_transformed/pit_modelcol.csv",
          quote = FALSE,
          row.names = FALSE)

# repeat the same but for prepit
some_files <-
  list.files("./data_transformed/pit_lag/",
             full.names = TRUE,
             recursive = TRUE)

## extract point estimates and put into one dataframe
tmp <- lapply(some_files, FUN = fix_files2)
tmp.df <- do.call(rbind.data.frame, tmp)
write.csv(tmp.df,
          file = "./data_transformed/pit_lag/combined_prepit.csv",
          quote = FALSE,
          row.names = FALSE)

sc <- read.csv("./data_transformed/scores_rename.csv") %>%
  dplyr::filter(Epiweek %in% c(43:53, 1:18), Model %in% model_list)
prept <- read.csv("./data_transformed/pit_lag/combined_prepit.csv")

sc$Model <- as.character(sc$Model)
prept$model_name <- as.character(prept$model_name)

combined_prepitscore <- prept %>%
  dplyr::left_join(
    sc,
    by = c(
      model_name = "Model",
      location = "Location",
      target = "Target",
      season = "Season",
      calendar_week = "Epiweek"
    )
  ) %>%
  dplyr::select(-"Multi.bin.score")

# transpose
comb3 <- combined_prepitscore %>%
  dplyr::select(-"Model.Week", -"Year", -"Score") %>%
  dplyr::filter(model_name != "equal-weights",
                model_name != "target-based-weights") %>%
  dcast(target + location + season + calendar_week ~ model_name, value.var =
          "cdf")

write.csv(comb3,
          file = "./data_transformed/prepit_modelcol.csv",
          quote = FALSE,
          row.names = FALSE)


# get cdf combined for each season

cdf_trans <- function(filename) {
  dat <- read.csv(filename)
  dis <- dat
  return(dis)
}

## get all model files without the metadata
for (i in 11:19) {
  some_files <-
    list.files(
      paste0("./data_transformed/cdf/", i, "/"),
      full.names = TRUE,
      recursive = TRUE
    )
  tmp <- lapply(some_files, FUN = cdf_trans)
  tmp.df <- do.call(rbind.data.frame, tmp)
  write.csv(
    tmp.df,
    file = paste0("./data_transformed/cdf/", i, "/cdf", i, ".csv"),
    quote = FALSE,
    row.names = FALSE
  )
}


for (i in 11:19) {
  for (j in 1:4) {
    data1 <-
      read.csv(file = paste0("./data_transformed/cdf/", i, "/cdf", i, ".csv")) %>%
      dplyr::select(-"forecast_week", -"year", -"value") %>%
      reshape2::dcast(
        location + target + bin_start_incl + bin_end_notincl + calendar_week + season ~ model_name,
        value.var = "cdf"
      ) %>%
      dplyr::select(-"target-based-weights", -"equal-weights")
    write.csv(
      data1,
      file = paste0("./data_transformed/cdf/", i, "/reduced_cdf", i, ".csv"),
      quote = FALSE,
      row.names = FALSE
    )
    data2 <-
      read.csv(file = paste0("./data_transformed/cdf/", i, "/cdf", i, ".csv")) %>%
      dplyr::select(-"forecast_week", -"year", -"cdf") %>%
      reshape2::dcast(
        location + target + bin_start_incl + bin_end_notincl + calendar_week + season ~ model_name,
        value.var = "value"
      ) %>%
      dplyr::select(-"target-based-weights", -"equal-weights")
    write.csv(
      data2,
      file = paste0("./data_transformed/cdf/", i, "/reduced_pdf", i, ".csv"),
      quote = FALSE,
      row.names = FALSE
    )
  }
}


# make empirical pdf table for each year
truths$valid.bin_start_incl <-
  as.numeric(as.character(truths$valid.bin_start_incl))

# pdf paths
pdf_paths <-
  paste0("./data_transformed/cdf/",
         11:19,
         "/reduced_pdf",
         11:19,
         ".csv")
data1 <- map_dfr(pdf_paths, function(path) {
  read.csv(file = path)
})

data2 <- truths %>%
  left_join(
    data1,
    by = c(
      season = "season",
      target = "target",
      location = "location",
      calendar.week = "calendar_week",
      valid.bin_start_incl = "bin_start_incl"
    )
  )
write.csv(
  data2,
  file = paste0("./data_transformed/pdf_modelcol.csv"),
  quote = FALSE,
  row.names = FALSE
)


# make empirical cdf table for each year
cdf_paths <-
  paste0("./data_transformed/cdf/",
         11:19,
         "/reduced_cdf",
         11:19,
         ".csv")
data3 <- map_dfr(cdf_paths, function(path) {
  read.csv(file = path)
})

data4 <- truths %>%
  left_join(
    data3,
    by = c(
      season = "season",
      target = "target",
      location = "location",
      calendar.week = "calendar_week",
      valid.bin_start_incl = "bin_start_incl"
    )
  )
write.csv(
  data4,
  file = paste0("./data_transformed/cdf_modelcol.csv"),
  quote = FALSE,
  row.names = FALSE
)

