library(tidyverse)
library(cdcfluview)
library(Matrix)
library(data.table)
library(reshape2)
library(gridExtra)
library(rmutil)
library(ranger)
library(xtable)
library(here)
library(stats)
library(FluSight)

extract_dis <- function(filename) {
  require(dplyr)
  require(FluSight)
  file_parts <- strsplit(filename, split = "/")[[1]]
  file_only <- file_parts[length(file_parts)]
  message(paste("file:", filename, Sys.time()))
  dat <- read_entry(filename)
  dis <- dat %>%
    dplyr::filter(!(type == "Point")) %>%
    dplyr::select(-type,-unit)
  dis$year <- as.numeric(substr(file_only, 6, 9))
  dis$calendar_week <- as.numeric(substr(file_only, 3, 4))
  dis$season <- ifelse(
    as.numeric(substr(file_only, 3, 4)) > 18,
    paste0(as.numeric(substr(file_only, 6, 9)), "/", as.numeric(substr(file_only, 6, 9)) +
             1),
    paste0(as.numeric(substr(file_only, 6, 9)) - 1, "/", as.numeric(substr(file_only, 6, 9)))
  )
  dis$model_name <-
    substr(file_only, 11, gregexpr("\\.", file_only)[[1]] - 1)
  return(dis)
}
extract_point_ests <- function(filename) {
  require(dplyr)
  require(FluSight)
  file_parts <- strsplit(filename, split = "/")[[1]]
  file_only <- file_parts[length(file_parts)]
  message(paste("file:", filename, Sys.time()))
  dat <- read_entry(filename)
  point_ests <- dat %>%
    dplyr::filter(!(type == "Point")) %>%
    generate_point_forecasts(method = "Expected Value") %>%
    ungroup() %>%
    dplyr::select(location, target, value)
  point_ests$year <- as.numeric(substr(file_only, 6, 9))
  point_ests$calendar_week <- as.numeric(substr(file_only, 3, 4))
  point_ests$model_name <-
    substr(file_only, 11, gregexpr("\\.", file_only)[[1]] - 1)
  return(point_ests)
}

ew_to_seasonweek <- function(EW, year, season_start_week = 30) {
  require(MMWRweek)
  num_days <-
    ifelse(MMWRweek::MMWRweek(as.Date(paste0(
      as.character(year), "-12-28"
    )))[2] == 53, 53, 52)
  return(ifelse(
    EW > season_start_week,
    EW - season_start_week,
    (EW + num_days) - season_start_week
  ))
}

generate_mean <- function(entry) {
  dat <- generate_point_forecasts(entry, method = "Expected Value")
  return(dat)
}

## get all model files without the metadata
some_files <-
  c(
    list.files(
      "model-forecasts/component-models/",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      "model-forecasts/cv-ensemble-models/",
      full.names = TRUE,
      recursive = TRUE
    )
  )
some_files <- some_files[-grep("metadata.txt", some_files)]
some_files <- some_files[-grep("models-for-2017-2018", some_files)]
some_files <- some_files[-grep("model-id-map", some_files)]
some_files <- some_files[-grep("complete-modelids", some_files)]

tmp <- lapply(some_files, FUN = extract_dis)
tmp.df <- do.call(rbind.data.frame, tmp)
write.csv(tmp.df,
          file = "data_code/allDis.csv",
          quote = FALSE,
          row.names = FALSE)


#-------for real time comparison---------#
## get all model files without the metadata
some_files2 <- c(
  list.files(
    "model-forecasts/real-time-component-models/",
    full.names = TRUE,
    recursive = TRUE
  ),
  list.files(
    "model-forecasts/real-time-ensemble-models/",
    full.names = TRUE,
    recursive = TRUE
  )
)
some_files2 <- some_files2[-grep("metadata.txt", some_files2)]
some_files2 <- some_files2[-grep("plots", some_files2)]

## extract point estimates and put into one data frame
tmp2 <- lapply(some_files2, FUN = extract_dis)

tmp2.df <- do.call(rbind.data.frame, tmp2)
write.csv(tmp2.df,
          file = "data_code/rtDis.csv",
          quote = FALSE,
          row.names = FALSE)