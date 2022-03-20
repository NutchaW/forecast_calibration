library(tidyverse)
library(purrr)
source("ensemble_code/BLPwork_functions_app/functions_app.R")
source("ensemble_code/BLPwork_functions_app/ls_pit_reli_funs.R")
#--------------------------------------------------cv build-----------------------------------------------------#
# produce a data frame with train and valid scores for mixture methods
# and select methods for the whole training set for each test season
model_name_list <- c("BMC2",
                     "BMC3",
                     "BMC4",
                     "BMC5",
                     "EW_BMC2",
                     "EW_BMC3",
                     "EW_BMC4",
                     "EW_BMC5")
targets <- c("1 wk ahead", "2 wk ahead", "3 wk ahead", "4 wk ahead")
test_s <- c(17:19)
#--------------------------------------------------EW_BLP----------------------------------------------------#
# make pdf for training
for (ts in test_s) {
  load(paste0("ensemble_code/train_params/EW_BLP", ts, ".rda"))
  if (ts == 17) {
    seas <- 11:16
  } else if (ts == 18) {
    seas <- 11:17
  } else if (ts == 19) {
    seas <- 11:18
  }
  pars <- get(paste0("EW_BLP", ts))[-1]
  for (i in 1:4) {
    tar <- targets[i]
    for (j in 1:4) {
      m_name <- model_name_list[5:length(model_name_list)][j]
      k <- as.numeric(as.character(strsplit(m_name, "")[[1]][7]))
      for (l in 1:length(seas)) {
        if (ts == 17) {
          pars_use <- pars[[j + (4 * (i - 1))]][[l]]
        } else {
          pars_use <- pars[[i + (4 * (j - 1))]][[l]]
        }
        frame_cv <- map_dfr(seas,
                            function(x)
                              make_pdfs(m_name, pars_use, k, x, tar))
        # save in a structure compatible with the function to calculate log score and pit
        write.csv(
          frame_cv,
          file = paste0(
            "ensemble_code/ensembles/cv_ensembles/",
            ts,
            "/",
            seas[l],
            "/",
            m_name,
            "/target",
            i,
            ".csv"
          ),
          row.names = FALSE
        )
      }
    }
  }
}


# calculate log score,pit
for (ts in test_s) {
  if (ts == 17) {
    seas <- 11:16
  } else if (ts == 18) {
    seas <- 11:17
  } else if (ts == 19) {
    seas <- 11:18
  }
  for (j in 1:4) {
    m_name <- model_name_list[5:length(model_name_list)][j]
    for (l in 1:length(seas)) {
      ls_frame <- ls_pit_calc(m_name, ts, seas[l], cv = TRUE, tar = NULL)
      write.csv(
        ls_frame,
        file = paste0(
          "ensembles/pit_ls_frame/cv/EW_BMC/",
          ts,
          "/",
          m_name,
          "_valid",
          l,
          ".csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# calculate cross-validated log score and pit
for (ts in test_s) {
  cv_ewbmc <- cv_frame("EW_BMC", ts)
  save(cv_ewbmc,
       file = paste0("ensembles/pit_ls_frame/cv/EW_BMC_", ts, ".rda"))
}
#------------------------------------------------------BLP--------------------------------------------------------------#
# make pdf for training
for (ts in test_s) {
  load(paste0("ensemble_code/train_params/BLP", ts, ".rda"))
  if (ts == 17) {
    seas <- 11:16
  } else if (ts == 18) {
    seas <- 11:17
  } else if (ts == 19) {
    seas <- 11:18
  }
  pars <- get(paste0("BLP", ts))[-1]
  for (i in 1:4) {
    tar <- targets[i]
    for (j in 1:4) {
      m_name <- model_name_list[j]
      k <- as.numeric(as.character(strsplit(m_name, "")[[1]][4]))
      for (l in 1:length(seas)) {
        if (ts == 17 && m_name == "BMC2") {
          pars_use <- pars[[1]][[i]][[l]]
        } else if (ts == 19 && m_name == "BMC5") {
          pars_use <- pars[[(j - 1) + i]][[l]]
        } else {
          pars_use <- pars[[j]][[l]][[i]]
        }
        frame_cv <- map_dfr(seas,
                            function(x)
                              make_pdfs(m_name, pars_use, k, x, tar))
        # save in a structure compatible with the function to calculate log score and pit
        write.csv(
          frame_cv,
          file = paste0(
            "ensemble_code/ensembles/cv_ensembles/",
            ts,
            "/",
            seas[l],
            "/",
            m_name,
            "/target",
            i,
            ".csv"
          ),
          row.names = FALSE
        )
      }
    }
  }
}


# calculate log score, pit
for (ts in test_s) {
  if (ts == 17) {
    seas <- 11:16
  } else if (ts == 18) {
    seas <- 11:17
  } else if (ts == 19) {
    seas <- 11:18
  }
  for (j in 1:4) {
    m_name <- model_name_list[j]
    for (l in 1:length(seas)) {
      ls_frame <- ls_pit_calc(m_name, ts, seas[l], cv = TRUE, tar = NULL)
      write.csv(
        ls_frame,
        file = paste0(
          "ensembles/pit_ls_frame/cv/BMC/",
          ts,
          "/",
          m_name,
          "_valid",
          l,
          ".csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# calculate cv
# calulate log score, pit
for (ts in test_s) {
  cv_bmc <- cv_frame("BMC", ts)
  save(cv_bmc,
       file = paste0("ensembles/pit_ls_frame/cv/BMC_", ts, ".rda"))
}

#--------------------------------------------------combine cv scores--------------------------------------------------#
## manual
runpre <- function(list_info, sea) {
  table <- list_info[[1]]
  tab <- table %>%
    dplyr::mutate(test_season = sea)
  return(tab)
}

for (s in 17:19) {
  load(paste0(
    "ensemble_code/ensembles/pit_ls_frame/cv/",
    "BMC",
    "_",
    s,
    ".rda"
  ))
  load(paste0(
    "ensemble_code/ensembles/pit_ls_frame/cv/",
    "EW_BMC",
    "_",
    s,
    ".rda"
  ))
  assign(paste0("bmc_", s), cv_bmc)
  assign(paste0("ewbmc_", s), cv_ewbmc)
}

bigtab <-
  rbind(
    runpre(bmc_17, "2016/2017"),
    runpre(bmc_18, "2017/2018"),
    runpre(bmc_19, "2018/2019"),
    runpre(ewbmc_17, "2016/2017"),
    runpre(ewbmc_18, "2017/2018"),
    runpre(ewbmc_19, "2018/2019")
  )
write.csv(bigtab, file = "ensemble_code/ensembles/pit_ls_frame/cv/cv_model.csv",
          row.names = FALSE)

sea_cv <- function(list_info) {
  table <- list_info[[2]]
  tab <- table %>%
    dplyr::mutate(valids = paste0("20", as.numeric(valid_season) - 1, "/20", valid_season))
  return(tab)
}

for (ts in test_s) {
  check1 <- get(paste0("bmc_", ts))
  check2 <- get(paste0("ewbmc_", ts))
  sea_tab <- rbind(sea_cv(check1), sea_cv(check2))
  write.csv(
    sea_tab,
    file = paste0(
      "ensemble_code/ensembles/pit_ls_frame/cv/tables/sea",
      ts,
      ".csv"
    ),
    row.names = FALSE
  )
}
#---------------------------------------------- training and test forecasts ---------------------------------------------#
# make pdf for training
for (ts in test_s) {
  if (ts == 17) {
    seas <- 11:17
  } else if (ts == 18) {
    seas <- 11:18
  } else if (ts == 19) {
    seas <- 11:19
  }
  for (i in 1:4) {
    tar <- targets[i]
    load(
      paste0(
        "ensemble_code/train_params/final_train/target",
        i,
        "_sea",
        ts,
        "_params.rda"
      )
    )
    for (j in 1:4) {
      pars_use <- unlist(get(paste0("target", i, "_sea", ts, "_params"))[j])
      m_name <-
        names(get(paste0("target", i, "_sea", ts, "_params"))[j])
      if (m_name == "BLP" | m_name == "EW_BLP") {
        k <- 1
      } else if (j == 2) {
        k <- as.numeric(as.character(strsplit(m_name, "")[[1]][4]))
      } else {
        k <- as.numeric(as.character(strsplit(m_name, "")[[1]][7]))
      }
      frame_cv <-
        map_dfr(seas, function(x)
          make_pdfs(m_name, pars_use, k, x, tar))
      # save in a structure compatible with the function to calculate log score and pit
      write.csv(
        frame_cv,
        file = paste0(
          "ensemble_code/ensembles/final_ensembles/",
          ts,
          "/target",
          i,
          "/",
          m_name,
          ".csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# calculate log score, pit
for (ts in test_s) {
  if (ts == 17) {
    seas <- 11:17
  } else if (ts == 18) {
    seas <- 11:18
  } else if (ts == 19) {
    seas <- 11:19
  }
  for (j in 1:4) {
    ls_frame <-
      map_dfr(seas, function(x)
        ls_pit_calc(NULL, seas, NULL, cv = FALSE, tar = j))
    write.csv(
      ls_frame,
      file = paste0(
        "ensemble_code/ensembles/pit_ls_frame/train_test/",
        ts,
        "/target",
        j,
        ".csv"
      ),
      row.names = FALSE
    )
  }
}
# build a mean log score data frame
test_s <- 17:19
frame <- map_dfr(test_s, function(x)
  mean_frame(x)) %>%
  dplyr::filter(!(model_name %in% c("EW_BMC3", "EW_BMC5")))
write.csv(
  frame,
  file = paste0(
    "ensemble_code/ensembles/pit_ls_frame/train_test/ls_tw_frame.csv"
  ),
  row.names = FALSE
)

### ---------------------------------------- build EW-TLP and TLP forecasts ---------------------------------------####
# train TLP
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
pdfdat <-
  read.csv("ensemble_code/data_transformed/pdf_modelcol.csv")
targets <- paste0(1:4, " wk ahead")
t_set <- c("2016/2017", "2017/2018", "2018/2019")
mods <- colnames(pdfdat[, 9:35])

tlp_params2 <- data.frame()
for (i in 1:4) {
  for (j in 2) {
    tlp_pars <- make_ensemble_emp(pdfdat, t_set[j], targets[i])
    small_set <-
      data.frame(cbind(
        mods,
        rep(t_set[j], length(mods)),
        rep(targets[i], length(mods)),
        tlp_pars
      ))
    names(small_set) <-
      c("component_model_id", "season", "target", "weight")
  }
  tlp_params2 <- rbind(tlp_params2, small_set)
}
tlp_pars <- rbind(tlp_params1, tlp_params2, tlp_params)
write.csv(
  tlp_pars,
  file = paste0("ensemble_code/train_params/TLP.csv"),
  row.names = FALSE
)

# make pdf for training
enss <- c("EW", "TLP")
for (ts in test_s) {
  if (ts == 17) {
    seas <- 11:17
  } else if (ts == 18) {
    seas <- 11:18
  } else if (ts == 19) {
    seas <- 11:19
  }
  for (i in 1:4) {
    tar <- targets[i]
    for (j in 1:length(enss)) {
      frame_cv <-
        map_dfr(seas, function(x)
          make_pdfs_fsn(enss[j], x, tar, ts))
      # save in a structure compatible with the function to calculate log score and pit
      write.csv(
        frame_cv,
        file = paste0(
          "ensemble_code/ensembles/final_ensembles/EW_TLP/",
          ts,
          "/target",
          i,
          "/",
          enss[j],
          ".csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# calculate ls,pit
for (ts in test_s) {
  if (ts == 17) {
    seas <- 11:17
  } else if (ts == 18) {
    seas <- 11:18
  } else if (ts == 19) {
    seas <- 11:19
  }
  for (j in 1:4) {
    ls_frame <- map_dfr(test_s, function(x)
      ls_pit_calc_fsn2(x, tar = j))
    write.csv(
      ls_frame,
      file = paste0(
        "ensemble_code/ensembles/pit_ls_frame/train_test/EW_TLP/",
        ts,
        "/target",
        j,
        ".csv"
      ),
      row.names = FALSE
    )
  }
}

# build a mean log score data frame
frame <- map_dfr(test_s, function(x)
  mean_frame_fsn(x))
write.csv(
  frame,
  file = paste0(
    "ensemble_code/ensembles/pit_ls_frame/train_test/ls_fsn_frame.csv"
  ),
  row.names = FALSE
)

