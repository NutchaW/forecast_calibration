#---------------------------------------------- input functions ------------------------------------------------#
ecdf_round <- function(q, values) {
  d <- round(q, 1)
  bin <- (d * 10) + 1
  return(cumsum(values)[bin])
}


cdf_to_pdf <- function(cdf_vector) {
  prob_vector <-
    c(cdf_vector[1], c(cdf_vector - dplyr::lag(cdf_vector))[-1])
  pdf <- prob_vector / sum(prob_vector)
  return(pdf)
}


##------------------------------------------------- Functions -----------------------------------------##
# build ensemble pdfs
make_pdfs <- function(ensemble_name, params, K, seas, tar) {
  # get component cdf and number of models
  cdf_frame <-
    read.csv(paste0("data_transformed/cdf/", seas, "/reduced_cdf", seas, ".csv")) %>%
    dplyr::filter(target == tar) %>%
    dplyr::group_by(location, calendar_week) %>%
    dplyr::arrange(location,
                   factor(calendar_week, levels = c(43:53, 1:18)),
                   bin_start_incl) %>%
    dplyr::ungroup()
  # turn cdf into matrix
  comp_cdf_matrix <- cdf_frame[, 7:ncol(cdf_frame)] %>%
    as.matrix()
  M <- ncol(comp_cdf_matrix)
  # frame for pdf
  ensemble_pdf <- cdf_frame[, 1:6]
  # build train and test cdf/pdf part
  if (grepl("EW_", ensemble_name, fixed = TRUE)) {
    # ew-bmc or ew-blp
    ## extract params
    omega <- rep(1 / M, M)
    munu <- params[1:(K * 2)]
    w <- params[(1 + (2 * K)):length(params)]
    # equally-weighed component models (ew ensemble)
    H <- comp_cdf_matrix %*% omega
    # build matrices for each k
    pbeta_mat <- matrix(NA, ncol = K, nrow = nrow(comp_cdf_matrix))
    for (k in 1:K) {
      assign(paste0("ab", k), c(munu[k] * munu[k + K], (1 - munu[k]) * munu[k +
                                                                              K]))
      pbeta_mat[, k] <-
        pbeta(H,
              shape1 = get(paste0("ab", k))[1],
              shape2 = get(paste0("ab", k))[2])
    }
    # turn pbeta_mat into a combined mixture (pbeta_mat %*% w)
    ensemble_pdf$cdf_vals <- pbeta_mat %*% w
    # then turn  a combined mixture into pdf using the cdf to pdf function
    ensemble_cpdf <- ensemble_pdf %>%
      dplyr::group_by(location, calendar_week) %>%
      dplyr::mutate(value = cdf_to_pdf(cdf_vals),
                    model_name = ensemble_name) %>%
      ungroup()
  } else {
    # BMC or blp
    ## extract params
    munu <- params[1:(K * 2)]
    omega <- params[1 + (K * 2):(length(params) - (K + 1))]
    w <- params[(1 + (2 * K) + (M * K)):length(params)]
    # build matrices for each k
    pbeta_mat <- matrix(NA, ncol = K, nrow = nrow(comp_cdf_matrix))
    for (k in 1:K) {
      assign(paste0("ab", k), c(munu[k] * munu[k + K], (1 - munu[k]) * munu[k +
                                                                              K]))
      assign(paste0("omega", k), c())
      for (m in 1:M) {
        assign(paste0("omega", k), c(get(paste0("omega", k)), omega[((m - 1) * K) +
                                                                      k]))
      }
      H <- comp_cdf_matrix %*% c(get(paste0("omega", k)))
      pbeta_mat[, k] <-
        pbeta(H,
              shape1 = get(paste0("ab", k))[1],
              shape2 = get(paste0("ab", k))[2])
    }
    # turn pbeta_mat into a combined mixture (pbeta_mat %*% w)
    ensemble_pdf$cdf_vals <- pbeta_mat %*% w
    # then turn  a combined mixture into pdf using the cdf to pdf function
    ensemble_cpdf <- ensemble_pdf %>%
      dplyr::group_by(location, calendar_week) %>%
      dplyr::mutate(value = cdf_to_pdf(cdf_vals),
                    model_name = ensemble_name) %>%
      ungroup()
  }
  return(ensemble_cpdf)
}

# function
# a function to calculate ls and pits
ls_pit_calc <-
  function(ensemble_name = NULL,
           test_seas,
           valid_seas = NULL,
           cv = TRUE,
           tar = NULL) {
    if (cv) {
      some_files <- list.files(
        paste0(
          "ensembles/cv_ensembles/",
          test_seas,
          "/",
          valid_seas,
          "/",
          ensemble_name,
          "/"
        ),
        full.names = TRUE,
        recursive = TRUE
      )
    } else{
      some_files <-
        list.files(
          paste0(
            "ensembles/final_ensembles/",
            test_seas,
            "/target",
            tar,
            "/"
          ),
          full.names = TRUE,
          recursive = TRUE
        )
    }
    tmp <- lapply(some_files, FUN = read.csv)
    frame <- do.call(rbind.data.frame, tmp)
    tars <- unique(frame$target)
    # read in truth
    truth <- read.csv("scores/target-multivals.csv") %>%
      dplyr::filter(Target %in% tars,
                    Calendar.Week %in% c(43:53, 1:18))
    truth$Season <- as.character(truth$Season)
    truth$Location <- as.character(truth$Location)
    truth$Target <- as.character(truth$Target)
    ##
    res_frame <- frame %>%
      dplyr::left_join(
        truth,
        by = c(
          "calendar_week" = "Calendar.Week",
          "season" = "Season",
          "target" = "Target",
          "location" = "Location"
        )
      ) %>%
      dplyr::select(-c("Model.Week", "Year")) %>%
      dplyr::group_by(season, calendar_week, target, location) %>%
      dplyr::filter(round(as.numeric(Valid.Bin_start_incl), 1) == round(as.numeric(bin_start_incl), 1)) %>%
      dplyr::mutate(ls = ifelse(value == 0 |
                                  log(value) < -10, -10, log(value))) %>%
      ungroup()
    if (cv) {
      res_frame <- res_frame %>%
        dplyr::mutate(valid_season = valid_seas)
      if (!("model_name" %in% colnames(res_frame))) {
        res_frame <- res_frame %>%
          dplyr::mutate(model_name = ensemble_name)
      }
    }
    return(res_frame)
  }

# save results in separate folders for EW and bmc
# build cv frame= continue here, make mean score that take res frame
cv_frame <- function(ew_ind, test_s) {
  some_files <- list.files(
    paste0("ensembles/pit_ls_frame/cv/",
           ew_ind, "/", test_s, "/"),
    full.names = TRUE,
    recursive = TRUE
  )
  tmp <- lapply(some_files, FUN = read.csv)
  frame <- do.call(rbind.data.frame, tmp)
  train_table <- frame %>%
    dplyr::filter(as.numeric(as.character(substr(
      season, nchar(season) - 2 + 1, nchar(season)
    ))) !=
      as.numeric(as.character(valid_season))) %>%
    dplyr::group_by(valid_season, target, model_name) %>%
    dplyr::mutate(mean_ls_train = mean(ls)) %>%
    dplyr::ungroup() %>%
    dplyr::select("target", "model_name", "valid_season", "mean_ls_train") %>%
    distinct()
  valid_table <- frame %>%
    dplyr::filter(as.numeric(as.character(substr(
      season, nchar(season) - 2 + 1, nchar(season)
    ))) ==
      as.numeric(as.character(valid_season))) %>%
    dplyr::group_by(valid_season, target, model_name) %>%
    dplyr::mutate(mean_ls_valid = mean(ls)) %>%
    dplyr::ungroup() %>%
    dplyr::select("target", "model_name", "valid_season", "mean_ls_valid") %>%
    distinct()
  # join table and calculatesd...
  table <- train_table %>%
    dplyr::left_join(valid_table)
  table2 <- table %>%
    dplyr::group_by(target, model_name) %>%
    dplyr::mutate(
      cv_mls_train = mean(mean_ls_train),
      cv_mls_valid = mean(mean_ls_valid)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c("valid_season", "mean_ls_train", "mean_ls_valid")) %>%
    distinct()
  tname <- table2 %>%
    dplyr::mutate(num_k = as.numeric(substr(
      model_name, nchar(model_name), nchar(model_name)
    ))) %>%
    dplyr::group_by(target) %>%
    dplyr::mutate(
      lower_bound = (max(cv_mls_valid) - sd(cv_mls_valid)),
      within_sd = ifelse(cv_mls_valid >= lower_bound, TRUE, FALSE)
    ) %>%
    dplyr::filter(within_sd == TRUE) %>%
    dplyr::filter(num_k == min(num_k)) %>%
    dplyr::ungroup() %>%
    distinct()
  #return(table)
  return(list(
    select = tname,
    cv_info = table,
    cv_frame = table2
  ))
}

mean_frame <- function(test_s) {
  some_files <-
    list.files(
      paste0("ensembles/pit_ls_frame/train_test/", test_s, "/"),
      full.names = TRUE,
      recursive = TRUE
    )
  tmp <- lapply(some_files, FUN = read.csv)
  frame <- do.call(rbind.data.frame, tmp)
  train_table <- frame %>%
    dplyr::mutate(tseason = test_s) %>%
    dplyr::filter(as.numeric(as.character(substr(
      season, nchar(season) - 2 + 1, nchar(season)
    ))) != test_s) %>%
    dplyr::group_by(target, model_name) %>%
    dplyr::mutate(mean_ls_train = mean(ls)) %>%
    dplyr::ungroup() %>%
    dplyr::select("target", "model_name", "mean_ls_train", "tseason") %>%
    distinct()
  valid_table <- frame %>%
    dplyr::mutate(tseason = test_s) %>%
    dplyr::filter(as.numeric(as.character(substr(
      season, nchar(season) - 2 + 1, nchar(season)
    ))) == test_s) %>%
    dplyr::group_by(target, model_name) %>%
    dplyr::mutate(mean_ls_test = mean(ls),
                  test_season = season) %>%
    dplyr::ungroup() %>%
    dplyr::select("target",
                  "model_name",
                  "mean_ls_test",
                  "test_season",
                  "tseason") %>%
    distinct()
  # join table and calculatesd...
  table <- train_table %>%
    dplyr::left_join(valid_table) %>%
    dplyr::select(-"tseason")
  return(table)
}


make_pdfs_fsn <- function(ensemble_name, season, tar, ts) {
  # get component cdf and number of models
  pdf_frame <-
    read.csv(paste0(
      "data_transformed/cdf/",
      season,
      "/reduced_pdf",
      season,
      ".csv"
    )) %>%
    dplyr::filter(target == tar) %>%
    dplyr::group_by(location, calendar_week) %>%
    dplyr::arrange(location,
                   factor(calendar_week, levels = c(43:53, 1:18)),
                   bin_start_incl) %>%
    dplyr::ungroup()
  cdf_frame <-
    read.csv(paste0(
      "data_transformed/cdf/",
      season,
      "/reduced_cdf",
      season,
      ".csv"
    )) %>%
    dplyr::filter(target == tar) %>%
    dplyr::group_by(location, calendar_week) %>%
    dplyr::arrange(location,
                   factor(calendar_week, levels = c(43:53, 1:18)),
                   bin_start_incl) %>%
    dplyr::ungroup()
  # turn cdf into matrix
  comp_cdf_matrix <- cdf_frame[, 7:ncol(cdf_frame)] %>%
    as.matrix()
  # turn pdf into matrix
  comp_pdf_matrix <- pdf_frame[, 7:ncol(pdf_frame)] %>%
    as.matrix()
  M <- ncol(comp_pdf_matrix)
  # frame for pdf
  ensemble_pdf <- pdf_frame[, 1:6]
  #build train and test cdf/pdf part
  if (ensemble_name == "EW") {
    omega <- rep(1 / M, M)
    # then turn  a combined mixture into pdf using the cdf to pdf function
    ensemble_cpdf <- ensemble_pdf %>%
      dplyr::mutate(
        value = comp_pdf_matrix %*% omega,
        cdf_vals = comp_cdf_matrix %*% omega,
        model_name = ensemble_name
      ) %>%
      ungroup()
  } else {
    # BMC or blp
    # change name for bma
    season_c <- paste0("20", ts - 1, "/20", ts)
    omega_frame <-
      read.csv("ensemble_code/train_params/TLP.csv") %>%
      dplyr::filter(target == tar,
                    season == season_c)
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-BMA"] <-
      "CUBMA"
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-EAKFC_SEIRS"] <-
      "CU_EAKFC_SEIRS"
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-EAKFC_SIRS"] <-
      "CU_EAKFC_SIRS"
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-EKF_SEIRS"] <-
      "CU_EKF_SEIRS"
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-EKF_SIRS"] <-
      "CU_EKF_SIRS"
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-RHF_SEIRS"] <-
      "CU_RHF_SEIRS"
    omega_frame$component_model_id[omega_frame$component_model_id == "CU-RHF_SIRS"] <-
      "CU_RHF_SIRS"
    omega_frame$component_model_id[omega_frame$component_model_id == "Delphi-BasisRegression"] <-
      "Delphi_BasisRegression"
    omega_frame$component_model_id[omega_frame$component_model_id == "Delphi-EmpiricalFuture"] <-
      "Delphi_EmpiricalFutures"
    omega_frame$component_model_id[omega_frame$component_model_id == "Delphi-EmpiricalTraj"] <-
      "Delphi_EmpiricalTrajectories"
    omega_frame$component_model_id[omega_frame$component_model_id == "Delphi-DeltaDensity1"] <-
      "Delphi_ExtendedDeltaDensity"
    omega_frame$component_model_id[omega_frame$component_model_id == "Delphi-DeltaDensity2"] <-
      "Delphi_MarkovianDeltaDensity"
    omega_frame$component_model_id[omega_frame$component_model_id == "Delphi-Uniform"] <-
      "Delphi_Uniform"
    omega_frame$component_model_id[omega_frame$component_model_id == "LANL-DBMplus"] <-
      "LANL_DBMplus"
    omega_frame$component_model_id[omega_frame$component_model_id == "Protea-Cheetah"] <-
      "Protea_Cheetah"
    omega_frame$component_model_id[omega_frame$component_model_id == "Protea-Kudu"] <-
      "Protea_Kudu"
    omega_frame$component_model_id[omega_frame$component_model_id == "Protea-Springbok"] <-
      "Protea_Springbok"
    omega_frame$component_model_id[omega_frame$component_model_id == "ReichLab-KCDE_PHB"] <-
      "ReichLab_kcde_backfill_post_hoc"
    omega_frame$component_model_id[omega_frame$component_model_id == "ReichLab-KCDE_NB"] <-
      "ReichLab_kcde_backfill_none"
    omega_frame$component_model_id[omega_frame$component_model_id == "ReichLab-KDE"] <-
      "ReichLab_kde"
    omega_frame$component_model_id[omega_frame$component_model_id == "ReichLab-SARIMA1"] <-
      "ReichLab_sarima_seasonal_difference_FALSE"
    omega_frame$component_model_id[omega_frame$component_model_id == "ReichLab-SARIMA2"] <-
      "ReichLab_sarima_seasonal_difference_TRUE"
    omega_frame$component_model_id[omega_frame$component_model_id == "FSNetwork-TW"] <-
      "target-based-weights"
    omega_frame$component_model_id[omega_frame$component_model_id == "FSNetwork-EW"] <-
      "equal-weights"
    omega_frame$component_model_id[omega_frame$component_model_id == "FluOutlook-Mech"] <-
      "FluOutlook_Mech"
    omega_frame$component_model_id[omega_frame$component_model_id == "FluOutlook-MechAug"] <-
      "FluOutlook_MechAug"
    omega_frame$component_model_id[omega_frame$component_model_id == "FluX-ARLR"] <-
      "FluX_ARLR"
    omega_frame$component_model_id[omega_frame$component_model_id == "FluX-FluX_LSTM"] <-
      "FluX_LSTM"
    omega_frame$component_model_id[omega_frame$component_model_id == "UA-EpiCos"] <-
      "UA_EpiCos"
    ## extract params matching names
    ordered_f <- omega_frame %>%
      dplyr::arrange(factor(component_model_id, levels = colnames(comp_pdf_matrix)))
    omega <- c(ordered_f$weight)
    # then turn  a combined mixture into pdf using the cdf to pdf function
    ensemble_cpdf <- ensemble_pdf %>%
      dplyr::mutate(
        value = comp_pdf_matrix %*% omega,
        cdf_vals = comp_cdf_matrix %*% omega,
        model_name = ensemble_name
      ) %>%
      ungroup()
  }
  return(ensemble_cpdf)
}

ls_pit_calc_fsn <- function(test_seas, tar) {
  some_files <- list.files(
    paste0(
      "ensembles/final_ensembles/EW_TLP/",
      test_seas,
      "/target",
      tar,
      "/"
    ),
    full.names = TRUE,
    recursive = TRUE
  )
  tmp <- lapply(some_files, FUN = read.csv)
  frame <- do.call(rbind.data.frame, tmp)
  tars <- unique(frame$target)
  # read in truth
  truth <- read.csv("scores/target-multivals.csv") %>%
    dplyr::filter(Target %in% tars,
                  Calendar.Week %in% c(43:53, 1:18))
  truth$Season <- as.character(truth$Season)
  truth$Location <- as.character(truth$Location)
  truth$Target <- as.character(truth$Target)
  ##
  res_frame <- frame %>%
    dplyr::left_join(
      truth,
      by = c(
        "calendar_week" = "Calendar.Week",
        "season" = "Season",
        "target" = "Target",
        "location" = "Location"
      )
    ) %>%
    dplyr::select(-c("Model.Week", "Year")) %>%
    dplyr::group_by(season, calendar_week, target, location) %>%
    dplyr::filter(round(as.numeric(Valid.Bin_start_incl), 1) == round(as.numeric(bin_start_incl), 1)) %>%
    dplyr::mutate(ls = ifelse(value == 0 |
                                log(value) < -10, -10, log(value))) %>%
    ungroup()
  return(res_frame)
}

mean_frame_fsn <- function(test_s) {
  some_files <-
    list.files(
      paste0("ensembles/pit_ls_frame/train_test/EW_TLP/", test_s, "/"),
      full.names = TRUE,
      recursive = TRUE
    )
  tmp <- lapply(some_files, FUN = read.csv)
  frame <- do.call(rbind.data.frame, tmp)
  train_table <- frame %>%
    dplyr::mutate(tseason = test_s) %>%
    dplyr::filter(as.numeric(as.character(substr(
      season, nchar(season) - 2 + 1, nchar(season)
    ))) != test_s) %>%
    dplyr::group_by(target, model_name) %>%
    dplyr::mutate(mean_ls_train = mean(ls)) %>%
    dplyr::ungroup() %>%
    dplyr::select("target", "model_name", "mean_ls_train", "tseason") %>%
    distinct()
  valid_table <- frame %>%
    dplyr::mutate(tseason = test_s) %>%
    dplyr::filter(as.numeric(as.character(substr(
      season, nchar(season) - 2 + 1, nchar(season)
    ))) == test_s) %>%
    dplyr::group_by(target, model_name) %>%
    dplyr::mutate(mean_ls_test = mean(ls),
                  test_season = season) %>%
    dplyr::ungroup() %>%
    dplyr::select("target",
                  "model_name",
                  "mean_ls_test",
                  "test_season",
                  "tseason") %>%
    distinct()
  # join table and calculated...
  table <- train_table %>%
    dplyr::left_join(valid_table) %>%
    dplyr::select(-"tseason")
  return(table)
}

ls_pit_calc_fsn2 <- function(test_seas, tar) {
  some_files1 <-
    read.csv(
      paste0(
        "/Users/ahctun_woon/git/end-of-2019-2020/cdc-flusight-ensemble/ensembles/final_ensembles/EW_TLP/",
        test_seas,
        "/target",
        tar,
        "/EW.csv"
      )
    )
  some_files2 <-
    read.csv(
      paste0(
        "/Users/ahctun_woon/git/end-of-2019-2020/cdc-flusight-ensemble/ensembles/final_ensembles/EW_TLP/",
        test_seas,
        "/target",
        tar,
        "/TLP.csv"
      )
    )
  frame <- rbind(some_files1, some_files2)
  tars <- unique(frame$target)
  # read in truth
  truth <- read.csv("scores/target-multivals.csv") %>%
    dplyr::filter(Target %in% tars,
                  Calendar.Week %in% c(43:53, 1:18))
  truth$Season <- as.character(truth$Season)
  truth$Location <- as.character(truth$Location)
  truth$Target <- as.character(truth$Target)
  ##
  res_frame <- frame %>%
    dplyr::left_join(
      truth,
      by = c(
        "calendar_week" = "Calendar.Week",
        "season" = "Season",
        "target" = "Target",
        "location" = "Location"
      )
    ) %>%
    dplyr::select(-c("Model.Week", "Year")) %>%
    dplyr::group_by(season, calendar_week, target, location) %>%
    dplyr::filter(round(as.numeric(Valid.Bin_start_incl), 1) == round(as.numeric(bin_start_incl), 1)) %>%
    dplyr::mutate(ls = ifelse(value == 0 |
                                log(value) < -10, -10, log(value))) %>%
    ungroup()
  return(res_frame)
}


#### function for transformations
# build ensemble pdfs
make_tran_pdfs <-
  function(ensemble_name, params, K, seas, tar, loc, wk) {
    # get component cdf and number of models
    cdf_frame <-
      read.csv(paste0(path, "/BLPdata_app/cdf/", seas, "/reduced_cdf", seas, ".csv")) %>%
      dplyr::filter(target == tar, location == loc, calendar_week == wk) %>%
      dplyr::arrange(bin_start_incl) %>%
      dplyr::ungroup()
    # turn cdf into matrix
    comp_cdf_matrix <- cdf_frame[, 7:ncol(cdf_frame)] %>%
      as.matrix()
    M <- ncol(comp_cdf_matrix)
    # frame for pdf
    ensemble_pdf <- cdf_frame[, 1:6]
    # build train and test cdf/pdf part
    if (grepl("EW_", ensemble_name, fixed = TRUE)) {
      # ew-bmc or ew-blp
      ## extract params
      omega <- rep(1 / M, M)
      munu <- params[1:(K * 2)]
      w <- params[(1 + (2 * K)):length(params)]
      # equally-weighed component models (ew ensemble)
      H <- comp_cdf_matrix %*% omega
      # build matrices for each k
      pbeta_mat <- matrix(NA, ncol = K, nrow = nrow(comp_cdf_matrix))
      for (k in 1:K) {
        assign(paste0("ab", k), c(munu[k] * munu[k + K], (1 - munu[k]) * munu[k +
                                                                                K]))
        pbeta_mat[, k] <-
          pbeta(H,
                shape1 = get(paste0("ab", k))[1],
                shape2 = get(paste0("ab", k))[2])
      }
    } else {
      # BMC or blp
      ## extract params
      munu <- params[1:(K * 2)]
      omega <- params[1 + (K * 2):(length(params) - (K + 1))]
      w <- params[(1 + (2 * K) + (M * K)):length(params)]
      # build matrices for each k
      pbeta_mat <- matrix(NA, ncol = K, nrow = nrow(comp_cdf_matrix))
      H_mat <- matrix(NA, ncol = K, nrow = nrow(comp_cdf_matrix))
      for (k in 1:K) {
        assign(paste0("ab", k), c(munu[k] * munu[k + K], (1 - munu[k]) * munu[k +
                                                                                K]))
        assign(paste0("omega", k), c())
        for (m in 1:M) {
          assign(paste0("omega", k), c(get(paste0("omega", k)), omega[((m - 1) * K) +
                                                                        k]))
        }
        H_mat[, k] <-
          H <- comp_cdf_matrix %*% c(get(paste0("omega", k)))
        pbeta_mat[, k] <-
          pbeta(H,
                shape1 = get(paste0("ab", k))[1],
                shape2 = get(paste0("ab", k))[2])
      }
    }
    return(list(H_mat, pbeta_mat, w))
  }


### ------------------- retrain TLP ------------------- ###

make_ensemble_emp <- function(component_pdfs, tseason, tar) {
  # generate some inputs
  ncols <- ncol(component_pdfs)
  # make cdfs
  pdf_frame <- component_pdfs %>%
    dplyr::filter(target == tar,
                  season != tseason)
  pdf_m <- as.matrix(pdf_frame[, 9:ncols])
  model_type <- "
    data {
    int<lower=0> n;
    int<lower=0> M;
    matrix[n,M] comp_likelihood;
    }
    parameters {
    simplex[M] omega;
    }

    model {
    vector[M] log_omega = log(omega);
    matrix[n, M] log_likelihood = log(comp_likelihood);
    vector[M] log_weighted_mixt_pdf;
    for (i in 1:n) {
    log_weighted_mixt_pdf = log_omega;
    for (m in 1:M) {
    log_weighted_mixt_pdf[m] += log_likelihood[i, m];
    }
    target += log_sum_exp(log_weighted_mixt_pdf);
    }
    }
    "
  w_vec_inits <-
    sample.int(100, size = ncol(pdf_m), replace = FALSE)
  init_vals <- list(omega = c(w_vec_inits / sum(w_vec_inits)))
  input_data <-
    list(n = nrow(pdf_m),
         M = ncol(pdf_m),
         comp_likelihood = pdf_m)
  
  model_obj <- stan_model(model_code = model_type)
  pars <- rstan::optimizing(
    model_obj,
    data = input_data,
    hessian = T,
    init = init_vals,
    seed = 1234,
    verbose = TRUE
  )
  #make a list of results
  params <- pars$par
  return(params)
}
