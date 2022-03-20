# load in saved params and list
n_list <- ls()
check <- c(lapply(n_list, function(x)
  get(x)))
names(check) <- n_list

save(check, file = "ensemble_code/train_params/EW_BLP18.rda")
# delete environment
n_list <- ls()
check <- c(lapply(n_list, function(x)
  get(x)))
names(check) <- n_list

save(check, file = "ensemble_code/train_params/EW_BLP19.rda")


#######---------------------- BLP------------_#####################
n_list <- ls()
check <- c(lapply(n_list, function(x)
  get(x)))
names(check) <- n_list

save(check, file = "ensemble_code/train_params/BLP17.rda")

n_list <- ls()
check <- c(lapply(n_list, function(x)
  get(x)))
names(check) <- n_list

save(check, file = "ensemble_code/train_params/BLP19.rda")



#######---------------------- final train-pars------------_#####################
# BMC2 and EW_BMC2 for  17 for all targets
for (i in 1:4) {
  assign(
    paste0("target", i, "_sea17_params"),
    list(
      BLP = BLP17[[1]][[i]],
      BMC2 = BMC2_17[[i]],
      EW_BLP = EW_BLP17[[1]][[i]],
      EW_BMC2 = EW_BMC2_17[[i]]
    )
  )
}

save(target1_sea17_params, file = "ensemble_code/train_params/final_train/target1_sea17_params.rda")
save(target2_sea17_params, file = "ensemble_code/train_params/final_train/target2_sea17_params.rda")
save(target3_sea17_params, file = "ensemble_code/train_params/final_train/target3_sea17_params.rda")
save(target4_sea17_params, file = "ensemble_code/train_params/final_train/target4_sea17_params.rda")

# BMC2 and EW_BMC2 for  18 for all targets
for (i in 1:4) {
  assign(
    paste0("target", i, "_sea18_params"),
    list(
      BLP = BLP18[[i]],
      BMC2 = BMC2_18[[i]],
      EW_BLP = EW_BLP18[[i]],
      EW_BMC2 = EW_BMC2_18[[i]]
    )
  )
}
target4_sea18_params <-
  list(
    BLP = BLP18[[4]],
    BMC2 = BMC2_18[[4]],
    EW_BLP = EW_BLP18[[4]],
    EW_BMC2 = EW_BMC2_18[[4]]
  )

save(target1_sea18_params, file = "ensemble_code/train_params/final_train/target1_sea18_params.rda")
save(target2_sea18_params, file = "ensemble_code/train_params/final_train/target2_sea18_params.rda")
save(target3_sea18_params, file = "ensemble_code/train_params/final_train/target3_sea18_params.rda")
save(target4_sea18_params, file = "ensemble_code/train_params/final_train/target4_sea18_params.rda")

# BMC2 and EW_BMC2 for  19 for all targets but EW-BMC5 for 3 wk ahead

for (i in c(1:4)) {
  assign(
    paste0("target", i, "_sea19_params"),
    list(
      BLP = BLP19[[i]],
      BMC2 = BMC2_19[[i]],
      EW_BLP = EW_BLP19[[i]],
      EW_BMC2 = EW_BMC2_19[[i]]
    )
  )
}
target3_sea19_params <-
  list(
    BLP = BLP19[[3]],
    BMC2 = BMC2_19[[3]],
    EW_BLP = EW_BLP19[[3]],
    EW_BMC2 = EW_BMC2_19[[3]]
  )

save(target1_sea19_params, file = "ensemble_code/train_params/final_train/target1_sea19_params.rda")
save(target2_sea19_params, file = "ensemble_code/train_params/final_train/target2_sea19_params.rda")
save(target3_sea19_params, file = "ensemble_code/train_params/final_train/target3_sea19_params.rda")
save(target4_sea19_params, file = "ensemble_code/train_params/final_train/target4_sea19_params.rda")

