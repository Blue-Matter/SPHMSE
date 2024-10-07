## code to prepare `DATASET` dataset goes here

RCM_hake_2023 <- readRDS("data-raw/RCM_hake_2023.rds")
usethis::use_data(RCM_hake_2023, overwrite = TRUE)
