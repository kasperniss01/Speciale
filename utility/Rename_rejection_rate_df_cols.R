rm(list = ls())

### Modifying naming conventions

naming_lookup <- c("rate_parametric_plugin" = "rates_parametric_plugin",
               "rate_parametric" = "rates_parametric",
               "rate_nonparametric" = "rate",
               "se_nonparametric" = "se")



for(file in   list.files("datasets") ){
  print(file)
  
  rds_obj <- readRDS(file.path("datasets", file))
  
  if(is.null(rds_obj$rejection_rate_df)) next
  
  rds_obj$rejection_rate_df <- rds_obj$rejection_rate_df %>% rename(any_of(naming_lookup))
  
  saveRDS(rds_obj, file.path("datasets", file))

}
