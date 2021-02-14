source("Import.R")
library(rollRegres)
library(mltools)


# Out-of-sample regression
do_regression <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ YIV")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  RMSFE <- rmse(preds = predicted, actuals = actual)
  return(RMSFE)
}


# Out-of-sample regression
do_regression2 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ YIV")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  SE <- (predicted-actual)^2
  SE <- SE[c(25,26,27,52,53,54,55,56,57)]
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}


# Out-of-sample regression
do_regression3 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ YIV")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  SE <- (predicted-actual)^2
  SE <- SE[-c(25,26,27,52,53,54,55,56,57)]
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}


# Out-of-sample regression
do_regression4 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ log_gdp")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  RMSFE <- rmse(preds = predicted, actuals = actual)
  return(RMSFE)
}

# Out-of-sample regression
do_regression5 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ log_gdp")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  SE <- (predicted-actual)^2
  SE <- SE[c(25,26,27,52,53,54,55,56,57)]
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}


# Out-of-sample regression
do_regression6 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ log_gdp")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  SE <- (predicted-actual)^2
  SE <- SE[-c(25,26,27,52,53,54,55,56,57)]
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}


# TRM1006 doesn't work in the regression for some reason
do_regression7 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ TRM1012")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  RMSFE <- rmse(preds = predicted, actuals = actual)
  return(RMSFE)
}


do_regression8 <- function(var)
{ # var - the name of the variable to be regressed
  res <-roll_regres(as.formula(paste0(var, "~ baa_aaa")), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  RMSFE <- rmse(preds = predicted, actuals = actual)
  return(RMSFE)
}

df_results1 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression)))
df_results2 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression2)))
df_results3 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression3)))
df_results4 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression4)))
df_results5 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression5)))
df_results6 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression6)))
df_results7 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression7)))
df_results8 <- as.data.frame(t(sapply(colnames(df)[29:40], do_regression8)))


df_resultss <- rbind(df_results1, df_results2, df_results3, df_results4, df_results5, df_results6, df_results7, df_results8) %>%
  round(2)

colnames(df_resultss, do.NULL = FALSE)
colnames(df_resultss) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
rownames(df_resultss) <- c("YIV", "YIV_Recessionary", "YIV_Expansionary" ,"Naive", "Naive_Recessionary", "Naive_Expansionary", "TRM", "CRS")


relative_rmsfes <- df_resultss %>%
  t() %>%
  as_tibble() %>%
  mutate(rRMSE_recess = YIV_Recessionary/YIV) %>%
  mutate(rRMSE_expans = YIV_Expansionary/YIV) %>%
  select(rRMSE_recess, rRMSE_expans) %>%
  round(2) %>%
  t()

colnames(relative_rmsfes) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")








