source("Import.R")
library("caret")
library("dplyr")

out_of_samp <- function(var1, var2, var3){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = TRUE)
  myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                 method = "rf",
                 ntree = 50,
                 trControl = mycontrol)
  dff <- myfit$pred
  SE <- (dff$pred-dff$obs)^2
  if (var3 == 5){
    SE <- SE[c(25,26,27,52,53,54,55,56,57)]
  }
  if (var3 == 6){
    SE <- SE[-c(25,26,27,52,53,54,55,56,57)]
  }
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
} 


rf_results1 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 4)))
rf_results2 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 5)))
rf_results3 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 6)))
#rf_results4 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "log_gdp", 4)))
#rf_results5 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "TRM0503 + TRM0506 + TRM1003 + TRM1006 + TRM1012", 4)))
#rf_results6 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "AAA + DBAA + baa_aaa", 4)))


rf_resultss <- rbind(rf_results1, rf_results2, rf_results3)

rf_resultss <- rf_resultss %>%
  round(2)

colnames(rf_resultss, do.NULL = FALSE)
colnames(rf_resultss) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
rownames(rf_resultss) <- c("RF_YIV", "RF_Recessionary", "RF_Expansionary")








