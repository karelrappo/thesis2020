source("Import.R")
library("caret")
library("dplyr")
library("rollRegres")
library("mltools")

out_of_samp <- function(var1, var2, var3, var4){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = TRUE)
  myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                 method = var4,
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

# SRT03M + gz_spr  + spy_logreturn

ols <- as.data.frame(as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng", 4, "lm"))))
rf <- as.data.frame(as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng", 4, "rf"))))

ols_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng", 5, "lm"))))
rf_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng", 5, "rf"))))

ols_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng", 6, "lm"))))
rf_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng", 6, "rf"))))


rf_resultsss <- rbind(ols, rf, ols_rec, rf_rec, ols_exp, rf_exp)
colnames(rf_resultsss, do.NULL = FALSE)
colnames(rf_resultsss) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
rownames(rf_resultsss) <- c("ols", "rf", "ols_rec", "rf_rec", "ols_exp", "rf_exp")


results <- function(var1){
  rf_resultss2 <- as.data.frame(t(rf_resultsss[var1,]))
  rf_resultss2$time <- c(1:12)
  rf_resultss2 <- melt(rf_resultss2 ,  id.vars = 'time', variable.name = 'series') 
  ggplot(rf_resultss2, aes((time), value)) + geom_line(aes(colour = series)) + theme_bw() + theme(legend.position = "bottom") + labs(x="Time horizon", y="RMSFE value", color="Model")
  
}

plot1 <- results(c("ols", "rf", "ols_rec", "rf_rec", "ols_exp", "rf_exp"))
plot1



#rf_results1 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 4)))
#rf_results2 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 5)))
#rf_results3 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 6)))
#rf_results4 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "log_gdp", 4)))
#rf_results5 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "TRM0503 + TRM0506 + TRM1003 + TRM1006 + TRM1012", 4)))
#rf_results6 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "AAA + DBAA + baa_aaa", 4)))



