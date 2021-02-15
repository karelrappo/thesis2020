source("Import.R")
library("caret")


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


df_results1 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 4)))
df_results2 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 5)))
df_results3 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "YIV", 6)))
#df_results4 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "log_gdp", 4)))
#df_results5 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "TRM0503 + TRM0506 + TRM1003 + TRM1006 + TRM1012", 4)))
#df_results6 <- as.data.frame(t(mapply(out_of_samp, colnames(df)[29:40], "AAA + DBAA + baa_aaa", 4)))


df_resultss <- rbind(df_results1, df_results2, df_results3)

df_resultss <- df_resultss %>%
  round(2)

colnames(df_resultss, do.NULL = FALSE)
colnames(df_resultss) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
rownames(df_resultss) <- c("RF_YIV", "RF_Recessionary", "RF_Expansionary")

ols_rmsfe <- read.csv("Data/RMSFE comparison.csv")
ols_rmsfe <- column_to_rownames(ols_rmsfe, "X")
df_resultss1 <- rbind(df_resultss, ols_rmsfe[1:3,])

df_resultss2 <- as.data.frame(t(df_resultss1[c(1,4),]))
df_resultss2$time <- c(1:12)
df_resultss2 <- melt(df_resultss2 ,  id.vars = 'time', variable.name = 'series')
ggplot(df_resultss2, aes(as.numeric(time), value)) + geom_line(aes(colour = series))

df_resultss3 <- as.data.frame(t(df_resultss1[c(2,5),]))
df_resultss3$time <- c(1:12)
df_resultss3 <- melt(df_resultss3 ,  id.vars = 'time', variable.name = 'series')
ggplot(df_resultss3, aes(time, value)) + geom_line(aes(colour = series))

df_resultss4 <- as.data.frame(t(df_resultss1[c(3,6),]))
df_resultss4$time <- c(1:12)
df_resultss4 <- melt(df_resultss4 ,  id.vars = 'time', variable.name = 'series')
ggplot(df_resultss4, aes(time, value)) + geom_line(aes(colour = series))






