source("Import.R")
library(dplyr)
library(reshape2)
library(stringr)
library(tidyverse)
library(summarytools)
library(kableExtra)
library(dynlm)
library(stats)
library(readr)
library(ggplot2)
library(MacroRF)
require(data.table)
library(zoo)
library(sandwich)
library(randomForest)
library(caret)
library(fbi)
library(lmtest)
library(purrr)
library(broom)
library(rollRegres)
library(modelr)
library(forcats)
library(tibble)
library(Metrics)
library(mltools)
library(reshape2)
library(lmtest)
library(olsrr)
library(pracma)
library(dyn)
library(gridExtra)

###############################################################################################
########################       LINEAR MODELS' RELATED SUMMARY STATISTICS   ####################
###############################################################################################
regrs <- function(mudelid)
{
  df_total <- df %>%
    select(H4, H6, H8, H10, H12, mudelid) %>%
    gather(Var,Value,  -mudelid) %>%
    nest(data=c(Value,  mudelid)) %>%
    mutate(model = map(data, ~lm(Value ~ ., data = .)),
           tidied = map(model, tidy),
           glanced = map(model, glance),
           augmented = map(model, augment),
           neweywest = map(model, ~tidy(coeftest(., vcov.=NeweyWest(., prewhite=FALSE))))) %>%
    select(-model, -data)
  return(df_total)
}


regr_results <- function(a){
  results1 <- a %>%
    select(-augmented,-tidied, -glanced) %>%
    unnest(cols = c(neweywest)) %>%
    select(-statistic) %>%
    pivot_longer(cols=-c(1:2), names_to = "mdeea", values_to = "Delay") %>%
    mutate(x = paste(term,mdeea,sep = "_"))%>%
    filter(!grepl('Intercept', term)) %>%
    select(-c(2:3)) %>%
    pivot_wider(names_from = Var, values_from = Delay) %>%
    column_to_rownames(var = "x")
  
results2 <- a %>%
    select(-augmented,-neweywest, -tidied) %>%
    unnest(cols = c(glanced)) %>%
    select(-df, -AIC, -BIC, -deviance, -nobs, -df.residual, -logLik, -statistic, -p.value, -sigma) %>%
    column_to_rownames(var = "Var")

results2 <- as.data.frame(t(results2))
  
results <- rbind(results1, results2)
remove(results1, results2)
  
  return(results)
}


###############################################################################################
########################     RMSFE calculation function     ################################################
###############################################################################################

out_of_samp <- function(var1, var2, type, var4){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = TRUE)
  if(var4=="rf"){
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   ntree = 50,
                   tuneGrid=data.frame(mtry=4),
                   trControl = mycontrol)
  }
  else{
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   ntree = 50,
                   trControl = mycontrol)
  }
  dff <- myfit$pred
  SE <- (dff$pred-dff$obs)^2
  
  ifelse(type=="recessionary",SE <- SE[c(25,26,27,52,53,54,55,56,57)],
         ifelse(type=="expansionary", SE <- SE[-c(25,26,27,52,53,54,55,56,57)],
                SE <- SE))
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}

##############################   OLS RMSFE-s  ##############################################################

df_results1 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV", "full", "lm")))
df_results2 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV", "recessionary", "lm")))
df_results3 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV", "expansionary", "lm")))
df_results4 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "log_gdp", "full", "lm")))
df_results5 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "log_gdp", "recessionary", "lm")))
df_results6 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "log_gdp", "expansionary", "lm")))
df_results7 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "TRM1012", "full", "lm")))
df_results8 <- as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "baa_aaa", "full", "lm")))

df_resultss <- rbind(df_results1, df_results2, df_results3, df_results4, df_results5, df_results6, df_results7, df_results8) %>%
  round(2)

colnames(df_resultss, do.NULL = FALSE)
colnames(df_resultss) <- c("H1","H2","H4","H8")
rownames(df_resultss) <- c("YIV", "YIV_Recessionary", "YIV_Expansionary" ,"Naive", "Naive_Recessionary", "Naive_Expansionary", "TRM", "CRS")

##############################   RELATIVE RMSFE  ##############################################################

relative_rmsfes <- df_resultss %>%
  t() %>%
  as_tibble() %>%
  mutate(rRMSE_recess = YIV_Recessionary/YIV) %>%
  mutate(rRMSE_expans = YIV_Expansionary/YIV) %>%
  select(rRMSE_recess, rRMSE_expans) %>%
  round(2) %>%
  t()

colnames(relative_rmsfes) <- c("H1","H2","H4","H8")
rownames(df_resultss)[rownames(df_resultss)=='YIV_Recessionary'] <- "YIV-Recess."
rownames(df_resultss)[rownames(df_resultss)=='YIV_Expansionary'] <- "YIV-Expans."
rownames(df_resultss)[rownames(df_resultss)=='Naive_Recessionary'] <- "Naive-Recess."
rownames(df_resultss)[rownames(df_resultss)=='Naive_Expansionary'] <- "Naive-Expans."


##############################   RF & OLS RMSFE results' comparison ##############################################################

ols <- as.data.frame(as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "full", "lm"))))
rf <- as.data.frame(as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "full", "rf"))))

ols_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "recessionary", "lm"))))
rf_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "recessionary", "rf"))))

ols_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "expansionary", "lm"))))
rf_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, c("F1", "F2", "F4", "F8"), "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "expansionary", "rf"))))


rf_resultsss <- rbind(ols, rf, ols_rec, rf_rec, ols_exp, rf_exp)
colnames(rf_resultsss, do.NULL = FALSE)
colnames(rf_resultsss) <- c("H1","H2","H4","H8")
rf_resultsss$period <- c("Linear model", "Random forest")
rf_resultsss$type <- c("Full sample","Full sample","Reccessionary period","Reccessionary period","Expansionary period","Expansionary period")
rf_resultsss$type2 <- c(1,2)
rf_resultsss <- rf_resultsss %>%
  pivot_longer(!c(type, period, type2), names_to = "variable", values_to = "value")


###############################################################################################
########################     Predicted and observed values calculation function     ###########
###############################################################################################

out_of_samp2 <- function(var1, var2, var4){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = TRUE)
  if(var4=="rf"){
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   ntree = 50,
                   tuneGrid=data.frame(mtry=4),
                   trControl = mycontrol)
  }
  else{
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   ntree = 50,
                   trControl = mycontrol)
  }
  dff <- myfit$pred
  return(list(dff$pred, dff$obs))
} 

act_vs_predicted <- function(var1, var2){
  act_vs_predicted <- mapply(out_of_samp2, c("F1", "F2", "F4", "F8"), "YIV", var1)
  dfff1 <- data.frame(matrix(unlist(act_vs_predicted), nrow=8, byrow=TRUE),stringsAsFactors=FALSE)
  colnames(dfff1) <- c(1:80)
  dfff1$type1 <- c("H1","H1", "H2","H2", "H4","H4", "H8", "H8")
  dfff1$type2 <- c(var2, "Actual")
  dfff1 <- dfff1 %>%
    pivot_longer(!c(type1, type2), names_to = "variable", values_to = "value")
}

lm <- act_vs_predicted("lm", "Predicted (Linear model)")
rf <- act_vs_predicted("rf", "Predicted (Random forest)")

pred_vs_actual_graph <- function(var){
  ggplot(var, aes(variable, value, group=factor(type2))) + geom_line(aes(color=factor(type2))) + theme_bw() + 
    theme(legend.position="bottom") + labs(x="Time horizon", y="GDP growth value", color="") + 
    facet_wrap(~type1) + scale_color_manual(values=c("blue", "red"))
}

###############################################################################################
########################       VARIABLE IMPORTANCE   ##########################################
###############################################################################################

variable_importance <- function(var){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = TRUE)
  myfit <- train(as.formula(paste0(var, "~ YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M")), data = df[1:sum(!is.na(df[var])),],
                 method = "rf",
                 ntree = 50,
                 tuneGrid = data.frame(mtry = 4),
                 trControl = mycontrol)
  output <- varImp(myfit)
  return(output)
  
}


###############################################################################################
########################          CSSFED            ##########################################
###############################################################################################

df_ts <- ts(df)  
get_statistics <- function(dep, indep, start=1, end=103, est_periods_OOS = 20) {
  
  # Creating vectors where to store values
  OOS_error_hist <- numeric(end - start - est_periods_OOS)
  #Only use information that is available up to the time at which the forecast is made
  j <- 0
    for (i in 21:(end-1)) {
      j <- j + 1
      #Get the actual ERP that you want to predict
      actual <- as.numeric(window(df_ts, i, i)[, dep])
      
      #1. Historical mean model
      OOS_error_hist[j] <- (actual - mean(window(df_ts, start, i-1)[, dep], na.rm=TRUE))^2
      
    }

  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = TRUE)
  
    myfit_rf <- train(as.formula(paste0(dep, "~", indep)), data = df[1:sum(!is.na(df[dep])),],
                   method = "rf",
                   ntree = 50,
                   tuneGrid=data.frame(mtry=4),
                   trControl = mycontrol)

    myfit_lm <- train(as.formula(paste0(dep, "~", indep)), data = df[1:sum(!is.na(df[dep])),],
                   method = "lm",
                   ntree = 50,
                   trControl = mycontrol)


  dff_rf <- myfit_rf$pred
  OOS_error_rf <- (dff_rf$pred-dff_rf$obs)^2
  
  dff_lm <- myfit_lm$pred
  OOS_error_lm <- (dff_lm$pred-dff_lm$obs)^2
  


  
  #### CREATE VARIABLES
  hist_lm <- cumsum(OOS_error_hist)-cumsum(OOS_error_lm)
  lm_rf <- cumsum(OOS_error_lm)-cumsum(OOS_error_rf)
  hist_rf <- cumsum(OOS_error_hist)-cumsum(OOS_error_rf)
  
  
  return(as.tibble(list(OOS_error_hist = OOS_error_hist,
                        OOS_error_lm = OOS_error_lm,
                        OOS_error_rf = OOS_error_rf,
                        hist_lm = hist_lm,
                        lm_rf = lm_rf,
                        hist_rf = hist_rf)))
}


CSSFED <- get_statistics( "F1", "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M")
CSSFED$Date <- df$Date[21:102]
CSSFED <- pivot_longer(CSSFED, cols = -c(7), names_to = "type", values_to = "values")


squared_error_plot <- function(var){

  squared_errors <- CSSFED %>%
    filter(type == var)

  plot <- ggplot(data = squared_errors, aes(x=Date,y=c(values))) + geom_line(aes(color=factor(type))) + 
    annotate("rect", xmin = as.Date("2001-04-01", "%Y-%m-%d"), xmax = as.Date("2001-10-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    annotate("rect", xmin = as.Date("2008-01-01", "%Y-%m-%d"), xmax = as.Date("2009-04-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    xlab("Time horizon") + ylab("CSSFED")
  
  return(plot)
}


CSSFED_plot <- function(var){
  
  cssfed_values <- CSSFED %>%
    filter(type == var)
  
  plot <- ggplot(data = cssfed_values, aes(x=Date,y=values)) + geom_point(colour='red') + 
    annotate("rect", xmin = as.Date("2001-04-01", "%Y-%m-%d"), xmax = as.Date("2001-10-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    annotate("rect", xmin = as.Date("2008-01-01", "%Y-%m-%d"), xmax = as.Date("2009-04-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    xlab("Time horizon") + ylab("CSSFED")
  
  return(plot)
}

