source("Import.R")


###############################################################################################
########################       LINEAR MODELS   ################################################
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
########################       OOS FORECASTS   ################################################
###############################################################################################

#var - dependent variable
#var2 - independent variable

out_of_sample <- function(var,var2,type)
{ 
  res <-roll_regres(as.formula(paste0(var,"~",var2)), df[1:91,], width=20L, do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  predicted <- unlist(res[4])
  predicted <- predicted[21:91]
  actual <- df[[var]][21:91]
  SE <- (predicted-actual)^2
  
  ifelse(type=="recessionary",SE <- SE[c(25,26,27,52,53,54,55,56,57)],
  ifelse(type=="expansionary", SE <- SE[-c(25,26,27,52,53,54,55,56,57)],
          SE <- SE))
  
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}


##############################   FORMAT OOS RESULTS  ##############################################################


df_results1 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "YIV", "full")))
df_results2 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "YIV", "recessionary")))
df_results3 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "YIV", "expansionary")))
df_results4 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "log_gdp", "full")))
df_results5 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "log_gdp", "recessionary")))
df_results6 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "log_gdp", "expansionary")))
df_results7 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "TRM1012", "full")))
df_results8 <- as.data.frame(t(mapply(out_of_sample, colnames(df)[29:40], "baa_aaa", "full")))

df_resultss <- rbind(df_results1, df_results2, df_results3, df_results4, df_results5, df_results6, df_results7, df_results8) %>%
  round(2)

colnames(df_resultss, do.NULL = FALSE)
colnames(df_resultss) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
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

colnames(relative_rmsfes) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
rownames(df_resultss)[rownames(df_resultss)=='YIV_Recessionary'] <- "YIV-Recess."
rownames(df_resultss)[rownames(df_resultss)=='YIV_Expansionary'] <- "YIV-Expans."
rownames(df_resultss)[rownames(df_resultss)=='Naive_Recessionary'] <- "Naive-Recess."
rownames(df_resultss)[rownames(df_resultss)=='Naive_Expansionary'] <- "Naive-Expans."


###############################################################################################
########################     RF FULL MODEL     ################################################
###############################################################################################

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


###############################################################################################
########################         MRF           ################################################
###############################################################################################
MRF.data <- df %>%
  select(-Date, -log_gdp) %>%
  relocate(GDP, YIV, dum)%>% 
  as.matrix()

mrf.output=MRF(data=MRF.data,y.pos=1,x.pos=2:5,S.pos=2:20,oos.pos=90:103)

predicted <- as_tibble(mrf.output$pred)
predicted <- predicted %>%
  mutate(GDP = df[90:103,]$GDP)

RMSE<-mltools::rmse(preds=predicted$value,actuals=predicted$GDP)
