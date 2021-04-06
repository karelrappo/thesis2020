source("Import.R")

set.seed(123)

##### List of dependent variables
# H1:H8 - Average quarterly year-on-year growth rates (as in the original paper)
# F1:F8 - Quarterly growth rates of GDP h-quarters ahead
# N1:N8 - Average quarterly growth rates of GDP h-quarters ahead

#### Choose which dependent variable to use for calculating OOS RMSFE-s
#dep <- c("H1", "H2", "H4", "H8")
dep <- c("F1", "F2", "F4", "F8")
#dep <- c("N1", "N2", "N4", "N8")

####Full model independent variables for different regressions
indep_replication <- "YIV + dum + log_gdp + TRM0503 + DGS3MO + SRT03M + VIX + AAA + housng"
indep_replication2 <- c("YIV", "dum", "log_gdp", "TRM0503", "DGS3MO", "SRT03M", "VIX", "AAA", "housng")
indep_RMSFE <- "YIV + dum + TRM1003 + TRM1006 + TRM1012 + TRM0506 + AAA + DBAA + baa_aaa + DGS3MO + SRT03M" 




########################     Replaces p-values with significance stars    ################################################

significance <- function(x){
  for (i in 1:nrow(x))
  { 
    for (k in 1:ncol(x))
    {
      if(grepl(pattern = "p.value", row.names(x)[i])==TRUE)
      {
        if(x[i,k]<0.01)
        { x[i-2,k] <-paste(x[i-2,k],"***")
        } 
        else if(x[i,k]<0.05)
        { x[i-2,k] <-paste(x[i-2,k],"**")
        }
        else if(x[i,k]<0.1)
        { x[i-2,k] <-paste(x[i-2,k],"*")
        } else {
          x[i-2,k] <-x[i-2,k]
        }
      }
    }
  }
  
  x <- x %>%
    rownames_to_column() %>%
    filter(!(grepl(pattern = "p.value",rowname)))
}


###############################################################################################
########################       LINEAR MODELS' RELATED SUMMARY OUTPUT     ######################
###############################################################################################
regrs <- function(indep_vars, interaction)
{
    if(interaction==T){
        df_total <- df %>%
          select(H1, H2, H4, H8, indep_vars) %>%
          gather(Var,Value,  -indep_vars) %>%
          nest(data=c(Value,  indep_vars)) %>%
          mutate(model = map(data, ~lm(Value ~ YIV + dum +YIV*dum, data = .)),
             tidied = map(model, tidy),
             glanced = map(model, glance),
             augmented = map(model, augment),
             neweywest = map(model, ~tidy(coeftest(., vcov.=NeweyWest(., prewhite=FALSE))))) %>%
          select(-model, -data)
  }
  else{
      df_total <- df %>%
        select(H1, H2, H4, H8, indep_vars) %>%
        gather(Var,Value,  -indep_vars) %>%
        nest(data=c(Value,  indep_vars)) %>%
        mutate(model = map(data, ~lm(Value ~ ., data = .)),
             tidied = map(model, tidy),
             glanced = map(model, glance),
             augmented = map(model, augment),
             neweywest = map(model, ~tidy(coeftest(., vcov.=NeweyWest(., prewhite=FALSE))))) %>%
      select(-model, -data)
    }
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
########################     RMSFE calculation function     ###################################
###############################################################################################


out_of_samp <- function(var1, var2, type, var4){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = "final")
  if(var4=="rf"){
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df,
                   method = var4,
                   ntree = 500,
                   tuneGrid = expand.grid(mtry = c(1:8)),
                   trControl = mycontrol,
                   na.action=na.pass)
    dff <- myfit$pred %>%
      arrange(rowIndex)

    dff <- dff %>%
      mutate(Date=lead(df$Date[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(dum=lead(df$dum[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(SE=(pred-obs)^2) %>%
      filter(row_number() <= n()- as.numeric(substring(var1,2,2)))
  }
  else{
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df,
                   method = var4,
                   trControl = mycontrol,
                  na.action=na.pass)
    dff <- myfit$pred %>%
      arrange(rowIndex)  %>%
      mutate(Date=lead(df$Date[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(dum=lead(df$dum[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(SE=(pred-obs)^2) %>%
      filter(row_number() <= n()- as.numeric(substring(var1,2,2)))
  }
  
  
    
    
  ifelse(type=="recessionary", SE <- filter(dff, dum==1),
         ifelse(type=="expansionary", SE <- filter(dff, dum==0),
                SE <- dff))
  SE <- SE$SE
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}



##############################   Calculating OLS RMSFE-s  #########################################

df_results1 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV", "full", "lm")))
df_results2 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV", "recessionary", "lm")))
df_results3 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV", "expansionary", "lm")))
df_results4 <- as.data.frame(t(mapply(out_of_samp, dep, "log_gdp", "full", "lm")))
df_results5 <- as.data.frame(t(mapply(out_of_samp, dep, "log_gdp", "recessionary", "lm")))
df_results6 <- as.data.frame(t(mapply(out_of_samp, dep, "log_gdp", "expansionary", "lm")))
df_results7 <- as.data.frame(t(mapply(out_of_samp, dep, "TRM1012", "full", "lm")))
df_results8 <- as.data.frame(t(mapply(out_of_samp, dep, "baa_aaa", "full", "lm")))
df_results9 <- as.data.frame(t(mapply(out_of_samp, dep, indep_replication, "full", "lm")))
df_results10 <- as.data.frame(t(mapply(out_of_samp, dep, indep_replication,"recessionary", "lm")))
df_results11 <- as.data.frame(t(mapply(out_of_samp, dep, indep_replication, "expansionary", "lm")))



df_resultss <- rbind(df_results1, df_results2, df_results3, df_results4, df_results5, df_results6, df_results7, df_results8, df_results9, df_results10, df_results11) %>%
  round(2)

rownames(df_resultss) <- c("YIV", "YIV_Recessionary", "YIV_Expansionary" ,"Naive", "Naive_Recessionary", "Naive_Expansionary", "TRM", "CRS","full", "full_rec","full_exp")

df_results <- transpose(df_resultss)
df_resultss2 <- melt(df_resultss)
df_resultss2$rowid <- c("YIV", "YIV_Recessionary", "YIV_Expansionary" ,"Naive", "Naive_Recessionary", "Naive_Expansionary", "TRM", "CRS","Full model", "Full - reccessionary", "Full - expansionary")

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

ols <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, indep_RMSFE, "full", "lm")))) %>%
  mutate(Specification="OLS",period="Full sample")
rf <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, indep_RMSFE, "full", "rf")))) %>%
  mutate(Specification="RF",period="Full sample")
ols_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, indep_RMSFE, "recessionary", "lm")))) %>%
  mutate(Specification="OLS", period="Reccessionary")
rf_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, indep_RMSFE, "recessionary", "rf")))) %>%
  mutate(Specification="RF", period="Reccessionary")
ols_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, indep_RMSFE, "expansionary", "lm"))))%>%
  mutate(Specification="OLS", period="Expansionary")
rf_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, indep_RMSFE, "expansionary", "rf"))))%>%
  mutate(Specification="RF",period="Expansionary")

rf_resultsss <- rbind(ols, rf, ols_rec, rf_rec, ols_exp, rf_exp)
rf_resultsss <- rf_resultsss %>%
  pivot_longer(!c(Specification, period,), names_to = "variable", values_to = "value")


###############################################################################################
########################     Predicted and observed values calculation function     ###########
###############################################################################################

out_of_samp2 <- function(var1, var2, var4){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = "final")
  if(var4=="rf"){
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df,
                   method = var4,
                   ntree = 500,
                   tuneGrid = expand.grid(mtry = c(1:8)),
                   trControl = mycontrol,
                   na.action=na.pass)
    dff <- myfit$pred %>%
      arrange(rowIndex) %>%
      mutate(Date=lead(df$Date[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(dum=lead(df$dum[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(SE=(pred-obs)^2) %>%
      filter(row_number() <= n()- as.numeric(substring(var1,2,2)))

  }
  else{
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df,
                   method = var4,
                   trControl = mycontrol,
                   na.action=na.pass)
    dff <- myfit$pred %>%
      arrange(rowIndex) %>%
      mutate(Date=lead(df$Date[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(dum=lead(df$dum[21:103], as.numeric(substring(var1,2,2)))) %>%
      mutate(SE=(pred-obs)^2) %>%
      filter(row_number() <= n()- as.numeric(substring(var1,2,2)))
  }
  
  return(data.frame(predicted=dff$pred, actuals=dff$obs,Date=dff$Date))
}

combiner <- function(type){
  combined <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("predicted", "actuals", "Date")
  colnames(combined) <- x
  for (i in dep){
    output <- out_of_samp2(i,indep_RMSFE,type) %>%
      mutate(variable=i)
      
    combined <- rbind(combined,output)
    
  }
  return(combined)
}

rf <-combiner(type="rf")
lm <-combiner(type="lm")

pred_vs_actual_graph <- function(var,lab){
  ggplot(var, aes(x=as.Date(Date))) + geom_line(aes(y=predicted, color="red")) + geom_line(aes(y=actuals, color="blue"))+ theme_bw() +
    theme(legend.position="bottom",legend.title = element_blank()) + labs(x="Time horizon", y=paste0("GDP predicted with ",lab)) + 
    facet_wrap(~variable) + scale_color_manual(labels = c("Actuals", "Predicted"), values = c("blue", "red")) +
    annotate("rect", xmin = as.Date("2001-04-01", "%Y-%m-%d"), xmax = as.Date("2001-10-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey") +
    annotate("rect", xmin = as.Date("2008-01-01", "%Y-%m-%d"), xmax = as.Date("2009-04-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")
}



###############################################################################################
########################       VARIABLE IMPORTANCE   ##########################################
###############################################################################################

variable_importance <- function(var){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = "final")
  myfit <- train(as.formula(paste0(var, "~", indep_RMSFE )), data = df,
                 method = "rf",
                 ntree = 500,
                 tuneGrid = expand.grid(mtry = c(1:8)),
                 trControl = mycontrol,
                 na.action=na.pass)
  output <- varImp(myfit)
  return(output)
  
}


###############################################################################################
########################          CSSFED            ##########################################
###############################################################################################

df_ts <- ts(df)  
get_statistics <- function(dep, indep_RMSFE, start=1, end=sum(!is.na(df[dep])), est_periods_OOS = 20) {
  
  # Creating vectors where to store values
  #Only use information that is available up to the time at which the forecast is made

  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = "final")
  
    myfit_rf <- train(as.formula(paste0(dep, "~", indep_RMSFE)), data = df,
                   method = "rf",
                   ntree = 500,
                   tuneGrid = expand.grid(mtry = c(1:8)),
                   trControl = mycontrol,
                   na.action=na.pass)

    myfit_lm <- train(as.formula(paste0(dep, "~", indep_RMSFE)), data = df,
                   method = "lm",
                   trControl = mycontrol,
                   na.action=na.pass)

  dff_rf <- myfit_rf$pred %>%
    arrange(rowIndex) %>%
    mutate(Date=lead(df$Date[21:103], as.numeric(substring(dep,2,2)))) %>%
    mutate(dum=lead(df$dum[21:103], as.numeric(substring(dep,2,2)))) %>%
    mutate(SE=(pred-obs)^2) %>%
    filter(row_number() <= n()- as.numeric(substring(dep,2,2)))
 
  OOS_error_rf <- (dff_rf$pred-dff_rf$obs)^2
  
  dff_lm <- myfit_lm$pred %>%
    arrange(rowIndex) %>%
    mutate(Date=lead(df$Date[21:103], as.numeric(substring(dep,2,2)))) %>%
    mutate(dum=lead(df$dum[21:103], as.numeric(substring(dep,2,2)))) %>%
    mutate(SE=(pred-obs)^2) %>%
    filter(row_number() <= n()- as.numeric(substring(dep,2,2)))
  
  OOS_error_lm <- (dff_lm$pred-dff_lm$obs)^2
  
  Date <- dff_lm$Date
  

  
  #### CREATE VARIABLES
  lm_rf <- cumsum(OOS_error_lm)-cumsum(OOS_error_rf)

  
  return(as.tibble(list(OOS_error_lm = OOS_error_lm,
                        OOS_error_rf = OOS_error_rf,
                        lm_rf = lm_rf,
                        Date=Date)))
}
#CSSFED results combine & create graph function
CSSFED_all <- function(dependent){
  output_combined <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("OOS_error_lm","OOS_error_rf","lm_rf", "Date")
  colnames(output_combined) <- x
  for (i in dependent){
  output <- get_statistics(i, indep_RMSFE) %>%
    mutate(Dependent=i)

   output_combined <- rbind(output_combined, output)
  }
  return(output_combined)
}
CSSFED <- CSSFED_all(dep) 
CSSFED <- pivot_longer(CSSFED, cols=!c(Date,Dependent), names_to = "type", values_to = "values")

squared_error_plot <- function(var, dependent){
  
  squared_errors <- CSSFED %>%
    filter(type %in% var)
  
  plot <- ggplot(data = squared_errors, aes(x=Date,y=values)) + geom_line(aes(color=factor(type))) + 
    annotate("rect", xmin = as.Date("2001-04-01", "%Y-%m-%d"), xmax = as.Date("2001-10-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    annotate("rect", xmin = as.Date("2008-01-01", "%Y-%m-%d"), xmax = as.Date("2009-04-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    xlab("Time horizon") + ylab(paste0("Squared Errors")) + theme_bw() + facet_wrap(~Dependent)
  return(plot)
}


CSSFED_plot <- function(var){
  
  cssfed_values <- CSSFED %>%
    filter(type == var)
  
  plot <- ggplot(data = cssfed_values, aes(x=Date,y=values)) + geom_point(colour='red') + 
    annotate("rect", xmin = as.Date("2001-04-01", "%Y-%m-%d"), xmax = as.Date("2001-10-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    annotate("rect", xmin = as.Date("2008-01-01", "%Y-%m-%d"), xmax = as.Date("2009-04-01",  "%Y-%m-%d"), ymin = -Inf, ymax = Inf,alpha = 0.4, fill = "grey")+
    xlab("Time horizon") +
    ylab(paste0("CSSFED for ",var)) + theme_bw() + facet_wrap(~Dependent) 
  
  return(plot)


}


