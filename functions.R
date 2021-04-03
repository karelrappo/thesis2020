source("Import.R")

######################## Select dataset, df by default ###########################################

set.seed(123)
##### List of dependent variables
# H1:H8 - Average quarerly year-on-year growth rates (as in the original paper)
# F1:F8 - Quarterly growth rates of GDP h-quarters ahead
# N1:N8 - Average quarterly growth rates of GDP h-quarters ahead

#dep <- c("H1", "H2", "H4", "H8")
dep <- c("F1", "F2", "F4", "F8")
#dep <- c("N1", "N2", "N4", "N8")


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
########################       LINEAR MODELS' RELATED SUMMARY STATISTICS   ####################
###############################################################################################
regrs <- function(mudelid, interaction)
{
    if(interaction==T){
        df_total <- df %>%
          select(H1, H2, H4, H8, mudelid) %>%
          gather(Var,Value,  -mudelid) %>%
          nest(data=c(Value,  mudelid)) %>%
          mutate(model = map(data, ~lm(Value ~ YIV + dum +YIV*dum, data = .)),
             tidied = map(model, tidy),
             glanced = map(model, glance),
             augmented = map(model, augment),
             neweywest = map(model, ~tidy(coeftest(., vcov.=NeweyWest(., prewhite=FALSE))))) %>%
          select(-model, -data)
  }
  else{
      df_total <- df %>%
        select(H1, H2, H4, H8, mudelid) %>%
        gather(Var,Value,  -mudelid) %>%
        nest(data=c(Value,  mudelid)) %>%
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
########################     RMSFE calculation function     ################################################
###############################################################################################

out_of_samp <- function(var1, var2, type, var4){
  mycontrol <- trainControl(method = "timeslice",
                            initialWindow = 20,
                            horizon = 1,
                            fixedWindow = TRUE, 
                            savePredictions = "final")
  if(var4=="rf"){
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   ntree = 500,
                   tuneGrid = expand.grid(mtry = c(1:8)),
                   trControl = mycontrol)
    dff <- myfit$pred %>%
      arrange(rowIndex)

    dff <- dff %>%
      mutate(dum=df$dum[21:103]) %>%
      mutate(SE=(pred-obs)^2)
  }
  else{
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   trControl = mycontrol)
    dff <- myfit$pred %>%
      arrange(rowIndex)  %>%
      mutate(dum=df$dum[21:103]) %>%
      mutate(SE=(pred-obs)^2)
  }
  
  
  ifelse(type=="recessionary", SE <- filter(dff, dum==1),
         ifelse(type=="expansionary", SE <- filter(dff, dum==0),
                SE <- dff))
  SE <- SE$SE
  RMSFE <- sqrt(mean(SE))
  return(RMSFE)
}



##############################   OLS RMSFE-s  ##############################################################

df_results1 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV", "full", "lm")))
df_results2 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV", "recessionary", "lm")))
df_results3 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV", "expansionary", "lm")))
df_results4 <- as.data.frame(t(mapply(out_of_samp, dep, "log_gdp", "full", "lm")))
df_results5 <- as.data.frame(t(mapply(out_of_samp, dep, "log_gdp", "recessionary", "lm")))
df_results6 <- as.data.frame(t(mapply(out_of_samp, dep, "log_gdp", "expansionary", "lm")))
df_results7 <- as.data.frame(t(mapply(out_of_samp, dep, "TRM1012", "full", "lm")))
df_results8 <- as.data.frame(t(mapply(out_of_samp, dep, "baa_aaa", "full", "lm")))
df_results9 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa + VIX + housng + SRT03M", "full", "lm")))
df_results10 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa + VIX + housng + SRT03M", "recessionary", "lm")))
df_results11 <- as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa + VIX + housng + SRT03M", "expansionary", "lm")))



df_resultss <- rbind(df_results1, df_results2, df_results3, df_results4, df_results5, df_results6, df_results7, df_results8, df_results9, df_results10, df_results11) %>%
  round(2)

rownames(df_resultss) <- c("YIV", "YIV_Recessionary", "YIV_Expansionary" ,"Naive", "Naive_Recessionary", "Naive_Expansionary", "TRM", "CRS","full", "full_rec","full_exp")

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

ols <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "full", "lm")))) %>%
  mutate(Specification="OLS",period="Full sample")
rf <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "full", "rf")))) %>%
  mutate(Specification="RF",period="Full sample")
ols_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "recessionary", "lm")))) %>%
  mutate(Specification="OLS", period="Reccessionary")
rf_rec <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "recessionary", "rf")))) %>%
  mutate(Specification="RF", period="Reccessionary")
ols_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "expansionary", "lm"))))%>%
  mutate(Specification="OLS", period="Expansionary")
rf_exp <- as.data.frame(as.data.frame(t(mapply(out_of_samp, dep, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M", "expansionary", "rf"))))%>%
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
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   ntree = 500,
                   tuneGrid = expand.grid(mtry = c(1:8)),
                   trControl = mycontrol)
    dff <- myfit$pred %>%
      arrange(rowIndex)

  }
  else{
    myfit <- train(as.formula(paste0(var1, "~", var2)), data = df[1:sum(!is.na(df[var1])),],
                   method = var4,
                   trControl = mycontrol)
    dff <- myfit$pred %>%
      arrange(rowIndex)
  }
  
  return(data.frame(predicted=dff$pred, actuals=dff$obs,Date=df$Date[21:103]))
}
type <- "lm"
i <- "H4"
combiner <- function(type){
  combined <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("predicted", "actuals", "Date")
  colnames(combined) <- x
  for (i in dep){
    output <- out_of_samp2(i,"YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M",type) %>%
      mutate(variable=i) %>%
      mutate(Date=lead(df$Date[21:103], as.numeric(substring(i,2,2)))) %>%
      filter(row_number() <= n()- as.numeric(substring(i,2,2)))
      
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
  myfit <- train(as.formula(paste0(var, "~ YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M")), data = df[1:sum(!is.na(df[var])),],
                 method = "rf",
                 ntree = 500,
                 tuneGrid = expand.grid(mtry = c(1:8)),
                 trControl = mycontrol)
  output <- varImp(myfit)
  return(output)
  
}


###############################################################################################
########################          CSSFED            ##########################################
###############################################################################################

df_ts <- ts(df)  
get_statistics <- function(dep, indep, start=1, end=sum(!is.na(df[dep])), est_periods_OOS = 20) {
  
  # Creating vectors where to store values
  OOS_error_hist <- numeric(end - start - est_periods_OOS)
  #Only use information that is available up to the time at which the forecast is made
  j <- 0
    for (i in 21:(end)) {
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
                            savePredictions = "final")
  
    myfit_rf <- train(as.formula(paste0(dep, "~", indep)), data = df[1:sum(!is.na(df[dep])),],
                   method = "rf",
                   ntree = 500,
                   tuneGrid = expand.grid(mtry = c(1:8)),
                   trControl = mycontrol)

    myfit_lm <- train(as.formula(paste0(dep, "~", indep)), data = df[1:sum(!is.na(df[dep])),],
                   method = "lm",
                   trControl = mycontrol)

  dff_rf <- myfit_rf$pred %>%
    arrange(rowIndex)
 
  OOS_error_rf <- (dff_rf$pred-dff_rf$obs)^2
  
  dff_lm <- myfit_lm$pred %>%
    arrange(rowIndex)
  
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
#CSSFED results combine & create graph function
CSSFED_all <- function(dependent){
  output_combined <- data.frame(matrix(ncol = 7, nrow = 0))
  x <- c("OOS_error_hist","OOS_error_lm","OOS_error_rf", "hist_lm","lm_rf","hist_rf", "Date")
  colnames(output_combined) <- x
  for (i in dependent){
  output <- get_statistics(i, "YIV + dum + DGS1 + TRM1012 + baa_aaa+ VIX + housng + SRT03M") %>%
    mutate(Date=lead(df$Date[21:103], as.numeric(substring(i,2,2)))) %>%
    mutate(Dependent=i) %>% 
    filter(row_number() <= n()- as.numeric(substring(i,2,2)))
  
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


