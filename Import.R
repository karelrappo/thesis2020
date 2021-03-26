library(stringr)
library(tidyverse)
library(summarytools)
library(kableExtra)
library(dynlm)
library(stats)
library(readr)
library(dplyr)
library(ggplot2)
library(data.table)
library(MacroRF)
require(data.table)
library(zoo)
library(rollRegres)
library(mltools)
library(caret)


start <- as.Date("01/03/1990", format="%d/%m/%Y")
end <- as.Date("01/10/2015", format="%d/%m/%Y")

#Import extracted GDP and YIV
dataset <- read_csv("data/gdp_yiv.csv")
dataset$Date <- as.Date(dataset$Date, format="%d.%m.%Y")


#Import Housing data and convert to quarterly
housng <- read_csv("data/HOUSNG.csv")
housng$DATE <- as.Date(housng$DATE, format="%d.%m.%Y")

housng_qrt <- housng %>%
  group_by(quarter=paste(quarters(DATE), lubridate::year(DATE))) %>%
  summarise(housng = sum(HOUSNG))

housng_qrt$DATE <- as.Date(as.yearqtr(housng_qrt$quarter, format = "Q%q %Y"))
housng_qrt <- housng_qrt %>%
  select(-quarter) %>%
  subset(DATE >= start & DATE <= end)


# Import VIX anx convert to quarterly
vixcurrent <- read_csv("data/vixcurrent.csv") #Current VIX methodology
vixarchive <- read_csv("data/vixarchive.csv") #Current applied retrospectively

vixcurrent$Date <- as.Date(vixcurrent$Date, format="%m/%d/%Y")
vixarchive$Date <- as.Date(vixarchive$Date, format="%d.%m.%Y")

vix <-rbind(vixarchive, vixcurrent)
vix <- vix %>%
  select(VIX='VIX Close', 'Date') %>%
  mutate(month=month(Date)) %>%
  mutate(year=year(Date))

vix_monthly <- aggregate(VIX ~ month + year, vix, mean)

vix_qrt<- vix_monthly  %>%
  mutate(DATE=as.Date(paste(year,month,"01", sep="/"))) %>%
  select(-month, -year)  %>%
  group_by(quarter = paste(quarters(DATE), lubridate::year(DATE))) %>%
  summarise(VIX = mean(VIX)) 


vix_qrt$DATE <- as.Date(as.yearqtr(vix_qrt$quarter, format = "Q%q %Y"))
vix_qrt <- vix_qrt %>%
  select(-quarter) %>%
  subset(DATE >= start & DATE <= end)

# Import credit rates, add frac=1 for month end. Data ends in 2010??
moodys.raw <- read_csv("Data/moodys.csv")
creditspreads<- subset(moodys.raw, DATE >= start & DATE <= end)
creditspreads$AAA <- as.numeric(creditspreads$AAA)
creditspreads$DBAA <- as.numeric(creditspreads$DBAA)
creditspreads <- creditspreads %>%
  mutate(baa_aaa=DBAA-AAA)


GZ.raw <- read_csv("Data/GZ_quarterly.csv")
GZ.raw$Date <- as.Date(as.yearqtr(GZ.raw$date, format = "%YQ%q"))
GZ <- subset(GZ.raw, Date >= start & Date <=end)
GZ <- GZ %>%
  select(Date, gz_spr)

# Import term rates
DGS.clean <- read_csv("Data/feds/DGS_combined.csv")

DGS.clean$DGS3MO <- as.numeric(DGS.clean$DGS3MO)
DGS.clean$DGS6MO <- as.numeric(DGS.clean$DGS6MO)

DGS.clean <- DGS.clean %>%
  mutate(TRM0503 = DGS5-DGS3MO) %>%
  mutate(TRM0506 = DGS5-DGS6MO) %>%
  mutate(TRM1003 = DGS10-DGS3MO) %>%
  mutate(TRM1006 = DGS10-DGS6MO) %>%
  mutate(TRM1012 = DGS10-DGS1) %>%
  mutate(SRT03M = DGS3MO-lag(DGS3MO, 4))

DGS.clean<- subset(DGS.clean, DATE >= start & DATE <= end)

#SPY
spy <- read_csv("data/SPY.csv")
spy <- spy %>%
  select(Date, Close) %>%
  group_by(quarter=paste(quarters(Date), lubridate::year(Date))) %>%
  summarise_at(c("Close"), sum)

spy$Date <- as.Date(as.yearqtr(spy$quarter, format = "Q%q %Y"))
spy <- spy%>%
  arrange(Date) %>%
  mutate(spy_logreturn=100*(log(Close)- log(lag(Close, 4)))) %>%
  select(-quarter, -Close) %>%
  subset(Date >= start & Date <= end)

#Import NBER recessions and create subsets
dummy <- read.csv("data/Dummy.csv")

# Make dataset with quarterly gdp
gdp_qoq <- read_csv("data/GDPC1_qoq.csv")
gdp_qoq <- gdp_qoq %>%
  subset(DATE >= start & DATE <= end)

dataset_quarterly <- dataset
dataset_quarterly$GDP <- unlist(lapply(gdp_qoq$GDPC1_PCH, as.numeric))

# combine
dataset_quarterly <- dataset_quarterly %>%
  mutate(DGS.clean[-1]) %>% 
  mutate(creditspreads[-1]) %>%
  mutate(dum=dummy$Reccession) %>%
  mutate(vix_qrt[-2]) %>%
  mutate(housng_qrt[-2]) %>%
  mutate(log_gdp=log(1+GDP/100)*100) %>%
  left_join(GZ) %>% 
  left_join(spy) %>%
  relocate(Date, YIV, GDP)

dataset <- dataset %>%
  mutate(DGS.clean[-1]) %>% 
  mutate(creditspreads[-1]) %>%
  mutate(dum=dummy$Reccession) %>%
  mutate(vix_qrt[-2]) %>%
  mutate(housng_qrt[-2]) %>%
  mutate(log_gdp=log(1+GDP/100)*100) %>%
  left_join(GZ) %>% 
  left_join(spy) %>%
  relocate(Date, YIV, GDP)


#For summary statistics
df_summary <- dataset %>%
  select("YIV", "GDP", "spy_logreturn","VIX","DBAA","AAA","baa_aaa","gz_spr", "housng" , "SRT03M", "TRM1003", "TRM1006", "TRM1012", "TRM0503", "TRM0506")


#standardize only independent + spy(since log) <-- SISESTA SIIA NEED MIDA POLE VAJA STANDARDIZEDA
standardizable_var <- setdiff(ls(dataset), c("GDP","spy_logreturn", "dum", "Date","log_gdp"))

dataset <- dataset %>%
  mutate_at(standardizable_var, scale)

dataset_quarterly <- dataset_quarterly %>%
  mutate_at(standardizable_var, scale)


#Create rolling averages
dataset <- dataset %>%
  arrange %>%
  mutate(H1=rollapply(log_gdp,2,FUN = function(df) mean(df[-2], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H2=rollapply(log_gdp,3,FUN = function(df) mean(df[-3], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H4=rollapply(log_gdp,5,FUN = function(df) mean(df[-5], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H6=rollapply(log_gdp, 7,FUN = function(df) mean(df[-7], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H8=rollapply(log_gdp, 9,FUN = function(df) mean(df[-9], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H10=rollapply(log_gdp, 11,FUN = function(df) mean(df[-11], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H12=rollapply(log_gdp, 13,FUN = function(df) mean(df[-13], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(F1=lead(log_gdp, n = 1L)) %>%
  mutate(F2=lead(log_gdp, n = 2L)) %>%
  mutate(F4=lead(log_gdp, n = 4L)) %>%
  mutate(F8=lead(log_gdp, n = 8L))

dataset_qoq <- dataset_quarterly %>%
  arrange %>%
  mutate(H1=rollapply(log_gdp,2,FUN = function(df) mean(df[-2], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H2=rollapply(log_gdp,3,FUN = function(df) mean(df[-3], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H4=rollapply(log_gdp,5,FUN = function(df) mean(df[-5], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H6=rollapply(log_gdp, 7,FUN = function(df) mean(df[-7], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H8=rollapply(log_gdp, 9,FUN = function(df) mean(df[-9], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H10=rollapply(log_gdp, 11,FUN = function(df) mean(df[-11], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(H12=rollapply(log_gdp, 13,FUN = function(df) mean(df[-13], na.rm = TRUE), fill = NA, align = "left" )) %>%
  mutate(F1=lead(log_gdp, n = 1L)) %>%
  mutate(F2=lead(log_gdp, n = 2L)) %>%
  mutate(F4=lead(log_gdp, n = 4L)) %>%
  mutate(F8=lead(log_gdp, n = 8L))


#Name differences with appendices
df <- dataset 
df.expansionary <- subset(dataset, dum== 0)
df.recessionary <- subset(dataset, dum== 1)
df_qoq <- dataset_qoq

#drop unnecessary shit for work proccesses
rm(list=setdiff(ls(), c("df", "df_summary", "df_qoq")))


statistics <- df_summary %>%
  descr(
    transpose = TRUE,
    stats = c("mean","sd","min","q1","med","q3","max","n.valid"))

statistics <- statistics %>%
  mutate(Variable=rownames(statistics)) %>%
  relocate(Variable)


