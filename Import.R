#Function for loading packages and installing them in case they are not already installed.
install_packages <- function(package){
  new.package <- package[!(package %in% installed.packages()[, "Package"])]
  if (length(new.package)) 
    install.packages(new.package, dependencies = TRUE)
  sapply(package, require, character.only = TRUE)
}
packages <- c("tidyverse", "summarytools", "kableExtra", "dynlm", "stats", "sandwich", "data.table", "zoo", "caret", "randomForest", "broom", "modelr", "lmtest", "gridExtra")
install_packages(packages)



start <- as.Date("01/03/1990", format="%d/%m/%Y")
end <- as.Date("01/10/2015", format="%d/%m/%Y")

#Import extracted GDP and YIV
dataset <- read_csv("data/gdp_yiv.csv") %>%
  rename(dum=Reccession)

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

# Import credit rates
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
DGS.clean <- read_csv("Data/DGS_combined.csv")

DGS.clean$DGS3MO <- as.numeric(DGS.clean$DGS3MO)
DGS.clean$DGS6MO <- as.numeric(DGS.clean$DGS6MO)

DGS.clean <- DGS.clean %>%
  mutate(TRM0503 = DGS5-DGS3MO) %>%
  mutate(TRM0506 = DGS5-DGS6MO) %>%
  mutate(TRM1003 = DGS10-DGS3MO) %>%
  mutate(TRM1006 = DGS10-DGS6MO) %>%
  mutate(TRM1012 = DGS10-DGS1) %>%
  mutate(SRT03M = DGS3MO-lag(DGS3MO, 1))

DGS.clean<- subset(DGS.clean, DATE >= start & DATE <= end)

#Import SPY
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


#Combine imported variables into a dataset
dataset <- dataset %>%
  mutate(DGS.clean[-1]) %>% 
  mutate(creditspreads[-1]) %>%
  mutate(vix_qrt[-2]) %>%
  mutate(housng_qrt[-2]) %>%
  mutate(log_gdp=log(1+GDP/100)*100) %>%
  left_join(GZ) %>% 
  left_join(spy) %>%
  relocate(Date, YIV, GDP)


#For summary statistics
df_summary <- dataset %>%
  select("YIV", "GDP", "VIX", "DBAA","AAA","baa_aaa", "housng" , "SRT03M", "TRM1003", "TRM1006", "TRM1012", "TRM0503", "TRM0506", "DGS3MO")


#Excluding variables in " " from standardization
standardizable_var <- setdiff(ls(dataset), c("GDP","spy_logreturn", "dum", "Date","log_gdp"))

dataset <- dataset %>%
  mutate_at(standardizable_var, scale)


#Defining expansionary and recessionary variables
df.expansionary <- subset(dataset, dum== 0)
df.recessionary <- subset(dataset, dum== 1)

#Importing different variations of dependent variable used in academia from excel
cleandata <- read_csv("data/RAWDATA.csv") %>%
  select(-1, -2)
df <- cbind(dataset, cleandata)


#Drop unnecessary data frames created in the process and leaving only the following
rm(list=setdiff(ls(), c("df", "df_summary", "df_qoq")))

# For summary statistics table
statistics <- df_summary %>%
  descr(
    transpose = TRUE,
    stats = c("mean","sd","min","q1","med","q3","max","n.valid"))

statistics <- statistics %>%
  mutate(Variable=rownames(statistics)) %>%
  relocate(Variable)

summary_table <- colnames(df_summary)
summary_table <- summary_table %>%
  as_tibble() %>%
  rename(Variable=value)
summary_table$Description <- c("5 - year Treasury Implied Volatility",
                               "Real gross domestic product",
                               "Returns of VIX index",
                               "BAA corporate bond yields",
                               "AAA corporate bond yields",
                               "Yield spread between BAA and AAA yields",
                               "New housing market starts",
                               "Changes in 3 month treasury yield",
                               "TRM1003 - 10 year and 3 month treasury yield spread",
                               "TRM1006 - 10 year and 6 month treasury yield spread",
                               "TRM1012 - 10 year and 1 year treasury yield spread",
                               "TRM0503 - 5 year and 3 month treasury yield spread",
                               "TRM0506 - 5 year and 6 month treasury yield spread",
                               "Three month corporate bond yield")


