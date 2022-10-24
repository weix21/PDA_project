## Xin and Chichun
## Porcess FI data for further analysis
## Date: 10/23/2022


# Load necessary packages 
library(tidyverse)
library(forcats)
library(ggpubr)
library(hydroTSM)
library(zoo)
library(stringr)


# Set working directory
setwd("/Users/weixin/Downloads/study/Brown/PDA/Project/data")


## Load data from New England
data_dir_1 <-  c("Salmonella enterica_New_England.csv","E.coli and Shigella_New_England.csv","Campylobacter jejuni_New_England.csv","Listeria monocytogenes_New_England.csv")
data_list_1  <-  lapply(data_dir_1, read.csv)
names(data_list_1) <- strsplit(data_dir_1,split=".csv")

data_1 <- do.call(rbind.data.frame, data_list_1)

## Format the State 
data_1$State <- fct_collapse(
  data_1$Location,
  'CT'= c("USA:CT", "USA: CT", "USA: Connecticut, Greenwich", "USA: Connecticut"),
  'MA'= c("USA:MA", "USA: MA", "USA: Massachusetts", "USA: Cambridge, MA, MIT", "USA:Boston", "USA: Massachusetts, Boston", 
          "USA: Boston", "USA: Cambridge, MA", "USA: Boston, MA", "USA: MA, Boston", "USA:Massachusetts, Beverly", "USA: Ipswich, MA"),
  'ME'= c("USA:ME", "USA: ME", "USA: Maine", "USA:Maine"),
  'NH'= c("USA:NH", "USA: New Hampshire", "USA: New Hampshire,Rindge", "USA: NH"),
  'RI'= c("USA:RI", "USA: Rhode Island", "USA: RI"),
  'VT'= c("USA:VT", "USA: VT")) 

data_1$Region <- "New England"

## Load data from CA
data_dir_2 <-  c("Salmonella enterica_CA.csv","E.coli and Shigella_CA.csv","Campylobacter jejuni_CA.csv","Listeria monocytogenes_CA.csv")
data_list_2  <-  lapply(data_dir_2, read.csv)
names(data_list_2) <- strsplit(data_dir_2,split=".csv")

data_2 <- do.call(rbind.data.frame, data_list_2)
data_2$State <- "CA"
data_2$Region <- "CA"

## Load data from TX
data_dir_3 <-  c("Salmonella enterica_TX.csv","E.coli and Shigella_TX.csv","Campylobacter jejuni_TX.csv","Listeria monocytogenes_TX.csv")
data_list_3  <-  lapply(data_dir_3, read.csv)
names(data_list_3) <- strsplit(data_dir_3,split=".csv")

data_3 <- do.call(rbind.data.frame, data_list_3)
data_3$State <- "TX"
data_3$Region <- "TX"

## Combine data and re-factorize
df <- rbind(data_1,data_2,data_3)


df$State <- factor(df$State, levels =c(sort(as.character(unique(data_1$State))),"CA","TX"))
df$Region <- factor(df$Region, levels =c("New England","CA","TX"))
df$X.Organism.group <- factor(df$X.Organism.group, 
                              levels = c("Salmonella enterica", "E.coli and Shigella", 
                                         "Campylobacter jejuni","Listeria monocytogenes"))


## Create time related variables
df$Date <- as.Date(substr(df$Create.date,1,10),format="%Y-%m-%d")
df$Season <- time2season(df$Date, out.fmt = "seasons")
df$Quarter <- as.yearqtr(df$Date,format = "%Y-%m-%d")

## Notice that there's typo in results from time2season()
df$Season[df$Season == "autumm"] <- "autumn"
df$Season <- str_to_title(df$Season)

## Create count variables
df <- df %>% add_count(Date,X.Organism.group,Region, name="nDate")
df <- df %>% add_count(Season,X.Organism.group,Region, name="nSeason")
df <- df %>% add_count(Quarter,X.Organism.group,Region, name="nQuarter")

df <- df %>% add_count(Isolation.source,X.Organism.group,Region, name="nIsolation.source")
df <- df %>% add_count(Serovar,X.Organism.group,Region, name="nSerovar")

## Format missing data as NA
df[df==""] <- NA
df$Serovar[df$Serovar=="not typed"] <- NA

## Save data
save(df, file = "Processed_data_10_23.RData")

