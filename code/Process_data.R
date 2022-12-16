# Load necessary packages 
library(tidyverse)
library(forcats)
library(ggpubr)
library(hydroTSM)
library(zoo)
library(stringr)
library(lubridate)


# Set working directory
setwd("PDA/Project/data")


## Load data from New England
data_dir_1 <-  c("Salmonella enterica_USA_1.csv","Salmonella enterica_USA_2.csv","E.coli and Shigella_USA.csv","Campylobacter jejuni_USA.csv","Listeria monocytogenes_USA.csv")
data_list_1  <-  lapply(data_dir_1, read.csv)
names(data_list_1) <- strsplit(data_dir_1,split=".csv")

data_1 <- do.call(rbind.data.frame, data_list_1)


location = strsplit(data_1$Location,split=":")
location = do.call(rbind, lapply(location, function(x) rev(`length<-`(x, 3))))

data_1 <- data_1[!is.na(location[,2])&location[,3]=="USA", ]

## Format the State 
data_1$State <- fct_collapse(
  data_1$Location,
  'CT'= c("USA:CT", "USA: CT", "USA: Connecticut, Greenwich", "USA: Connecticut"),
  'MA'= c("USA:MA", "USA: MA", "USA: Massachusetts", "USA: Cambridge, MA, MIT", "USA:Boston", "USA: Massachusetts, Boston", 
          "USA: Boston", "USA: Cambridge, MA", "USA: Boston, MA", "USA: MA, Boston", "USA:Massachusetts, Beverly", "USA: Ipswich, MA"),
  'ME'= c("USA:ME", "USA: ME", "USA: Maine", "USA:Maine"),
  'NH'= c("USA:NH", "USA: New Hampshire", "USA: New Hampshire,Rindge", "USA: NH"),
  'RI'= c("USA:RI", "USA: Rhode Island", "USA: RI"),
  'VT'= c("USA:VT", "USA: VT"),
  'CA'= c("USA: CA", "USA: California", "USA: California, Alameda County", "USA: Califronia", 
          "USA: Irvine, CA","USA: Merced, CA","USA: Palo Alto, California","USA:Sacramento CA",
          "USA:Ronald Reagan UCLA Medical Center","USA:Ocean Beach, California",
          "USA:CA-UC Davis","USA:CA-Contra Costa","USA:CA-Contra Costa","USA:CA",
          "USA: UCSD Medical Center, San Diego, California","USA: San Diego, CA",
          "USA: San Diego","USA: Merced, CA","USA: La Jolla, California",
          "USA: San Francisco", "USA:San Francisco"),
  'TX'= c("USA: San Antonio","USA: Dallas, Texas","USA: Houston","USA: Houston, TX","USA: Texas",
          "USA: Texas Panhandle", "USA: TX","USA:Dallas TX","USA:Houston","USA:Texas, Amarillo",
          "USA:Texas, Cactus","USA:Texas, Falls City", "USA:Texas, Harlingen", "USA:Texas, Herefore", 
          "USA:Texas,Corpus Christi", "USA:TX"),
  other_level = "others") 

data_1$Region <- fct_collapse(
  data_1$State,
  "New England" = c("MA", "CT", "ME", "NH", "RI", "VT"),
  other_level = "others") 


df <- data_1


## Factorize State, Region and Species 
df$State <- factor(df$State, levels =c(sort(c("MA", "CT", "ME", "NH", "RI", "VT")),"others"))
df$Region <- factor(df$Region, levels =c("New England","others"))
df$X.Organism.group <- factor(df$X.Organism.group, 
                              levels = c("Salmonella enterica", "E.coli and Shigella", 
                                         "Campylobacter jejuni","Listeria monocytogenes"))


## Remove abnormal points
df <- df[-which(df$Collection_YM=="2022-11-01"),]


## Adjust for population
data_2 <- read.csv("us_pop.csv")
data_2$Collection_YM <- data_2$DATE
data_2$DATE <- substr(data_2$DATE,1,4)
data_2 <- data_2[data_2$DATE>=2010,]
data_2 <- data_2 %>% select(-c('DATE'))
names(data_2) <- c("pop","Collection_YM")
data_2$prop <- 1/ data_2$pop*1000

df = merge(df,data_2,by=Collection_YM,all.x = T)



## Create time related variables
df$Create_date <- as.Date(substr(df$Create.date,1,10),format="%Y-%m-%d")
df$Create_YM <- as.Date(paste(substr(df$Create.date,1,7),01,sep="-"),format="%Y-%m-%d")
df$Create_Month <- month(as.POSIXlt(df$Create_date, format = "%m-%d-%y"))
df$Create_Season <- time2season(df$Create_date, out.fmt = "seasons")
df$Create_Quarter <- as.yearqtr(df$Create_date,format = "%Y-%m-%d")
df$Create_Year <- substr(df$Create.date,1,4)

df$Collection_Date <- as.Date(substr(df$Collection.date,1,10),format="%Y-%m-%d")
df$Collection_YM <- as.Date(paste(substr(df$Collection.date,1,7),01,sep="-"),format="%Y-%m-%d")
df$Collection_Month <- month(as.POSIXlt(df$Collection_YM, format = "%m-%d-%y"))
df$Collection_Season <- time2season(df$Collection_YM, out.fmt = "seasons")
df$Collection_Quarter <- as.yearqtr(df$Collection_YM,format = "%Y-%m-%d")
df$Collection_Year <- substr(df$Collection.date,1,4)



## Notice that there's typo in results from time2season()
df$Create_Season[df$Create_Season == "autumm"] <- "autumn"
df$Create_Season <- str_to_title(df$Create_Season)

df$Collection_Season[df$Collection_Season == "autumm"] <- "autumn"
df$Collection_Season <- str_to_title(df$Collection_Season)



## Create count variables
df <- df %>% add_count(Collection_Date,X.Organism.group,Region, name="nDate")
df <- df %>% add_count(Collection_Season,X.Organism.group,Region, name="nSeason")
df <- df %>% add_count(Collection_Quarter,X.Organism.group,Region, name="nQuarter")


## Format missing data as NA
df[df==""] <- NA
df$Serovar[df$Serovar=="not typed"] <- NA



# factor: Month, Quarter, Pathogen
df$Collection_Month_Factor <- factor(df$Collection_Month,levels = c(1:12))
df$Collection_Season <- factor(df$Collection_Season,levels = c("Spring","Summer","Autumn","Winter"))
df$X.Organism.group <- factor(df$X.Organism.group,levels = c("Salmonella enterica","E.coli and Shigella","Listeria monocytogenes","Campylobacter jejuni"))




# recategory isolation sources
df$Isolation.source.cdc <- df$Isolation.source
df$Isolation.source.cdc[df$Isolation.type == "clinical"] <- "clinical"
df$Isolation.source.cdc[grep("chicken|Chicken|CHICKEN|poultry|turkey|TURKEY|Turkey|POULTRY|duck|Duck|WINGS|Wings|Breast|
                             BREASTS|Drums|Thighs|THIGHS|Poultry|Gallus|hen|Meleagris|goose|meleagris",df$Isolation.source)] <- "poultry"
df$Isolation.source.cdc[grep("salmon|cod|tuna|surimi|fish|caviar|catfish|fillet|shrimp|Fish|fish",df$Isolation.source)] <- "fish"
df$Isolation.source.cdc[grep("crab|lobster",df$Isolation.source)] <- "crustacens"
df$Isolation.source.cdc[grep("oyster|sponge|clam",df$Isolation.source)] <- "mollusks"
df$Isolation.source.cdc[grep("Milk|milk|cheese|butter|cream|yogurt|cow|Cheese",df$Isolation.source)] <- "diary"
df$Isolation.source.cdc[grep("Egg|egg",df$Isolation.source)] <- "eggs"
df$Isolation.source.cdc[grep("beef|cattle|Beef|bovine|BEEF|Bovine|BOVINE|bov|Bos taurus|bos taurus",df$Isolation.source)] <- "beef"
df$Isolation.source.cdc[grep("plover|tern|cormorant|iguana|gull|siskin|mouse|marmoset|macaque|Hedgedog|rat|Mouse|heron|snake|bird|rabbit|deer|kangaroo|
                             bango|zebra|lizard|turtle|monkey|Corvus|Odocoileus|Quiscalus|Phasianus|Loxia|monkey|raccoon|Bison|bison|Cardinalis|hedgehog|Alces alces|Anas|Antelope|Aythya|Canis|canis|gecko|giraffe",df$Isolation.source)] <- "wild"
df$Isolation.source.cdc[grep("pig|pork|Pork|PORK|Pig|swine|prosciutto|salami|bologna|porcine|Porcine|ham|Sus scrofa|Sus|Swine",df$Isolation.source)] <- "pork"
df$Isolation.source.cdc[grep("lamb|Lamb|goat|Ovis|Ovine|ovine|Capra|Caprinae|Equus|Horse|horse|equus|equine|Equine",df$Isolation.source)] <- "mutton-horse"
df$Isolation.source.cdc[grep("soy|Soy|bean|edamame",df$Isolation.source)] <- "beans"
df$Isolation.source.cdc[grep("flour|Noodles|noodle|canola|bran|cereal|cornmeal|wheat|rice|Cereal|lasagna|biscuit|dough",df$Isolation.source)] <- "oils-sugar-wheat"
df$Isolation.source.cdc[grep("hazelnut|Chia|Coconut|coconut|cantaloupe|nut|peanut|avocado|almond|pistachio|mango|peach|seed|cashew|cherry|blueberry|kiwi|
                             nectarine|pecan|apple|fruit|grape|lemon|papaya|Apple|watermelon|strawberry|melon",df$Isolation.source)] <- "fruit-nut"
df$Isolation.source.cdc[grep("vegetable|pepper|lettuce|romaine|Romaine|onion|broccoli|salad|Beans|sprout|Spinach|Squash|Avocado|
                                  cardamo|spinach|tarragon|cilantro|kale|arugula|parsley|celery|root|leaf|avocado|carrot|tomato|mashroom|squash|chard|
                                  mushroom|cucumber|chili|potato|ginger|cinnamon|Avocado|asparagus|cabbage|mushroom|Cabbage|pea|corn",df$Isolation.source)] <- "vegetables"
df$Isolation.source.cdc[grep("Cat|cat|Dog|dog|Canine|canine|Feline|feline|pet|Pet",df$Isolation.source)] <- "pet"
df$Isolation.source.cdc[grep("meat|Meat|sausage|Sausage",df$Isolation.source.cdc)] <- "other-meat"
df$Isolation.source.cdc[grep("Urine|urine|swab|Swab|Blood|blood",df$Isolation.source.cdc)] <- "unknown-sample-source"
df$Isolation.source.cdc[!df$Isolation.source.cdc %in% c("poultry","fish","crustacens","diary","eggs","mollusks","beef","wild","pork","mutton",
                                                                  "beans","oils-sugar-wheat","fruit-nut","vegetables","dog-cat","unknown-sample-source","clinical")] <- "others"




# create new category for isolation source
df$Isolation.source.simple <- case_when(df$Isolation.source.cdc %in% c("pork","beef","mutton-horse") ~ "livestock",
                                             df$Isolation.source.cdc %in% c("fish","mollusks","crustacens") ~ "seafood",
                                             df$Isolation.source.cdc %in% c("fruit-nut","vegetables","beans") ~ "vege-fruit",
                                             df$Isolation.source.cdc %in% c("poultry","eggs") ~ "poultry",
                                             df$Isolation.source.cdc %in% c("diary") ~ "diary",
                                             df$Isolation.source.cdc %in% c("pet") ~ "pet",
                                             df$Isolation.source.cdc %in% c("oil-sugar-wheat","clinical","wild","others","unknown-sample-source") ~ "others")
table(df$Isolation.source.simple)




# AMR Phenotype
genotypes <- strsplit(df$AMR.genotypes,split=",")
clean_genotypes <- list()
clean_genotypes_categories <- list()

for(i in 1:length(genotypes)){
  if(is.na(genotypes[[i]])){
    clean_genotypes[[i]]  <- NA
    clean_genotypes_categories[[i]] <- NA
  }else{
    gene_vec <- c()
    type_vec <- c()
    for(j in 1:length(genotypes[[i]])){
      gene_vec[j] <- strsplit(genotypes[[i]][j],split="=")[[1]][1]
      type_vec[j] <- strsplit(genotypes[[i]][j],split="=")[[1]][2]
    }
    clean_genotypes[[i]]  <- gene_vec
    clean_genotypes_categories[[i]]  <- type_vec
  }
}



genotypes_Season <- c()
genotypes_Month <- c()
genotypes_Organ <- c()
genotypes_Strain <- c()
genotypes_Region <- c()
for(i in 1:nrow(df)){
  genotypes_Season <- c(genotypes_Season,rep(df$Collection_Season[i],length(clean_genotypes[[i]])))
  genotypes_Month <- c(genotypes_Month,rep(df$Collection_Month[i],length(clean_genotypes[[i]])))
  genotypes_Organ <- c(genotypes_Organ,rep(as.character(df$X.Organism.group[i]),length(clean_genotypes[[i]])))
  genotypes_Strain <- c(genotypes_Strain,rep(as.character(df$Strain[i]),length(clean_genotypes[[i]])))
  genotypes_Region <- c(genotypes_Region,rep(as.character(df$Region[i]),length(clean_genotypes[[i]])))
  
}
AMR_df <- data.frame(gene = unlist(clean_genotypes),
                     type = unlist(clean_genotypes_categories),
                     organism = genotypes_Organ,
                     strain = genotypes_Strain,
                     season = genotypes_Season,
                     month = genotypes_Month,
                     region = genotypes_Region
)






# Top 3 AMR genotype Union for all spieces
a=table(AMR_df$organism,AMR_df$gene)
top3_AMR <- c()
for(i in 1:4){
  top3_AMR <- c(top3_AMR,names(sort(a[i,],decreasing = T)[1:3]))
}
top3_AMR <- unique(top3_AMR) 

df <- cbind(df,setNames(lapply(top3_AMR, function(x) x=0), top3_AMR))
for(i in 1:nrow(df)){
  for(j in 1:length(top3_AMR)){
    if(top3_AMR[j] %in% clean_genotypes[[i]]){
      df[i,top3_AMR[j]] <- 1
    }
  }
}


# create binary variable for interested food sources
sources <- c("diary","livestock","poultry","seafood","vege-fruit") 
df <- cbind(df,setNames(lapply(sources, function(x) x=0), sources))
for(i in 1:nrow(df)){
  for(j in 1:length(sources)){
    if(sources[j] %in% df$Isolation.source.simple[i]){
      df[i,sources[j]] <- 1
    }
  }
}

# Rename variable for Species
names(df)[names(df) == 'X.Organism.group'] <- 'Species'


## Save data
save(df, file = "Processed_data.RData")

