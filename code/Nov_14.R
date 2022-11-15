# load package
library(tidyverse)
library(ggplot2)
library(forcats)
library(ggpubr)
library(hydroTSM)
library(zoo)
library(stringr)
library(lubridate)

# set working direction
setwd("~/OneDrive - Brown University/PHP 2550/Final Project/Data")
load("Processed_data_11_11.RData")

# Only use the data collected after 2000 with Month record
df_2000 <- df[!is.na(as.numeric(df$Collection_Month)) & as.numeric(df$Collection_Year) >= 2000,] 

# Factorize string variables
df_2000$Collection_Month_Factor <- factor(df_2000$Collection_Month,levels = c(1:12))
df_2000$Collection_Season <- factor(df_2000$Collection_Season,levels = c("Spring","Summer","Autumn","Winter"))
df_2000$X.Organism.group <- factor(df_2000$X.Organism.group,levels = c("Salmonella enterica","E.coli and Shigella","Listeria monocytogenes","Campylobacter jejuni"))

df_2000$Isolation.source.cdc <- df_2000$Isolation.source
df_2000$Isolation.source.cdc[df_2000$Isolation.type == "clinical"] <- "clinical"
df_2000$Isolation.source.cdc[grep("chicken|Chicken|poultry|turkey|Turkey|POULTRY|duck|Duck",df_2000$Isolation.source)] <- "poultry"
df_2000$Isolation.source.cdc[grep("salmon|cod|tuna|surimi|fish|caviar|catfish|fillet",df_2000$Isolation.source)] <- "fish"
df_2000$Isolation.source.cdc[grep("crab|lobster",df_2000$Isolation.source)] <- "crustacens"
df_2000$Isolation.source.cdc[grep("oyster|sponge|clam",df_2000$Isolation.source)] <- "mollusks"
df_2000$Isolation.source.cdc[grep("Milk|milk|cheese|butter|cream|yogurt|cow",df_2000$Isolation.source)] <- "diary"
df_2000$Isolation.source.cdc[grep("Egg|egg",df_2000$Isolation.source)] <- "eggs"
df_2000$Isolation.source.cdc[grep("beef|cattle|Beef|bovine",df_2000$Isolation.source)] <- "beef"
df_2000$Isolation.source.cdc[grep("plover|tern|cormorant|iguana|gull|siskin|mouse|marmoset|macaque|Hedgedog|rat|Mouse|heron|snake|bird|rabbit|deer|kangaroo|bango|zebra|lizard",df_2000$Isolation.source)] <- "wild"
df_2000$Isolation.source.cdc[grep("pig|pork|Pork|Pig|swine|prosciutto|salami|bologna|porcine|Porcine|ham",df_2000$Isolation.source)] <- "pork"
df_2000$Isolation.source.cdc[grep("lamb|Lamb|goat",df_2000$Isolation.source)] <- "mutton"
df_2000$Isolation.source.cdc[grep("soy|Soy|bean|edamame",df_2000$Isolation.source)] <- "beans"
df_2000$Isolation.source.cdc[grep("flour|Noodles|noodle|canola|bran|cereal|cornmeal|wheat",df_2000$Isolation.source)] <- "oils-sugar-wheat"
df_2000$Isolation.source.cdc[grep("hazelnut|Chia|Coconut|coconut|cantaloupe|nut|peanut|avocado|almond|pistachio|mango|peach|seed|cashew|cherry|blueberry|kiwi|nectarine|pecan|apple",df_2000$Isolation.source)] <- "fruit-nut"
df_2000$Isolation.source.cdc[grep("vegetable|pepper|lettuce|romaine|Romaine|onion|broccoli|salad|Beans|sprout|Squash|
                                  cardamo|spinach|tarragon|cilantro|kale|arugula|parsley|celery|root|leaf|avocado|carrot|tomato|mashroom|squash|chard|
                                  mushroom",df_2000$Isolation.source)] <- "vegetables"
df_2000$Isolation.source.cdc[grep("Cat|cat|Dog|dog|Canine|canine|Feline|feline",df_2000$Isolation.source)] <- "dog-cat"
df_2000$Isolation.source.cdc[grep("meat|Meat|sausage|Sausage",df_2000$Isolation.source.cdc)] <- "other-meat"
df_2000$Isolation.source.cdc[grep("Urine|urine|swab|Swab|Blood|blood",df_2000$Isolation.source.cdc)] <- "unknown-sample-source"
df_2000$Isolation.source.cdc[!df_2000$Isolation.source.cdc %in% c("poultry","fish","crustacens","diary","eggs","mollusks","beef","wild","pork","mutton","beans","oils-sugar-wheat","fruit-nut","vegetables","dog-cat","unknown-sample-source","clinical")] <- "others"

## Explore the realtionship between region and number of cases on species
df_2000 %>% select(Region, X.Organism.group) %>% table()

df_2000 %>% ggplot(aes(X.Organism.group)) +
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme_minimal() +
  coord_flip() +
  facet_wrap(~as.factor(Region), nrow=3)

## Explore the relationship between region, season and number of cases on species 
# Count of total cases for season and region
df_2000 %>% ggplot(.,aes(x=Collection_Season)) + 
  geom_bar(stat="count", width=0.7) +
  labs(x = "Season",y = "Count") + 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(~Region)

# Count of pathogens for season and region
df_2000 %>% ggplot(.,aes(x=Collection_Season,fill = X.Organism.group)) + 
  geom_bar(stat="count", width=0.7)+
  labs(x = "Season",y = "Count") + 
  guides(fill=guide_legend(title="Organism")) + 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(Region ~ X.Organism.group)

# df_2000  %>% ggplot(.,aes(x=factor(Collection_Season, levels = c("Spring","Summer","Autumn","Winter")))) + 
#   geom_bar(stat="count", width=0.7)+
#   xlab("Season")+ 
#   theme(axis.text.x = element_text(size=8, angle = 30)) +
#   facet_grid(as.factor(Region) ~ as.factor(X.Organism.group))

# Propoption of pathogens for season and region
df_2000 %>% 
  ggplot(.,aes(x=Collection_Season, fill = X.Organism.group)) + 
  geom_bar(position = "fill", width=0.7) +
  labs(x = "Season",y = "Proportion") + 
  guides(fill=guide_legend(title="Organism")) + 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(~Region)

## Explore the relationship between region, month and number of cases on species
# Count of total cases for month and region
df_2000 %>% ggplot(.,aes(x=Collection_Month_Factor)) + 
  geom_bar(stat="count", width=0.7) +
  labs(x = "Month",y = "Count") + 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(~Region)

# Count of pathogens for season and region
df_2000 %>% ggplot(.,aes(x=Collection_Month_Factor,fill = X.Organism.group)) + 
  geom_bar(stat="count", width=0.7)+
  labs(x = "Month",y = "Count") + 
  guides(fill=guide_legend(title="Organism")) + 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(Region ~ X.Organism.group)

# Propoption of pathogens for season and region
df_2000 %>% 
  ggplot(.,aes(x=Collection_Month_Factor, fill = X.Organism.group)) + 
  geom_bar(position = "fill", width=0.7) +
  labs(x = "Month",y = "Proportion") + 
  guides(fill=guide_legend(title="Organism")) + 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(~Region)


# Top sources for each Organism in each region
df_2000 <- df_2000 %>% add_count(Isolation.source.cdc,X.Organism.group,Region, name="nIsolation.organ.source.cdc")
df_2000 <- df_2000 %>% add_count(Isolation.source.cdc,Collection_Season,Region, name="nIsolation.season.source.cdc")
Season_type <- c("Spring","Summer","Autumn","Winter")
Region_type <- c("New England","CA","TX")
Organism_group <- c("Salmonella enterica","E.coli and Shigella","Listeria monocytogenes","Campylobacter jejuni")
p <- NULL
for(i in 1:3){
  for(j in 1:4){
    p[[4*(i-1)+j]] <- df_2000 %>% 
      filter(!Isolation.source.cdc %in% c("others","unknown-sample-source"),Collection_Season == Season_type[j], Region == Region_type[i]) %>%    
      mutate(rank = dense_rank(-nIsolation.season.source.cdc)) %>% arrange(desc(nIsolation.season.source.cdc)) %>% filter(rank<= 5) %>%
      ggplot(.,aes(x=factor(Isolation.source.cdc,levels=rev(unique(Isolation.source.cdc))))) + geom_bar(stat="count", width=0.7,fill = "red3") + labs(x = Region_type[i], y = Season_type[j]) + ylim(c(0,450)) +
      theme(axis.text.x = element_text(angle = -45),legend.position = NULL)
  }
}

commonplot1 <- ggarrange(plotlist=p, nrow = 3, ncol=4)
commonplot1



# time series counting for pathogens
df_2000 <- df_2000 %>% add_count(Collection_YM,X.Organism.group,Region, name="nCollection_YM")
df_2000 %>% ggplot(.,aes(x=Collection_YM,y=nCollection_YM)) + 
  geom_line(colour = 'blue') + 
  theme(axis.text.x = element_text(size=8, angle = 45)) +
  facet_grid(as.factor(Region) ~ as.factor(X.Organism.group))



df_2000 <- df_2000 %>% add_count(Collection_Quarter,X.Organism.group,Region, name="nCollection_Quarter")
df_2000 %>% ggplot(.,aes(x=Collection_Quarter,y=nCollection_Quarter)) + 
  geom_line(colour = 'blue')  + 
  theme(axis.text.x = element_text(size=7, angle = 30)) + 
  facet_grid(as.factor(Region) ~ as.factor(X.Organism.group))


# AMR genotype 
genotypes <- strsplit(df_2000$AMR.genotypes,split=",")
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
for(i in 1:nrow(df_2000)){
  genotypes_Season <- c(genotypes_Season,rep(df_2000$Collection_Season[i],length(clean_genotypes[[i]])))
  genotypes_Month <- c(genotypes_Month,rep(df_2000$Collection_Month[i],length(clean_genotypes[[i]])))
  genotypes_Organ <- c(genotypes_Organ,rep(as.character(df_2000$X.Organism.group[i]),length(clean_genotypes[[i]])))
  genotypes_Strain <- c(genotypes_Strain,rep(as.character(df_2000$Strain[i]),length(clean_genotypes[[i]])))
  genotypes_Region <- c(genotypes_Region,rep(as.character(df_2000$Region[i]),length(clean_genotypes[[i]])))
  
}
AMR_df <- data.frame(gene = unlist(clean_genotypes),
                     type = unlist(clean_genotypes_categories),
                     organism = genotypes_Organ,
                     strain = genotypes_Strain,
                     season = genotypes_Season,
                     month = genotypes_Month,
                     region = genotypes_Region
                     )

# Top AMR for each Organism in each region
AMR_df$month <- factor(AMR_df$month,levels = c(1:12))
AMR_df$season <- factor(AMR_df$season,levels = c(1:4),labels = c("Spring","Summer","Autumn","Winter"))
AMR_df$region <- factor(AMR_df$region)
AMR_df <- AMR_df %>% add_count(gene,season,region, name="ngene_season") %>%
  add_count(gene,month,region, name="ngene_month")

Season_type <- c("Spring","Summer","Autumn","Winter")
Region_type <- c("New England","CA","TX")
p <- NULL
for(i in 1:3){
  for(j in 1:4){
    p[[4*(i-1)+j]] <- AMR_df %>% 
      filter(season == Season_type[j], region == Region_type[i]) %>%    
      mutate(rank = dense_rank(-ngene_season)) %>% arrange(desc(ngene_season)) %>% filter(rank<= 5) %>%
      ggplot(.,aes(x=factor(gene,levels=rev(unique(gene))))) + geom_bar(stat="count", width=0.7,fill = "mediumblue") + labs(x = Region_type[i], y = Season_type[j]) + ylim(c(0,850)) +
      theme(axis.text.x = element_text(angle = -45),legend.position = NULL)
  }
}

commonplot2 <- ggarrange(plotlist=p, nrow = 3, ncol=4)
commonplot2

# Strain
strain <- df_2000$Strain
length(grep("20CT",strain))
length(grep("21CT",strain))
length(grep("22CT",strain))
length(grep("BIOML",strain))
length(grep("CFSA",strain))/nrow(df_2000)
length(grep("PSU_",strain))
length(grep("NY_",strain))
length(grep("FDA",strain))/nrow(df_2000)
length(grep("PNUS",strain))
length(grep("HEC",strain))
length(grep("BEC",strain))
length(grep("CVM",strain))
length(grep("RO_",strain))
length(grep("MJ",strain))
length(grep("CDPH",strain))
length(grep("UMDU",strain))
length(grep("RM",strain))/nrow(df_2000)

png(file = "/Users/chichuntan/OneDrive - Brown University/PHP 2550/Final Project/Sources_AMR.png",width = 1000,height = 800)
cowplot::plot_grid(commonplot1,commonplot2,nrow = 2,labels = c("Isolation Sources","Detected AMR genotype"))
dev.off()
