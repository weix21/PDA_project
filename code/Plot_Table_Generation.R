# load package
library(hydroTSM)
library(zoo)
library(stringr)
library(lubridate)
library(glmnet)
source("PDA/Project/code/analysis.R")


# set working direction
setwd("PDA/Project/data")
load("Processed_data.RData")

# screen out the case without AMR genotype records
df <- df[!is.na(df$AMR.genotypes),]

# Filter to clean data 
df <- df %>% dplyr::select(-c("Strain", "Outbreak", "Isolate.identifiers", "Serovar", "Isolate", "SNP.cluster", "Host", "Host.disease",
                              "Min.same", "Min.diff", "BioSample", "Assembly", "AMR.genotypes.core", "AMRFinderPlus.version",
                              "AMRFinderPlus.analysis.type", "BioProject", "AST.phenotypes", "Computed.types", "Source.type"))


# creating the count for cases
count_df <- df %>%
  group_by(Species,Region) %>%
  summarise(count = n(),
            "tet(O)" = sum(`tet(O)`),
            "blaOXA-193" = sum(`blaOXA-193`),
            "50S_L22_A103V" = sum(`50S_L22_A103V`),
            "acrF" = sum(acrF),
            "blaEC" = sum(blaEC),
            "lin" = sum(lin),
            "fosX" = sum(`fosX`),
            "abc-f" = sum(`abc-f`),
            "mdsB" = sum(`mdsB`),
            "mdsA" = sum(`mdsA`),
            "aph(3'')-Ib" = sum(`aph(3'')-Ib`),
            "diary" = sum(`diary`),
            "livestock" = sum(`livestock`),
            "poultry" = sum(`poultry`),
            "seafood" = sum(`seafood`),
            "vege-fruit" = sum(`vege-fruit`))

# relabel region
count_df$Region <- factor(count_df$Region,levels = c("others","New England"), labels = c("Other Regions","New England"))


# Figure 1: Species Count comparison
fig1a <- count_df %>%
  group_by(Region) %>%
  summarise(count = sum(count)) %>%
  ggplot(aes(x = Region, y = count, fill = Region)) +
  geom_bar(stat = "identity")  +
  geom_text(aes(label = count), vjust = -0.5,size = 7) +
  scale_fill_manual(values = c("black","red")) +
  ggtitle("A") +
  # labs(fill = "Pathogens") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 40),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(size = 20),
        legend.title =  element_text(size = 25),
        legend.text  =  element_text(size = 25),
        legend.position = "none")
fig1a

# create the label of percentage in the plot
count_prop <- count_df %>%
  group_by(Region,Species) %>%
  summarise(count = sum(count)) %>%
  group_by(Region) %>%
  mutate(prop = round(count/sum(count),2)) %>%
  ungroup() %>%
  arrange(Species) %>%
  plyr::ddply( "Region",
               transform, 
               label_ypos=cumsum(prop) - 0.6*prop)

count_prop$Species <- factor(count_prop$Species, levels = levels(count_prop$Species)[c(4,3,2,1)])
count_prop$prop_string <- paste0(count_prop$prop*100,"%")
  
fig1b <- 
count_prop %>%
  ggplot(aes(x = Region,y = prop, fill = Species)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#c2a5cf","#a6dba0","#fdae61","#a6611a")) +
  geom_text(aes(y = label_ypos,label = prop_string),  color="black", size = 7) +
  ggtitle("B") + 
  labs(fill = "Pathogens") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 40),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(size = 20),
        legend.title =  element_text(size = 25),
        legend.text  =  element_text(size = 25))
fig1b

png("~/OneDrive - Brown University/PHP 2550/Final Project/Final Report/Figure1.png",width = 1400,height = 800)
ggarrange(fig1a,fig1b,nrow =1)
dev.off()

## Figure 2: correlation of food sources and FI cases
cases = get_all_cases(df,"others")
us_corr <- cor(cases$imputed)
cases = get_all_cases(df,"New England")
ne_corr <- cor(cases$imputed)

colnames(us_corr)[1] <- rownames(us_corr)[1] <- "FI cases"
colnames(ne_corr)[1] <- rownames(ne_corr)[1] <- "FI cases"
us_corr <- us_corr[1:6,1:6]
us_corr <- reshape::melt(us_corr)

fig7a <- ggplot(us_corr, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label = round(value,2)), color = "black", size = 7) +
  coord_fixed() +
  ggtitle("Other Regions") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20,angle = -90),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

ne_corr <- ne_corr[1:6,1:6]
ne_corr <- reshape::melt(ne_corr)

fig7b <- ggplot(ne_corr, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label = round(value,2)), color = "black", size = 7) +
  coord_fixed() +
  ggtitle("New England") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20,angle = -90),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

png("PDA/Project/Figure/Cases_food_corr.png",width = 1200,height = 800)
ggarrange(fig7a,fig7b, ncol = 2)
dev.off()

### Prediction Performance
### all FI cases
# US ARIMA model evaluation
cases = get_all_cases(df,"others")
results1 = list()
mode_list=c("none", "harmonic", "food", "HF", "AMR", "HAMR")
for(i in 1:6){
  results1[[i]] <- ARIMA_used(cases$imputed,mode=mode_list[i])
}

us_result_arima <- results1

us_aic <- rep(0,6)
us_long_term_error <- rep(0,6)
us_short_term_error <- rep(0,6)
us_long_term_corr <- rep(0,6)
us_short_term_corr <- rep(0,6)
for(i in 1:6){
  us_aic[i]<- results1[[i]]$model$aic
  us_long_term_error[i] <- sqrt(mean((results1[[i]]$prediction$mean - window(cases$imputed,start=c(2020,1))[,1])^2))
  us_short_term_error[i] <- sqrt(mean(head((results1[[i]]$prediction$mean-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  us_short_term_corr[i] <- cor(results1[[i]]$prediction$mean[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  us_long_term_corr[i] <- cor(results1[[i]]$prediction$mean,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  checkresiduals(us_result_arima[[i]]$model)
}


# New England ARIMA model evalation
cases = get_all_cases(df,"New England")
for(i in 1:6){
  results1[[i]] = ARIMA_used(cases$imputed,mode=mode_list[i])
}

ne_result_arima <- results1

ne_aic <- rep(0,6)
ne_long_term_error <- rep(0,6)
ne_short_term_error <- rep(0,6)
ne_long_term_corr <- rep(0,6)
ne_short_term_corr <- rep(0,6)
for(i in 1:6){
  ne_aic[i]<- results1[[i]]$model$aic
  ne_long_term_error[i] <- sqrt(mean((results1[[i]]$prediction$mean - window(cases$imputed,start=c(2020,1))[,1])^2))
  ne_short_term_error[i] <- sqrt(mean(head((results1[[i]]$prediction$mean-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  ne_short_term_corr[i] <- cor(results1[[i]]$prediction$mean[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  ne_long_term_corr[i] <- cor(results1[[i]]$prediction$mean,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  checkresiduals(ne_result_arima[[i]]$model)
}


# US GLM model evaluation
cases = get_all_cases(df,"others")
results2 = list()
mode_list=c("none", "harmonic", "food", "HF", "AMR", "HAMR")
for(i in 1:6){
  results2[[i]] <- TSGLM_used(cases$imputed, past_obs = c(1:2,12), mode=mode_list[i]) 
}
us_result_glm <- results2


us_aic_glm <- rep(0,6)
us_long_term_error_glm  <- rep(0,6)
us_short_term_error_glm  <- rep(0,6)
us_long_term_corr_glm <- rep(0,6)
us_short_term_corr_glm <- rep(0,6)

for(i in 1:6){
  us_aic_glm[i]<- AIC(results2[[i]]$model)
  us_long_term_error_glm[i] <- sqrt(mean((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2))
  us_short_term_error_glm[i] <- sqrt(mean(head((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  us_short_term_corr_glm[i] <- cor(results2[[i]]$prediction$pred[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  us_long_term_corr_glm[i] <- cor(results2[[i]]$prediction$pred,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  pit(us_result_glm[[i]]$model,main = "PIT Negative Binomial")
}


# New England GLM model evaluation
cases = get_all_cases(df,"New England")
results2 = list()
mode_list=c("none", "harmonic", "food", "HF", "AMR", "HAMR")
for(i in 1:6){
  results2[[i]] <- TSGLM_used(cases$imputed, past_obs = c(1), mode=mode_list[i]) 
}
ne_result_glm <- results2

ne_aic_glm <- rep(0,6)
ne_long_term_error_glm  <- rep(0,6)
ne_short_term_error_glm  <- rep(0,4)
ne_long_term_corr_glm <- rep(0,6)
ne_short_term_corr_glm <- rep(0,6)
for(i in 1:6){
  ne_aic_glm[i]<- AIC(results2[[i]]$model)
  ne_long_term_error_glm[i] <- sqrt(mean((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2))
  ne_short_term_error_glm[i] <- sqrt(mean(head((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  ne_short_term_corr_glm[i] <- cor(results2[[i]]$prediction$pred[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  ne_long_term_corr_glm[i] <- cor(results2[[i]]$prediction$pred,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  pit(ne_result_glm[[i]]$model,main = "PIT Negative Binomial")
}


# Figure 3: Forecast of 2020-2022 FI cases attributed to all four primary pathogens by SARIMA and GLM-based model
us_forcast_arima <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                               obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
us_forcast_arima$AR[121:155] <- as.numeric(us_result_arima[[1]]$prediction$mean)
us_forcast_arima$HAR[121:155] <- as.numeric(us_result_arima[[2]]$prediction$mean)
us_forcast_arima$FAR[121:155] <- as.numeric(us_result_arima[[3]]$prediction$mean)
us_forcast_arima$HARF[121:155] <- as.numeric(us_result_arima[[4]]$prediction$mean)
us_forcast_arima$AMR[121:155] <- as.numeric(us_result_arima[[5]]$prediction$mean)
us_forcast_arima$HAMR[121:155] <- as.numeric(us_result_arima[[6]]$prediction$mean)

us_forcast_glm <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                             obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
us_forcast_glm$AR[121:155] <- as.numeric(us_result_glm[[1]]$prediction$pred)
us_forcast_glm$HAR[121:155] <- as.numeric(us_result_glm[[2]]$prediction$pred)
us_forcast_glm$FAR[121:155] <- as.numeric(us_result_glm[[3]]$prediction$pred)
us_forcast_glm$HARF[121:155] <- as.numeric(us_result_glm[[4]]$prediction$pred)
us_forcast_glm$AMR[121:155] <- as.numeric(us_result_glm[[5]]$prediction$pred)
us_forcast_glm$HAMR[121:155] <- as.numeric(us_result_glm[[6]]$prediction$pred)

fig2a <- ggplot(us_forcast_arima,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_arima$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, ARIMA") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig2a

fig2b <- ggplot(us_forcast_glm,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_glm$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, GLM") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig2b


ne_forcast_arima <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                               obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
ne_forcast_arima$AR[121:155] <- as.numeric(ne_result_arima[[1]]$prediction$mean)
ne_forcast_arima$HAR[121:155] <- as.numeric(ne_result_arima[[2]]$prediction$mean)
ne_forcast_arima$FAR[121:155] <- as.numeric(ne_result_arima[[3]]$prediction$mean)
ne_forcast_arima$HARF[121:155] <- as.numeric(ne_result_arima[[4]]$prediction$mean)

ne_forcast_glm <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                             obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
ne_forcast_glm$AR[121:155] <- as.numeric(ne_result_glm[[1]]$prediction$pred)
ne_forcast_glm$HAR[121:155] <- as.numeric(ne_result_glm[[2]]$prediction$pred)
ne_forcast_glm$FAR[121:155] <- as.numeric(ne_result_glm[[3]]$prediction$pred)
ne_forcast_glm$HARF[121:155] <- as.numeric(ne_result_glm[[4]]$prediction$pred)
ne_forcast_glm$HARF[151] <- NA

fig2c <- ggplot(ne_forcast_arima,aes(x = time)) +
  geom_vline(xintercept = as.numeric(ne_forcast_arima$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("New England, ARIMA") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig2c

fig2d <- ggplot(ne_forcast_glm,aes(x = time)) +
  geom_vline(xintercept = as.numeric(ne_forcast_glm$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("New England, GLM") +
  ylim(c(0,80)) + 
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig2d

png("PDA/Project/Figure/Figure2.png",width = 1500,height = 800)
ggarrange(fig2a,fig2b,fig2c,fig2d,nrow =2,ncol = 2,common.legend = T,legend = "bottom")
dev.off()



### Species specific FI cases
### Salmonella Enterica in the U.S.
### ARIMA evaluation
cases = get_specific_cases(df,"others","Salmonella enterica")
results1 = list()
mode_list=c("none", "harmonic", "food", "HF")
for(i in 1:4){
  results1[[i]] <- ARIMA_used(cases$imputed,mode=mode_list[i])
}

us_result_arima_sal <- results1

us_aic_sal <- rep(0,4)
us_long_term_error_sal <- rep(0,4)
us_short_term_error_sal <- rep(0,4)
us_long_term_corr_sal <- rep(0,4)
us_short_term_corr_sal <- rep(0,4)
for(i in 1:4){
  us_aic_sal[i]<- results1[[i]]$model$aic
  us_long_term_error_sal[i] <- sqrt(mean((results1[[i]]$prediction$mean - window(cases$imputed,start=c(2020,1))[,1])^2))
  us_short_term_error_sal[i] <- sqrt(mean(head((results1[[i]]$prediction$mean-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  us_short_term_corr_sal[i] <- cor(results1[[i]]$prediction$mean[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  us_long_term_corr_sal[i] <- cor(results1[[i]]$prediction$mean,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  checkresiduals(us_result_arima_sal[[i]]$model)
}



### GLM evaluation
cases = get_specific_cases(df,"others","Salmonella enterica")
results2 = list()
mode_list=c("none", "harmonic", "food", "HF")
for(i in 1:4){
  results2[[i]] <- TSGLM_used(cases$imputed, past_obs = c(1:2,12), mode=mode_list[i]) 
}
us_result_glm_sal <- results2


us_aic_glm_sal <- rep(0,4)
us_long_term_error_glm_sal  <- rep(0,4)
us_short_term_error_glm_sal  <- rep(0,4)
us_long_term_corr_glm_sal <- rep(0,4)
us_short_term_corr_glm_sal <- rep(0,4)

for(i in 1:4){
  us_aic_glm_sal[i]<- AIC(results2[[i]]$model)
  us_long_term_error_glm_sal[i] <- sqrt(mean((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2))
  us_short_term_error_glm_sal[i] <- sqrt(mean(head((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  us_short_term_corr_glm_sal[i] <- cor(results2[[i]]$prediction$pred[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  us_long_term_corr_glm_sal[i] <- cor(results2[[i]]$prediction$pred,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  pit(us_result_glm_sal[[i]]$model,main = "PIT Negative Binomial")
}



### Listeria in the US
### ARIMA evaluation
cases = get_specific_cases(df,"others","Listeria monocytogenes")
results1 = list()
mode_list=c("none", "harmonic", "food", "HF")
for(i in 1:4){
  results1[[i]] <- ARIMA_used(cases$imputed,mode=mode_list[i])
}

us_result_arima_lis <- results1

us_aic_lis <- rep(0,4)
us_long_term_error_lis <- rep(0,4)
us_short_term_error_lis <- rep(0,4)
us_long_term_corr_lis <- rep(0,4)
us_short_term_corr_lis <- rep(0,4)
for(i in 1:4){
  us_aic_lis[i]<- results1[[i]]$model$aic
  us_long_term_error_lis[i] <- sqrt(mean((results1[[i]]$prediction$mean - window(cases$imputed,start=c(2020,1))[,1])^2))
  us_short_term_error_lis[i] <- sqrt(mean(head((results1[[i]]$prediction$mean-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  us_short_term_corr_lis[i] <- cor(results1[[i]]$prediction$mean[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  us_long_term_corr_lis[i] <- cor(results1[[i]]$prediction$mean,window(cases$imputed,start=c(2020,1))[,1])
}


for(i in 1:6){
  checkresiduals(us_result_arima_lis[[i]]$model)
}

### GLM evaluation
cases = get_specific_cases(df,"others","Listeria monocytogenes")
results2 = list()
mode_list=c("none", "harmonic", "food", "HF")
for(i in 1:4){
  results2[[i]] <- TSGLM_used(cases$imputed, past_obs = c(1:2,12), mode=mode_list[i]) 
}
us_result_glm_lis <- results2


us_aic_glm_lis <- rep(0,4)
us_long_term_error_glm_lis  <- rep(0,4)
us_short_term_error_glm_lis  <- rep(0,4)
us_long_term_corr_glm_lis <- rep(0,4)
us_short_term_corr_glm_lis <- rep(0,4)

for(i in 1:4){
  us_aic_glm_lis[i]<- AIC(results2[[i]]$model)
  us_long_term_error_glm_lis[i] <- sqrt(mean((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2))
  us_short_term_error_glm_lis[i] <- sqrt(mean(head((results2[[i]]$prediction$pred-window(cases$imputed,start=c(2020,1))[,1])^2,12)))
  us_short_term_corr_glm_lis[i] <- cor(results2[[i]]$prediction$pred[1:12],window(cases$imputed,start=c(2020,1))[1:12,1])
  us_long_term_corr_glm_lis[i] <- cor(results2[[i]]$prediction$pred,window(cases$imputed,start=c(2020,1))[,1])
}


## Model diagnosis
for(i in 1:6){
  pit(us_aic_glm_lis[[i]]$model,main = "PIT Negative Binomial")
}

# Figure 4: Forecast of 2020-2022 Salmonella enterica FI cases by SARIMA and GLM-based model
us_forcast_arima_sal <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                                   obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
us_forcast_arima_sal$AR[121:155] <- as.numeric(us_result_arima_sal[[1]]$prediction$mean)
us_forcast_arima_sal$HAR[121:155] <- as.numeric(us_result_arima_sal[[2]]$prediction$mean)
us_forcast_arima_sal$FAR[121:155] <- as.numeric(us_result_arima_sal[[3]]$prediction$mean)
us_forcast_arima_sal$HARF[121:155] <- as.numeric(us_result_arima_sal[[4]]$prediction$mean)

us_forcast_glm_sal <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                                 obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
us_forcast_glm_sal$AR[121:155] <- as.numeric(us_result_glm_sal[[1]]$prediction$pred)
us_forcast_glm_sal$HAR[121:155] <- as.numeric(us_result_glm_sal[[2]]$prediction$pred)
us_forcast_glm_sal$FAR[121:155] <- as.numeric(us_result_glm_sal[[3]]$prediction$pred)
us_forcast_glm_sal$HARF[121:155] <- as.numeric(us_result_glm_sal[[4]]$prediction$pred)


fig3a <- ggplot(us_forcast_arima_sal,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_arima_sal$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, Salmonella, ARIMA") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig3a

fig3b <- ggplot(us_forcast_glm_sal,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_glm_sal$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, Salmonella, GLM") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig3b

png("PDA/Project/Figure/US_Salmonella.png",width = 1500,height = 400)
ggarrange(fig3a,fig3b, ncol = 2,common.legend = T,legend = "bottom")
dev.off()

 
# Figure 5 Forecast of 2020-2022 Listeria monocytogenes FI cases by SARIMA and GLM-based model.
us_forcast_arima_lis <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                                  obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
us_forcast_arima_lis$AR[121:155] <- as.numeric(us_result_arima_lis[[1]]$prediction$mean)
us_forcast_arima_lis$HAR[121:155] <- as.numeric(us_result_arima_lis[[2]]$prediction$mean)
us_forcast_arima_lis$FAR[121:155] <- as.numeric(us_result_arima_lis[[3]]$prediction$mean)
us_forcast_arima_lis$HARF[121:155] <- as.numeric(us_result_arima_lis[[4]]$prediction$mean)

us_forcast_glm_lis <- data.frame(time = as.Date(as.yearmon(time(window(cases$imputed,start=c(2010,1))[,1]))),
                                 obs = as.numeric(window(cases$imputed,start=c(2010,1))[,1]))
us_forcast_glm_lis$AR[121:155] <- as.numeric(us_result_glm_lis[[1]]$prediction$pred)
us_forcast_glm_lis$HAR[121:155] <- as.numeric(us_result_glm_lis[[2]]$prediction$pred)
us_forcast_glm_lis$FAR[121:155] <- as.numeric(us_result_glm_lis[[3]]$prediction$pred)
us_forcast_glm_lis$HARF[121:155] <- as.numeric(us_result_glm_lis[[4]]$prediction$pred)


fig5a <- ggplot(us_forcast_arima_lis,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_arima_lis$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, Listeria, ARIMA") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig5a

fig5b <- ggplot(us_forcast_glm_lis,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_glm_lis$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAR, color = "Harmonic"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = HARF, color = "Harmonic + Sources"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, Listeria, GLM") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Harmonic', 'Sources','Harmonic + Sources'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77',"Harmonic" = "#e6ab02", "Sources" = "#0868ac", "Harmonic + Sources" = "#f03b20")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig5b

png("PDA/Project/Figure/US_Lister.png",width = 1500,height = 400)
ggarrange(fig5a,fig5b, ncol = 2,common.legend = T,legend = "bottom")
dev.off()



# Figure S3. The top 5 most common Isolation Sources of Foodborne Illness cases reported in New England, CA and TX during different seasons
fig6a <- ggplot(us_forcast_arima,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_arima$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = AMR, color = "AMR Genotype"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAMR, color = "Harmonic + AMR Genotype"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, ARIMA") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Sources', 'AMR Genotype','Harmonic + AMR Genotype'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77', "Sources" = "#0868ac", "AMR Genotype" = "#7570b3", "Harmonic + AMR Genotype" = "#e7298a")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig6a

fig6b <- ggplot(us_forcast_glm,aes(x = time)) +
  geom_vline(xintercept = as.numeric(us_forcast_glm$time[121]),
             linetype=2, colour="#08519c") +
  scale_x_date(date_labels = "%Y-%m-%d",date_breaks = "1 year") + 
  geom_line(aes(y = obs,color = "Observations"),size = 1) +
  geom_line(aes(y = AR,color = "No Covariate"),size = 1,alpha = 0.4) +
  geom_line(aes(y = FAR, color = "Sources"),size = 1) +
  geom_line(aes(y = AMR, color = "AMR Genotype"),size = 1,alpha = 0.4) +
  geom_line(aes(y = HAMR, color = "Harmonic + AMR Genotype"),size = 1,alpha = 0.4) +
  labs(y = "Detection Count") +
  ggtitle("Other Regions, GLM") +
  scale_color_manual(breaks=c('Observations', 'No Covariate','Sources', 'AMR Genotype','Harmonic + AMR Genotype'),
                     values = c("Observations" = "#252525",'No Covariate' = '#1b9e77', "Sources" = "#0868ac", "AMR Genotype" = "#7570b3", "Harmonic + AMR Genotype" = "#e7298a")) +
  theme(plot.title = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,size = 10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20))
fig6b

png("PDA/Project/Figure/US_AMR_Gene.png",width = 1500,height = 400)
ggarrange(fig6a,fig6b, ncol = 2,common.legend = T,legend = "bottom")
dev.off()


## Figure 2: correlation of food sources and FI cases
cases = get_all_cases(df,"others")
us_corr <- cor(cases$imputed)
cases = get_all_cases(df,"New England")
ne_corr <- cor(cases$imputed)

colnames(us_corr)[1] <- rownames(us_corr)[1] <- "FI cases"
colnames(ne_corr)[1] <- rownames(ne_corr)[1] <- "FI cases"
us_corr <- us_corr[1:6,1:6]
us_corr <- reshape::melt(us_corr)

fig7a <- ggplot(us_corr, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label = round(value,2)), color = "black", size = 7) +
  coord_fixed() +
  ggtitle("Other Regions") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20,angle = -90),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

ne_corr <- ne_corr[1:6,1:6]
ne_corr <- reshape::melt(ne_corr)

fig7b <- ggplot(ne_corr, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label = round(value,2)), color = "black", size = 7) +
  coord_fixed() +
  ggtitle("New England") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20,angle = -90),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

png("PDA/Project/Figure/Cases_food_corr.png",width = 1200,height = 800)
ggarrange(fig7a,fig7b, ncol = 2)
dev.off()



