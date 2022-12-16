# Import required packages 
library(ggpubr)
library(lme4)
library(MASS)
library(AER)
library(tidyverse)
library(imputeTS)
library(forecast)
library(tscount)



# GLM based time series analysis
## cases: time series of monthly count of cases with covariates 
## train: last time point of training set
## test: First time point of test set
## prepoints: How many time points ahead will be used for prediction
## seed: fix random seed
## past_obs: which past observations will be included for autoregression
## past_mean: which past distributions' mean will be included for autoregression
## link: link function
## distr: conditional distribution
## K: period of harmonic term (Fourier transformation)
## mod: "none", "harmonic", "food", "HF", "AMR", "HAMR"
## None covariates included; only include harmonic term; only include food source; 
## include harmonic term with food source; only include AMR genotype; include harmonic term with AMR genotype;
TSGLM_used <- function(cases, train = c(2019,12), test = c(2020,1), prepoints=35, seed=1234,
                       past_obs = NULL,past_mean = NULL, link= "log", distr="nbinom",
                       K=1,mode=c("none", "harmonic", "food", "HF", "AMR", "HAMR")){
  
  set.seed(seed)
  cases1 <- window(cases,end=train)
  cases2 <- window(cases,start=test)
  cases2 <- head(cases2,prepoints)
  mode <- switch(mode,"none" = 0, "harmonic" = 1, "food" = 2, "HF" = 3, "AMR" = 4, "HAMR" = 5)
  
  if(mode==0){
    print("none additional terms")
    m1 <- tsglm(cases1[,1], model = list(past_obs = past_obs, past_mean = past_mean), 
                link = link, distr = distr)
    
    prediction <- predict(m1,n.ahead = prepoints)
    
  }
  if(mode==1){
    print("additional fourier terms")
    m1 <- tsglm(cases1[,1], model = list(past_obs = past_obs, past_mean = past_mean), 
                xreg =  cbind(fourier(cases1[,1],K=K)), link = link, distr = distr)
    
    prediction <- predict(m1,n.ahead = prepoints,
                          newxreg = cbind(fourier(cases2[,1],K=K)))
  }
  if(mode==2){
    print("additional food terms")
    m1 <- tsglm(cases1[,1], model = list(past_obs = past_obs, past_mean = past_mean), 
                xreg =  cbind(cases1[,2:6]), link = link, distr = distr)
    
    prediction <- predict(m1,n.ahead = prepoints,
                          newxreg = cbind(cases2[,2:6]))
  }
  if(mode==3){
    print("additional fourier&food terms")
    m1 <- tsglm(cases1[,1], model = list(past_obs = past_obs, past_mean = past_mean), 
                xreg =  cbind(fourier(cases1[,1],K=K),cases1[,2:6]), link = link, distr = distr)
    
    prediction <- predict(m1,n.ahead = prepoints,
                          newxreg = cbind(fourier(cases2[,1],K=K),cases2[,2:6]))
  }
  if(mode==4){
    print("additional AMR terms")
    m1 <- tsglm(cases1[,1], model = list(past_obs = past_obs, past_mean = past_mean), 
                xreg =  cbind(cases1[,7:10]), link = link, distr = distr)
    
    prediction <- predict(m1,n.ahead = prepoints,
                          newxreg = cbind(cases2[,7:10]))
  }
  if(mode==5){
    print("additional fourier&AMR terms")
    m1 <- tsglm(cases1[,1], model = list(past_obs = past_obs, past_mean = past_mean), 
                xreg =  cbind(fourier(cases1[,1],K=K),cases1[,7:10]), link = link, distr = distr)
    
    prediction <- predict(m1,n.ahead = prepoints,
                          newxreg = cbind(fourier(cases2[,1],K=K),cases2[,7:10]))
  }   
  return(list(model = m1, prediction = prediction))
}



## SARIMA model 
## cases: time series of monthly count of cases with covariates 
## train: last time point of training set
## test: First time point of test set
## prepoints: How many time points ahead will be used for prediction
## seed: fix random seed
## K: period of harmonic term (Fourier transformation)
## mod: "none", "harmonic", "food", "HF", "AMR", "HAMR"
## None covariates included; only include harmonic term; only include food source; 
## include harmonic term with food source; only include AMR genotype; include harmonic term with AMR genotype;
ARIMA_used <- function(cases, train=c(2019,12), test=c(2020,1), prepoints=35, seed=1234, 
                       K=1, mode=c("none", "harmonic", "food", "HF", "AMR", "HAMR")){
  
  set.seed(seed)
  
  cases1 <- window(cases,end=train)
  cases2 <- window(cases,start=test)
  cases2 <- head(cases2,prepoints)
  mode <- switch(mode,"none" = 0, "harmonic" = 1, "food" = 2, "HF" = 3, "AMR" = 4, "HAMR" = 5)
  
  if(mode==0){
    print("none additional terms")
    m1 <- auto.arima(cases1[,1],seasonal =TRUE)
    predict <- forecast(m1, h=prepoints)
    
  }
  if(mode==1){
    print("additional fourier terms")
    m1 <- auto.arima(cases1[,1],xreg = cbind(fourier(cases1[,1],K=K)),seasonal =TRUE)
    predict <- forecast(m1, xreg = cbind(fourier(cases2[,1],K=K)), h=prepoints)
    
  }
  if(mode==2){
    print("additional food terms")
    m1 <- auto.arima(cases1[,1],xreg = cbind(unclass(cases1[,2:6])),seasonal =TRUE)
    predict <- forecast(m1, xreg = cbind(unclass(cases2[1:prepoints,2:6])), h=prepoints)
  }
  if(mode==3){
    print("additional fourier&food terms")
    m1 <- auto.arima(cases1[,1],xreg = cbind(fourier(cases1[,1],K=K),unclass(cases1[,2:6])),seasonal =TRUE)
    predict <- forecast(m1, xreg = cbind(fourier(cases2[,1],K=K),unclass(cases2[,2:6])), h=prepoints)
  }
  if(mode==4){
    print("additional AMR terms")
    m1 <- auto.arima(cases1[,1],xreg = cbind(unclass(cases1[,7:10])),seasonal =TRUE)
    predict <- forecast(m1, xreg = cbind(unclass(cases2[1:prepoints,7:10])), h=prepoints)
  }
  if(mode==5){
    print("additional fourier&AMR terms")
    m1 <- auto.arima(cases1[,1],xreg = cbind(fourier(cases1[,1],K=K),unclass(cases1[,7:10])),seasonal =TRUE)
    predict <- forecast(m1, xreg = cbind(fourier(cases2[,1],K=K),unclass(cases2[,7:10])), h=prepoints)
  }
  return(list(model=m1,prediction=predict))
}



# Obtain cases for specific species
get_specific_cases <- function(df, region, species){
  new_df <- df %>% dplyr::select(c("Species","Collection_YM","Region","tet(O)",
                                   "blaOXA-193","50S_L22_A103V","acrF","blaEC",
                                   "mdtM","lin","fosX","abc-f","mdsB","mdsA","aph(3'')-Ib",
                                   "diary","livestock","poultry","seafood","vege-fruit")) %>%
    group_by(Collection_YM,Region,Species) %>%
    summarise(across(everything(), sum)) 
  
  new_df2 <- df %>% dplyr::select(c("Species","Collection_YM","Region")) %>%
    group_by(Collection_YM,Region,Species) %>%
    summarise(count=n()) 
  
  new_df <- merge(new_df,new_df2)
  
  new_df <- new_df %>% filter(Region == region & Species == species)
  
  missing_df <- data.frame(Collection_YM = sort(unique(df$Collection_YM)))
  
  new_df <- full_join(new_df,missing_df)
  
  new_df <- new_df[order(new_df$Collection_YM),]
  
  new_df <- new_df %>% dplyr::select(c("count","diary","livestock","poultry","seafood","vege-fruit",
                                       "tet(O)","blaEC","lin","mdsA")) 
  
  raw_cases <- ts(new_df,start=c(2010, 1), end=c(2022, 11), frequency=12)
  if(sum(is.na(raw_cases))!=0){
    cases <- round(na_interpolation(raw_cases, option = "linear"))
  }else{
    cases <- raw_cases
  }
  return(list(raw=raw_cases,imputed=cases))
}



# Obtain overall FI cases
get_all_cases <- function(df,region){
  new_df <- df %>% dplyr::select(c("Collection_YM","Region","tet(O)",
                                   "blaOXA-193","50S_L22_A103V","acrF","blaEC",
                                   "mdtM","lin","fosX","abc-f","mdsB","mdsA","aph(3'')-Ib",
                                   "diary","livestock","poultry","seafood","vege-fruit")) %>%
    group_by(Collection_YM,Region) %>%
    summarise(across(everything(), sum)) 
  
  new_df2 <- df %>% dplyr::select(c("Collection_YM","Region")) %>%
    group_by(Collection_YM,Region) %>%
    summarise(count=n()) 
  
  new_df <- merge(new_df,new_df2)
  
  new_df <- new_df %>% filter(Region == region)
  
  missing_df <- data.frame(Collection_YM = sort(unique(df$Collection_YM)))
  new_df <- full_join(new_df,missing_df)
  
  new_df <- new_df[order(new_df$Collection_YM),]
  
  new_df <- new_df %>% dplyr::select(c("count","diary","livestock","poultry","seafood","vege-fruit",
                                       "tet(O)","blaEC","lin","mdsA")) 
  
  raw_cases <- ts(new_df,start=c(2010, 1), end=c(2022, 11), frequency=12)
  if(sum(is.na(raw_cases))!=0){
    cases <- round(na_interpolation(raw_cases, option = "linear"))
  }else{
    cases <- raw_cases
  }
  return(list(raw=raw_cases,imputed=cases))
}
