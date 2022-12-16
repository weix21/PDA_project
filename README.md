## PHP 2550 Final Project: United States regional primary foodborne illness pathogen activity prediction using time series models with food sources 
Xin Wei & Chichun Tan, Fall 2022, Brown University School of Public Health

![](https://github.com/weix21/PDA_project/blob/main/figure/26911666558978_.pic_hd.jpg)

### Introduction 
In recent years, foodborne illness (FI) has become a concerning public health problem in
the U.S. Food poisoning could happen when one was unaware, lead to severe symptoms and even
threatened the victims’ life. Therefore, a FI outbreak could burden local public health system and cause
numerous economics loss. Therefore, developing a reliable forecasting tool for FI outbreaks was in high
demand which could enhance the ability of public health department to make a quick response to a
potential outbreak. We developed two time series models: Seasonal ARIMA and GLM-based time series regressions model, to capture the periodicity in FI pathogens' activities and implemented food sources information to improve FI outbreak predictions.


### Data

We used the isolation data from NCBI Pathogen Detection (https://www.ncbi.nlm.nih.gov/pathogens/). We pull out cases in six states in New England, and the rest cases from other regions in U.S, from the Species *Salmonella enterica*, *E.coli and Shigella*, *Campylobacter jejuni* and *Listeria monocytogenes*. A set of cleaned data are stored at [Processed Data](https://github.com/weix21/PDA_project/tree/main/data/Processed%20data).

### Model 

Seasonal ARIMA

![](https://github.com/weix21/PDA_project/blob/main/figure/SARIMA.png)

GLM-based Time Series Regression 

![](https://github.com/weix21/PDA_project/blob/main/figure/glm_based.png)


### Figure

![](https://github.com/weix21/PDA_project/blob/main/figure/Figure1.png)

*The detections of four primary FI pathogens from 2010 to 2022.*

![](https://github.com/weix21/PDA_project/blob/main/figure/Cases_food_corr.png)

*Associations between FI cases and food sources from 2010 to 2022*

![](https://github.com/weix21/PDA_project/blob/main/figure/Figure2.png)

*Forecast of 2020-2022 FI cases attributed to all four primary pathogens by SARIMA and
GLM-based model*

### Conclusion 

Overall, food source could improve the prediction performance for future FI
outbreaks. Instead of making prediction based on the periodicity of pathogen’s viability, food sources
data enhanced the ability for a prediction model in capturing the trends of FI cases. Our analysis showed
SARIMA and GLM had comparable performance on prediction, and both of them provided predictions
more reliable in a short term than a long term.

### Code

* Code for data preprocessing [Process Code](https://github.com/weix21/PDA_project/blob/main/code/Process_data.R)
* Functions for modelling and forecasting [Model Code](https://github.com/weix21/PDA_project/blob/main/code/analysis.R)
* Code for generating the figure [Figure Code](https://github.com/weix21/PDA_project/blob/main/code/Plot_Table_Generation.R)


