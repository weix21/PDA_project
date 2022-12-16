## PHP 2550 Final Project: United States regional primary foodborne illness pathogen activity prediction using time series models with food sources 
Xin Wei & Chichun Tan, Fall 2022, Brown University School of Public Health

![](https://github.com/weix21/PDA_project/blob/main/figure/26911666558978_.pic_hd.jpg)

### Introduction 
In recent years, foodborne illness (FI) has become a concerning public health problem in
the U.S. Food poisoning could happen when one was unaware, lead to severe symptoms and even
threatened the victims’ life. Therefore, a FI outbreak could burden local public health system and cause
numerous economics loss. Therefore, developing a reliable forecasting tool for FI outbreaks was in high
demand which could enhance the ability of public health department to make a quick response to a
potential outbreak. Therefore 


### Data

We used the isolation data from NCBI Pathogen Detection (https://www.ncbi.nlm.nih.gov/pathogens/). We pull out cases in six states in New England, as well as Califonia and Texa, from the Species *Salmonella enterica*, *E.coli and Shigella*, *Campylobacter jejuni* and *Listeria monocytogenes*. A set of cleaned data are stored at [Processed Data](https://github.com/weix21/PDA_project/tree/main/data/Processed%20data).

### Figure

![](https://github.com/weix21/PDA_project/blob/main/figure/Figure1.png)

*The detections of four primary FI pathogens from 2010 to 2022.*

![](https://github.com/weix21/PDA_project/blob/main/figure/Cases_food_corr.png)

*Associations between FI cases and food sources from 2010 to 2022*

![](https://github.com/weix21/PDA_project/blob/main/figure/Figure2.png)

*Forecast of 2020-2022 FI cases attributed to all four primary pathogens by SARIMA and
GLM-based model*

### Code

* Code for generating the figure [Figure Code](https://github.com/weix21/PDA_project/blob/main/code/Exploratory_analysis.R)
* Code for data preprocessing [Process Code](https://github.com/weix21/PDA_project/blob/main/code/Process_data.R)
* Code for modelling and forecasting [Model Code]()
