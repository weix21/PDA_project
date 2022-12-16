## PHP 2550 Final Project: United States regional primary foodborne illness pathogen activity prediction using time series models with food sources 
Xin Wei & Chichun Tan, Fall 2022, Brown University School of Public Health

![](https://github.com/weix21/PDA_project/blob/main/figure/26911666558978_.pic_hd.jpg)
### Introduction 

In recent years, foodborne illness has become a concerning public health problem in the U.S. A large proportion of foodborne illness cases were attributed to biological reasons, for example, different types of pathogens. According to the recent reports from the foodborne disease outbreak monitoring and coordinating systems in CDC and FDA , Salmonella, Escherichia coli (usually called E.coli, especially the type produces Shiga Toxin), Listeria monocytogenes and Campylobactor have been the four primary pathogen since the beginning of this century. Capturing the potential patterns of the outbreaks for these four pathogens could help public health workers prevent outbreaks and respond to emergency more effectively. As residents in New England, we concern about food safety problems in this area all the time. The states in New England share similar climate and dietary habit, indicating that they might take similar risk in foodborne illness issues. However, few studies discussed about the foodborne illness in New England as an entity. In this project, we will use data from Pathogen Detection in National Center for Biotechnology Information (NCBI) as a primary source to review and investigate the patterns of foodborne illness due to the four priority pathogens.

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
