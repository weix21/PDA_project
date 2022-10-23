library(tidyverse)
library(forcats)



setwd("/Users/weixin/Downloads/study/Brown/PDA/Project/data")

data_dir <-  c("Salmonella enterica.csv","E.coli and Shigella.csv","Campylobacter jejuni.csv","Listeria monocytogenes.csv")
data_list  <-  lapply(data_dir, read.csv)
names(data_list) <- strsplit(data_dir,split=".csv")

data = do.call(rbind.data.frame, data_list)

data$Location <- fct_collapse(
  data$Location,
  'CT'= c("USA:CT", "USA: CT", "USA: Connecticut, Greenwich", "USA: Connecticut"),
  'MA'= c("USA:MA", "USA: MA", "USA: Massachusetts", "USA: Cambridge, MA, MIT", "USA:Boston", "USA: Massachusetts, Boston", 
          "USA: Boston", "USA: Cambridge, MA", "USA: Boston, MA", "USA: MA, Boston", "USA:Massachusetts, Beverly", "USA: Ipswich, MA"),
  'ME'= c("USA:ME", "USA: ME", "USA: Maine", "USA:Maine"),
  'NH'= c("USA:NH", "USA: New Hampshire", "USA: New Hampshire,Rindge", "USA: NH"),
  'RI'= c("USA:RI", "USA: Rhode Island", "USA: RI"),
  'VT'= c("USA:VT", "USA: VT")) 


data %>% select(Location, X.Organism.group) %>% table()

data$Date <- as.Date(substr(data$Create.date,1,10),format="%Y-%m-%d")
