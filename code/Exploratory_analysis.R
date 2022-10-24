## Xin and Chichun
## Exploratory analysis
## Date: 10/23/2022

## Load processed data
load("/Users/weixin/Downloads/study/Brown/PDA/Project/Processed_data_10_23.RData")


## Explore the realtionship between region and number of cases on species
## Make figure 1
df %>% select(Region, X.Organism.group) %>% table()
df %>% ggplot(aes(X.Organism.group)) +
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme_minimal() +
  coord_flip() +
  facet_wrap(~as.factor(Region), nrow=3)
ggsave("figure_1.pdf")



## Explore the realtionship between region, season and number of cases on species
## Make figure 2
df %>% ggplot(.,aes(x=factor(Season, levels = c("Spring","Summer","Autumn","Winter")))) + 
  geom_bar(stat="count", width=0.7)+
  xlab("Season")+ 
  theme(axis.text.x = element_text(size=8, angle = 30)) +
  facet_grid(as.factor(Region) ~ as.factor(X.Organism.group))
ggsave("figure_2.pdf", width = 7, height = 5)


## Explore the realtionship between region, day and number of cases on species
## Make figure 3.1
df %>% ggplot(.,aes(x=Date,y=nDate)) + 
  geom_line(colour = 'blue') + 
  theme(axis.text.x = element_text(size=8, angle = 45)) +
  facet_grid(as.factor(Region) ~ as.factor(X.Organism.group))
ggsave("figure_3_1.pdf", width = 7, height = 5)


## Explore the realtionship between region, quarter and number of cases on species
## Make figure 3.2
df %>% ggplot(.,aes(x=Quarter,y=nQuarter)) + 
  geom_line(colour = 'blue')  + 
  theme(axis.text.x = element_text(size=7, angle = 30)) + 
  facet_grid(as.factor(Region) ~ as.factor(X.Organism.group))
ggsave("figure_3_2.pdf", width = 7, height = 5)


## Identify the common sources of FI
Region_type <- levels(df$Region)
Organism_group <- levels(df$X.Organism.group)

p <- NULL
for(i in 1:3){
  for(j in 1:4){
    p[[4*(i-1)+j]] <- df %>% drop_na(Isolation.source) %>% 
      filter(X.Organism.group == Organism_group[j], Region == Region_type[i])%>%    
      mutate(rank = dense_rank(-nIsolation.source)) %>% arrange(desc(nIsolation.source)) %>% filter(rank<=3) %>%
      ggplot(.,aes(x=factor(Isolation.source,levels=rev(unique(Isolation.source))))) + geom_bar(stat="count", width=0.7) +
      labs(x = "")
  }
}

commonplot <- ggarrange(plotlist=p, nrow = 3, ncol=4)