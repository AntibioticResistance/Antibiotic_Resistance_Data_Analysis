setwd("~/Documents/GitHub/Antibiotic_Resistance_Data_Analysis")

library(reshape2)
library(dplyr)
library(ggplot2)

#read "DivIndices.csv" 
DivIndices <- read.csv(file = "Output/DivIndices.csv")

#to remove the "_" in the name of the datasets and replace with empty space, 
#it is useful to fix the name of the dataset with labeller function in facet_grid
DivIndices$Dataset <- gsub("_", " ", DivIndices$Dataset)
DivIndices$Dataset <- gsub("-", " ", DivIndices$Dataset)

#to add asterisk instead of pvalue
DivIndices$star<-"n.s."
DivIndices$star[DivIndices$GSIpvalue > .05]  <- "n.s."
DivIndices$star[DivIndices$GSIpvalue <= .05]  <- " *"
DivIndices$star[DivIndices$GSIpvalue <= .01]  <- " **"
DivIndices$star[DivIndices$GSIpvalue <= .001] <- "***"

DatasetList = unique(DivIndices$Dataset)

#Remove rows with no resistance
DivIndices<-DivIndices[DivIndices$NumRes>1,]

DivIndices$FracRes <- DivIndices$NumRes / (DivIndices$NumRes + DivIndices$NumSus)

######################################################################################################## 
#to make barplot using ggplot

#selecct 6 required columns from dataframe; "DivIndices"
df1 <- data.frame(DivIndices$Dataset ,DivIndices$Drug,DivIndices$star, DivIndices$GSIvalueS, DivIndices$GSIvalueR, DivIndices$FracRes) 

#make a column with info about how high the p-value indicator should be above the bars
df1$heightofstars<- 0
for (i in 1:nrow(df1)){df1$heightofstars[i] = max(df1$DivIndices.GSIvalueS[i], df1$DivIndices.GSIvalueR[i])}

#merge two columns into one using melt function in order to plot side-by-side in ggplot
df2 <- melt(df1, id=c("DivIndices.Dataset", "DivIndices.Drug", "DivIndices.star", "heightofstars", "DivIndices.FracRes"))
#change name of the columns
colnames(df2)[1:7]<-c("Dataset","Drug","star", "heigthofstar", "FracRes", "Response", "SI_values")

#We need the stars only for half of the bars
df2$star[df2$Response == "DivIndices.GSIvalueR"] <- ""

#to change the order of the barplots
df2$Dataset<-factor(df2$Dataset, levels = 
                      c("Yamaji 1999", "Yamaji 2016", "Addams Sapper", "Kallonen BSAC", "Kallonen CUH", "Wurster" , "Manara", "Galloway" ))

DFForBarplots <- df2 %>%
  #arrange(factor(Dataset, levels = c("Yamaji 1999", "Yamaji 2016", "Addams Sapper", "Kallonen" , "Wurster" , "Manara", "Galloway" )), FracRes) %>%
  arrange(Dataset,FracRes,Drug) %>%
  # 3. Add order column of row numbers
  mutate(order = row_number())

#to make barplot using ggplot
png(paste0("Output/", "Figure2A_GiniSimpson", ".png"), width = 350, height = 150, units='mm', res = 300)
P<-ggplot(DFForBarplots) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill= "#F5F5F5",alpha= 0.5,
            colour = "black",lwd = 0.7) + 
  geom_bar(aes(x=as.factor(order), y=SI_values,fill=Response),stat='identity', position='dodge')+ 
  ylim(0,1)+
  ylab(" Gini Simpson Index")+ xlab("Drugs")+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Susceptible", "Resistant"))+ #to add color and legend manually
  facet_grid( ~ Dataset, scales = "free_x", space = "free", labeller = label_wrap_gen(10))+ 
  theme_bw() + 
  theme(panel.border = element_blank(),
        text = element_text(size=12),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.1, size = 10),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 9, face = "bold"))+
  scale_x_discrete(breaks= seq(2 , max(DFForBarplots$order), by = 2) ,
                   labels=DFForBarplots$Drug[seq(2 , max(DFForBarplots$order), by = 2)])+
  geom_text(aes(label=star, x=as.factor(order), y=heigthofstar), vjust=0., hjust = 0.2, size = 3)+
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(colour="white", size=10, 
                                                             face="bold"))

P
dev.off()


#to make barplot using ggplot
png(paste0("Output/", "Figure2B_FracRes", ".png"), width = 350, height = 150, units='mm', res = 300)

P<-ggplot(DFForBarplots[DFForBarplots$Response == "DivIndices.GSIvalueS",]) +
  geom_rect(aes( xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill= "#F5F5F5",alpha= 0.5,
            colour = "black",size = 0.7) + 
  geom_bar(aes(x=as.factor(order), y=FracRes),stat='identity', position='dodge')+ 
  ylim(0,1)+
  ylab("Fraction Resistant")+ xlab("Drugs")+
  facet_grid( ~ Dataset, scales = "free_x", space = "free", labeller = label_wrap_gen(10))+ 
  theme_bw() + 
  theme(panel.border = element_blank(),
        text = element_text(size=12),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.1, size = 10),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 9, face = "bold"))+
  scale_x_discrete(breaks= DFForBarplots$order[DFForBarplots$Response == "DivIndices.GSIvalueS"] ,
                   labels=DFForBarplots$Drug[DFForBarplots$Response == "DivIndices.GSIvalueS"])

P 

dev.off()

#########################################################################################################






