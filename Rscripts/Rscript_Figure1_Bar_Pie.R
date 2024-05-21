#getwd()
setwd("~/Documents/GitHub/Antibiotic_Resistance_Data_Analysis")
library(RColorBrewer)
library(ggplot2)
library(ggpubr) #to use ggarrange
#to create list of all csv files #Nov 2022 Only use Manara and Wurster
ResDoc<- c("Wurster_ST_Data.csv", "Manara_ST_Data.csv")

#Datawrangling is done in another R script / Figure 1 is made after fig 2! 
source("RScripts/Rscript_Figure2_GSI_Datasets.R") 

#Make list of all ST so that we can assign consistent colors
ListST <- c()
for (a in 1:length(ResDocs)){
  DataSet <- read.csv(paste0("Data/", ResDocs[a]), stringsAsFactors = F)
  ListST <- c(ListST, unique(DataSet$ST))
}
ListST <- sort(unique(ListST))
#I want to remove a bunch of ST that are not really sequence types
ListToRemove <- c("-", "Minor", "ND " ,  "OC39" , "U1"  ,  "U10"  , "U11",   "U12",   "U13" ,  "U14" ,  "U15" ,  "U2"  , 
                  "U3"  ,  "U4" ,   "U5"  ,  "U6"  ,  "U7" ,   "U8"  ,  "U9")
ListST <- ListST[which(!ListST%in%ListToRemove)]
length(ListST) ##We have 324 sequence types

#to create specific set of colors 
set.seed(1) #so that if we run it again the same colors assigned
myPalette  = sample(colorRampPalette(brewer.pal(8, "Set2"))(length(ListST)),replace = FALSE)
names(myPalette) <- ListST #set names, useful later for consistent colors

for (X in 1:3){
  if (X == 1) {ResDoc = "Wurster_ST_Data.csv";  DrugName = "Penicillin"}
  if (X == 2) {ResDoc = "Wurster_ST_Data.csv"; DrugName = "Oxacillin"}
  if (X == 3) {ResDoc = "Kallonen_CUH_ST_Data.csv" ; DrugName = "Tobramycin"}
  
  #Creates a for loop that loops through each dataset
  DataSet <- read.csv(paste0("Data/", ResDoc), stringsAsFactors = F)
  DataSet$ST <- as.character(DataSet$ST) #make sure ST is read as character 
  ## List of ST to remove (because unclear if they are ST)
  DataSet<-DataSet[which(!DataSet$ST %in% ListToRemove),] #May 2024 Pleuni removed bunch of ST categories. 
  DataSet <- DataSet[order(DataSet$ST, decreasing = FALSE),] #change the order (easier for colors later)
  
  #gsub is to put in the chart, to identify it and remove the extension "ST_Data.csv"
  StudyName <-gsub("_ST_Data.csv",'',ResDoc, ignore.case = FALSE,perl = FALSE, fixed = FALSE,useBytes = FALSE)
  if (StudyName == "Kallonen_CUH") StudyName <- "Kallonen CUH"  
  
  DataSet$ST2 <-gsub("[a-zA-Z]" , "", as.character(DataSet$ST)) #Get rid of the characters of the ST
  DataSet$ST2<-sub("^$", 999, DataSet$ST2)
  
  ###FIRST MAKE THE BARPLOT 
  P<-ggplot(DFForBarplots[DFForBarplots$Dataset == StudyName & DFForBarplots$Drug == DrugName,]) +
    #geom_rect(aes( xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill= "#F5F5F5",alpha= 0.5,
    #           colour = "black",size = 0.7) + 
    geom_bar(aes(x=as.factor(order), y=SI_values,fill=Response),stat='identity', position='dodge')+ 
    ylim(0,1)+
    ylab("Gini Simpson Index")+ xlab("")+
    scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Susceptible", "Resistant"))+ #to add color and legend manually
    #facet_grid( ~ Dataset, scales = "free_x", space = "free", labeller = label_wrap_gen(10))+ 
    theme_classic() + 
    scale_x_discrete(breaks= seq(1 , max(DFForBarplots$order), by = 1) , labels=rep(c("Susc.","Res."), max(DFForBarplots$order)/2))
  P <- P + geom_text(aes(label=star, x=as.numeric(as.factor(order))+0.3, y=heigthofstar), vjust=-0.2, hjust = -0.1, size = 4 )
  P <- P + theme(legend.position = "none") + ggtitle(paste(DrugName, "(", StudyName, ")")) + 
    theme(plot.title = element_text(size = 10))
  
  ### Create the first pie chart (SUSCEPTIBLE)
  SubData = subset(DataSet, Drug == DrugName)
  
  Splot = ggplot(SubData, aes(x=Drug, y=NumSusc, fill=ST2))+
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0)+
    scale_fill_manual(values=myPalette)+
    theme_minimal()+ #we changed this from blank_theme to theme_bw() because we get an error
    theme(legend.position = "none")+
    theme(axis.text.x = element_blank(), 
          axis.title.y=element_blank(), 
          axis.title.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() )+
    labs(title = "Diversity susceptible strains")+
    theme(plot.title = element_text (size = 10, hjust = 0.5))
  
  ### Create the second pie chart (RESISTANT)
  
  Rplot = ggplot(SubData, aes(x=Drug, y=NumRes, fill=ST2))+
    geom_bar(width = 1, stat = "identity") + 
    scale_y_continuous(labels = NULL)+
    coord_polar("y", start=0)+
    scale_fill_manual(values=myPalette)+
    theme_minimal()+ #we changed this from blank_theme to theme_bw() because we get an error
    theme(legend.position = "none")+
    theme(axis.text.x = element_blank(), 
          axis.title.y=element_blank(),
          axis.title.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() )+
    labs(title = "Diversity resistant strains")+
    theme(plot.title = element_text (size = 10, hjust = 0.5))
  
  
  if (ResDoc == "Wurster_ST_Data.csv" & DrugName == "Penicillin"){
    WuPen_Pplot <- P ; WuPen_Splot <- Splot ; WuPen_Rplot<- Rplot}
  if (ResDoc == "Wurster_ST_Data.csv" & DrugName == "Oxacillin"){
    WuOxa_Pplot <- P ; WuOxa_Splot <- Splot ; WuOxa_Rplot<- Rplot}
  if (ResDoc == "Kallonen_CUH_ST_Data.csv" & DrugName == "Tobramycin"){
    KalTob_Pplot <- P ; KalTob_Splot <- Splot ; KalTob_Rplot<- Rplot}
}

#### Arrange the 9 figures into one: 

png("Output/Fig_1_Bar_PieCharts.png", width = 180, height = 180, units='mm', res = 300)
ggarrange(WuPen_Pplot, WuPen_Splot, WuPen_Rplot, 
          WuOxa_Pplot, WuOxa_Splot, WuOxa_Rplot, 
          KalTob_Pplot, KalTob_Splot, KalTob_Rplot,  
          ncol = 3, nrow = 3, labels = c("A","B","C","D","E","F", "G","H", "I" ))

dev.off()

