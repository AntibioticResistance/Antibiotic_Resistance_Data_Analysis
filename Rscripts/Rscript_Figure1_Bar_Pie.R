#getwd()
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
library(RColorBrewer)
library(ggpubr) #to use ggarrange
#to create list of all csv files #Nov 2022 Only use Manara and Wurster
ResDoc<- c("Wurster_ST_Data.csv", "Manara_ST_Data.csv")

#Datawrangling is done in another R script / Figure 1 is made after fig 2! 
source("RScripts/Rscript_Figure2_GSI_Datasets.R") 
#Now we have the data in DFForBarplots
#DFForBarplots<-DFForBarplots[DFForBarplots$Dataset%in%c("Wurster", "Manara"),]
#to select specific set of colors 

set.seed(30)
myPalette = rep(sample(colorRampPalette(brewer.pal(12, "Paired"))(100),100))


for (X in 1:3){
  if (X == 1) {ResDoc = "Wurster_ST_Data.csv";  DrugName = "Penicillin"}
  if (X == 2) {ResDoc = "Wurster_ST_Data.csv"; DrugName = "Oxacillin"}
  if (X == 3) {ResDoc = "Manara_ST_Data.csv" ; DrugName = "Oxacillin"}
  if (X == 3) {ResDoc = "Galloway_ST_Data.csv" ; DrugName = "Vancomycin"}
  if (X == 3) {ResDoc = "Kallonen_CUH_ST_Data.csv" ; DrugName = "Tobramycin"}
  
#Creates a for loop that loops through each dataset
DataSet <- read.csv(paste0("Data/", ResDoc), stringsAsFactors = F)
DataSet<-DataSet[which(DataSet$ST!="Minor"),] #Oct 2021 Pleuni removed the minor category. It is not one ST, but many, if we count it as one it messes up our diversity calculations). 
DataSet<-DataSet[which(DataSet$ST!="-"),]
#gsub is to put in the chart, to identify it and remove the extension "ST_Data.csv"
StudyName <-gsub("_ST_Data.csv",'',ResDoc, ignore.case = FALSE,perl = FALSE, fixed = FALSE,useBytes = FALSE)
if (StudyName == "Kallonen_CUH") StudyName <- "Kallonen CUH"  

DataSet$ST2 <-gsub("[a-zA-Z]" , "", as.character(DataSet$ST)) #Get rid of the characters of the ST
DataSet$ST2<-sub("^$", 999, DataSet$ST2)
  
#to make pie chart for each data sets with susceptible population followed by resistant population
colourCount = myPalette[1+as.numeric(DataSet$ST2[DataSet$Drug==DrugName])%%100]   #Oct 2021 Pleuni changed this line  
#we need to figure out how to assign color to character in ST column (eg "Minor" in Yamaji data sets, may be assign grey color to NA/character ST)

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
    
    #creates two pie charts side by side
    SubData = subset(DataSet, Drug == DrugName)
    
    Splot = ggplot(SubData, aes(x=Drug, y=NumSusc, fill=ST2))+
      geom_bar(width = 1, stat = "identity") + 
      coord_polar("y", start=0)+
      scale_fill_manual(values=colourCount)+
      theme_minimal()+ #we changed this from blank_theme to theme_bw() because we get an error
      theme(legend.position = "none")+
      theme(axis.text.x = element_blank(), 
            axis.title.y=element_blank(), 
            axis.title.x=element_blank(), 
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank() )+
      labs(title = "Diversity susceptible strains")+
      theme(plot.title = element_text (size = 10, hjust = 0.5))
    
    Rplot = ggplot(SubData, aes(x=Drug, y=NumRes, fill=ST2))+
      geom_bar(width = 1, stat = "identity") + 
      scale_y_continuous(labels = NULL)+
      coord_polar("y", start=0)+
      scale_fill_manual(values=colourCount)+
      theme_minimal()+ #we changed this from blank_theme to theme_bw() because we get an error
      theme(legend.position = "none")+
      theme(axis.text.x = element_blank(), 
            axis.title.y=element_blank(),
            axis.title.x=element_blank(), 
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank() )+
      labs(title = "Diversity resistant strains")+
      theme(plot.title = element_text (size = 10, hjust = 0.5))
    
    
    #To do: remove labels, NOV2022
    # less important: add better labels NOV2022
    
    if (ResDoc == "Wurster_ST_Data.csv" & DrugName == "Penicillin"){
      WuPen_Pplot <- P ; WuPen_Splot <- Splot ; WuPen_Rplot<- Rplot}
    if (ResDoc == "Wurster_ST_Data.csv" & DrugName == "Oxacillin"){
      WuOxa_Pplot <- P ; WuOxa_Splot <- Splot ; WuOxa_Rplot<- Rplot}
    if (ResDoc == "Manara_ST_Data.csv" & DrugName == "Oxacillin"){
      MaPen_Pplot <- P ; MaPen_Splot <- Splot ; MaPen_Rplot<- Rplot}
    if (ResDoc == "Kallonen_CUH_ST_Data.csv" & DrugName == "Tobramycin"){
      KalTob_Pplot <- P ; KalTob_Splot <- Splot ; KalTob_Rplot<- Rplot}
}
  

png("Output/Fig_1_Bar_PieCharts.png", width = 180, height = 180, units='mm', res = 300)
    ggarrange(WuPen_Pplot, WuPen_Splot, WuPen_Rplot, 
              WuOxa_Pplot, WuOxa_Splot, WuOxa_Rplot, 
              KalTob_Pplot, KalTob_Splot, KalTob_Rplot,  
              ncol = 3, nrow = 3, labels = c("A","B","C","D","E","F", "G","H", "I" ))
              
dev.off()

