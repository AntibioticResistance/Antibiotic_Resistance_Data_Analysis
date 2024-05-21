setwd("~/Documents/GitHub/Antibiotic_Resistance_Data_Analysis")

library(RColorBrewer)
myPalette = brewer.pal(6, "Pastel1")
OverviewDatasets<-read.csv(file = "Output/OverviewDatasets.csv")


#creating list of unique drugs
ListOfDrugs <- unique(OverviewDatasets$Drug)
for (drug in ListOfDrugs){
  #to select specific drug with specific information 
  Pathogen <- subset(OverviewDatasets, Drug == drug, select = c("Dataset","Pathogen","Drug","GSIvalueR","GSIvalueS","GSIpvalue"))
  #creating list of each datasets
  ListOfDatasets <- unique(Pathogen$Dataset)
  for (i in 1:length(ListOfDatasets)){
    dataset = ListOfDatasets[i]
  png(filename=paste0("Output/Barplots/SingleBarplot_",drug,dataset,".png"), width = 150, height = 200, units='mm', res = 300)
    
  #to make plot 
  #using rbind function 
  R<-rbind(Pathogen$GSIvalueS[i],Pathogen$GSIvalueR[i])
  par(mfrow=c(1,1))
  par(mar = c(5.1, 4.1, 4.1, 7.1), xpd = TRUE)
  barplot(R,las=1,ylim=c(0,1),  col=c("#00BFC4","#F8766D"),cex.main=1.0, cex.lab = 1.0,
          main="Diversity of Resistant and Susceptible strains", sub=drug, 
          ylab="Gini Simpson Index",names.arg = dataset,
          beside=T)
  legend("topright", c("Susceptible","Resistant"), fill = c("#00BFC4","#F8766D"),inset = c(-.35, 0.45), bty = "n", cex=1)

    #to add level of significance
    if (Pathogen$GSIpvalue[i]<0.0005) {
      text(x=2, y=Pathogen$GSIvalueS[i] +0.02, "______", srt = 0)
      text(x=2, y=Pathogen$GSIvalueS[i] +0.04, "***", srt = 0)}
    else if (Pathogen$GSIpvalue[i]<0.005) {
      text(x=2, y=Pathogen$GSIvalueS[i] +0.02, "______", srt = 0)
      text(x=2,y= Pathogen$GSIvalueS[i] +0.04, "**", srt = 0)}
    else if (Pathogen$GSIpvalue[i]<0.05) {
      text(x=2, y=Pathogen$GSIvalueS[i] +0.02, labels = "______", srt = 0)
      text(x=2, y=Pathogen$GSIvalueS[i] +0.04, "*", srt = 0)}
    
  
  dev.off()
  } 
}



