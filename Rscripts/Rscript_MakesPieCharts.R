#getwd()
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
#Make sure Data Working Directory is in DataAnalysis
#install.packages("RColorBrewer")
library(RColorBrewer)
#to create list of all csv files
ResDocs<-list.files("Data/", pattern ="_ST_Data.csv")
#ResDocs<-list.files("Data/", pattern ="PubMLST_ST_Data.csv") #To only look at the PubMLST data

#to select specific set of colors 
myPalette = rep(sample(colorRampPalette(brewer.pal(12, "Paired"))(100),100))

#Creates a for loop that loops through each dataset
for (a in 1:length(ResDocs)){
  DataSet <- read.csv(paste0("Data/", ResDocs[a]), stringsAsFactors = F)
  DataSet<-DataSet[which(DataSet$ST!="Minor"),] #Oct 2021 Pleuni removed the minor category. It is not one ST, but many, if we count it as one it messes up our diversity calculations). 
  #gsub is to put in the chart, to identify it and remove the extension "ST_Data.csv"
  StudyName <-gsub("_ST_Data.csv",'',ResDocs[a], ignore.case = FALSE,perl = FALSE, fixed = FALSE,useBytes = FALSE)
  
  DataSet$ST2 <-gsub("[a-zA-Z]" , "", as.character(DataSet$ST)) #Get rid of the characters of the ST
  DataSet$ST2<-sub("^$", 999, DataSet$ST2)
  
  DrugList<-unique(DataSet$Drug)
  #to make pie chart for each data sets with susceptible population followed by resistant population
  
  colourCount = myPalette[as.numeric(DataSet$ST2[DataSet$Drug==DrugList[1]])%%100]   #Oct 2021 Pleuni changed this line  
  #we need to figure out how to assign color to character in ST column (eg "Minor" in Yamaji data sets, may be assign grey color to NA/character ST)
  
  for (i in 1:length(DrugList)){
    #creates png files
    png(paste0("Output/PieCharts/piechart_", StudyName, "_", DrugList[i], ".png"),  width = 350, height = 150, units='mm', res = 300)
        
    #creates two pie charts side by side
    par(mfrow=c(1,2))
    #to subset each data sets to slecect only susceptible population
    S = subset(DataSet, Drug == DrugList[i], select = "NumSusc")
    S<-S[S$NumSusc != 0, ] #get rids of overlaps (zeros), "!= "is used to get rid of the overlaps that are created due to no susceptible
    S <- as.numeric(unlist(S))
    
    #to subset each data sets to slecect only resistant population
    R = subset(DataSet, Drug == DrugList[i], select = "NumRes")
    R<-R[R$NumRes != 0, ] #get rids of overlaps (zeros), "!= "is used to get rid of the overlaps that are created due to no resistance
    R <- as.numeric(unlist(R))
    # "&" is used to combined two conditions so that it creates pie chart only if there are resistant population (meets both conditions)
    if (sum(S)!= 0 && sum(R)!= 0){
      #to include labels only if it is greater than 5% (0.05)
      label_vector = unique(DataSet$ST)
      proportion_S = S/sum(S)
      label_vector[which(proportion_S < 0.05)] = " " 
      #to plot pie chart
      pie(S, clockwise=TRUE, labels = label_vector, main = paste(StudyName,DrugList[i], "Susc"), col = colourCount)
      }
    #to create pie chart for resistant population
    if (sum(R)!= 0){
      label_vector = unique(DataSet$ST)
      proportion_R = R/sum(R)
      label_vector[which(proportion_R < 0.05)] = " "
      pie(R, clockwise=TRUE, labels = label_vector, main = paste(StudyName,DrugList[i], "Res"), col = colourCount)
      }
    
    dev.off()
  }
}

