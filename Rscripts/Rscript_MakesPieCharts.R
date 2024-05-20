#######################
####
## Prepare to make Pie charts
####
#######################

setwd("~/Documents/GitHub/Antibiotic_Resistance_Data_Analysis")
#Make sure Data Working Directory is in DataAnalysis
#install.packages("RColorBrewer")
library(RColorBrewer)
#create list of all csv files we will read
ResDocs<-list.files("Data/", pattern ="_ST_Data.csv")

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

#######################
####
## Loop through each dataset, and then for each dataset loop through all drugs 
####
#######################

for (a in 1:length(ResDocs)){
  DataSet <- read.csv(paste0("Data/", ResDocs[a]), stringsAsFactors = F) #read data
  DataSet$ST <- as.character(DataSet$ST) #make sure ST is read as character 
  ## List of ST to remove (because unclear if they are ST)
  DataSet<-DataSet[which(!DataSet$ST %in% ListToRemove),] #May 2024 Pleuni removed bunch of ST categories. 
  DataSet <- DataSet[order(DataSet$ST, decreasing = FALSE),] #change the order (easier for colors later)
  
  StudyName <-gsub("_ST_Data.csv",'',ResDocs[a], ignore.case = FALSE,perl = FALSE, fixed = FALSE,useBytes = FALSE)
  DrugList<-unique(DataSet$Drug) #get a list of the drugs in this dataset
  #to make pie chart for each data sets with susceptible population followed by resistant population
  
  for (i in 1:length(DrugList)){ ##Loop over all drugs for the dataset
    #creates png files
    png(paste0("Output/PieCharts/piechart_", StudyName, "_", DrugList[i], ".png"),  width = 350, height = 150, units='mm', res = 300)
        
    par(mfrow=c(1,2))   #creates two pie charts side by side
    #to subset each data sets to select only susceptible population
    S = subset(DataSet, Drug == DrugList[i], select = c("ST","NumSusc"))
    S<-S[S$NumSusc != 0, ] #get rids of zeros because we don't want them in the pie chart

    #to subset each data sets to slecect only resistant population
    R = subset(DataSet, Drug == DrugList[i], select = c("ST","NumRes"))
    R<-R[R$NumRes != 0, ] #get rids of zeros because we don't want them in the pie chart
    # "&" is used to combined two conditions so that it creates pie chart only if there are resistant population (meets both conditions)
    if (sum(S$NumSusc)!= 0 && sum(R$NumRes)!= 0){
      #to include ST labels only if it is greater than 5% (0.05)
      label_vector = S$ST
      proportion_S = S$NumSusc/sum(S$NumSusc)
      label_vector[which(proportion_S < 0.05)] = " " 
      #get the colors from myPalette
      color_vectorS = myPalette[which(names(myPalette) %in% S$ST)] 
      #to plot the pie chart
      pie(S$NumSusc, clockwise=TRUE, labels = label_vector, main = paste(StudyName,DrugList[i], "Susc"), col = color_vectorS)
      }
    #to create pie chart for resistant population
    if (sum(R$NumRes)!= 0){
      #to include ST labels only if it is greater than 5% (0.05)
      label_vector = R$ST
      proportion_R = R$NumRes/sum(R$NumRes)
      label_vector[which(proportion_R < 0.05)] = " "
      #get the colors from myPalette
      color_vectorR = myPalette[which(names(myPalette) %in% R$ST)]
      #to plot the pie chart
      pie(R$NumRes, clockwise=TRUE, labels = label_vector, main = paste(StudyName,DrugList[i], "Res"), col = color_vectorR)
      }
    
    dev.off() #close the png file
  }
}

