
getwd()
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
#WST1 <-read.csv
#to read the .csv file
Staph = read.csv(file = "Data/PatientLevel_Data/Wurster_S_aureus_Patient_Data_Revised.csv", header=TRUE, sep = ",")
#to check if the packages we need are present
library(ggplot2)
library(reshape2)
library(dplyr)
#install.packages(cowplot)  
library(colorspace)
library(grDevices)
#install.packages(colormap)
library(data.table)
library(RColorBrewer)
#WURSTER = name of the author
#to extract the cloumns we require; here we only require Sequence Type and Drug Names
#16th colum = Sequence Type
#6-15 = Drug Names
WURSTER = select(Staph, 16, 6:15)
#to get Sequence Type(ST) from .csv file
ST = Staph$Sequence.Type
#make list of the Sequence Type (each ST is represented only once and sorted)
UniqeStrains = sort(unique(ST))
UniqeStrains
drugnames = names(WURSTER)[2:11]#refers two the second through 11th, this is where the names of drug starts
#to create new dataframe for Simpson's Diversity Index (SI)
S_count <- data.frame("UniqueStrains" = UniqeStrains, "1" = 0, "2" = 0, "3" = 0, "4" = 0, "5" = 0, "6" = 0, "7" = 0, "8" = 0, "9" = 0, "10" = 0)
names(S_count)[2:11] <- drugnames #to assign drug names to the column
R_count = S_count
#to check drug name corresponding to the column
#drugColumn = 4
#the for loop for all the drugs starts here, therefore make sure you indentate for everything till the first for loop closes
for (drugColumn in  2:11) {
  drugName = names(WURSTER)[drugColumn]
  drugR = ST[which(WURSTER[, drugColumn]=="R")]
  drugS = ST[which(WURSTER[, drugColumn]=="S")]
  #to make empty list to enlist  all the ST
  drugRcount = c()
  drugScount = c()
  #for loop to count number of each Sequence Type
  for (i in UniqeStrains) {
    drugRcount<-c(drugRcount,length(which(drugR==i)))
    drugScount<-c(drugScount,length(which(drugS==i)))
  }
  #to include the Resistant and Susceptible counts in data frame
  S_count[which(names(S_count) == drugName)] <- drugScount
  R_count[which(names(R_count) == drugName)] <- drugRcount
  #to create csv files with S_count and R_count and is stored in Output folder   , "NumRes", "NumSusc"
  write.csv(S_count, "Output/S_count.csv")
  write.csv(R_count, "Output/R_count.csv")
}

#to change the format of the data sets 
R_Data = reshape2::melt(R_count, id.vars = c("UniqueStrains")) 
S_Data = reshape2::melt(S_count, id.vars = c("UniqueStrains")) 

#to merge resistance and susceptible counts
ST_Data = merge(R_Data, S_Data, by = c("UniqueStrains", "variable"), all.x=T)

#to change column name 
colnames(ST_Data)[1:4] <- c("ST", "Drug", "NumRes", "NumSusc")
ST_Data$Pathogen <- "S. aureus"

#to write into csv file (we can write only once otherwise it will keep modifying the data set in dropbox)
#write.csv(ST_Data, "Data/Wurster_ST_Data.csv") 




 

