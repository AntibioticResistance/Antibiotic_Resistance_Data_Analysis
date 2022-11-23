
#to set up working directory
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")

#to read the .csv file
Enterococcus = read.csv(file = "Data/PatientLevel_Data/Galloway_Patient_Data.csv", header=TRUE, sep = ",", fileEncoding="latin1")

#to check if the packages we need are present
#library(ggplot2)
library(reshape2)
library(dplyr)
#install.packages(cowplot)  
#library(colorspace)
#library(grDevices)
#install.packages(colormap)
#library(data.table)
#library(RColorBrewer)

#Galloway = name of the author
#to extract the cloumns we require; here we only require Sequence Type and Drug Names
#8th colum = Sequence Type
#4-7 = Drug Names
GALLOWAY = select(Enterococcus, 8, 4:7)

#to get Sequence Type(ST) from .csv file
ST = GALLOWAY$Sequence.Type

#make list of the Sequence Type (each ST is represented only once and sorted) and drug names
UniqeStrains = sort(unique(ST))
drugnames = names(GALLOWAY)[2:5]#refers 2-5 columns which represnts antibiotics used

#to create new dataframe for Simpson's Diversity Index (SI)
S_count <- data.frame("UniqueStrains" = UniqeStrains, "1" = 0, "2" = 0, "3" = 0, "4" = 0)
names(S_count)[2:5] <- drugnames #to assign drug names to the column
R_count = S_count

#the for loop for all the drugs starts here, therefore make sure you indentation for everything till the first for loop closes at the end 
for (drugColumn in  2:5) {
  drugName = names(GALLOWAY)[drugColumn]
  drugR = ST[which(GALLOWAY[, drugColumn]=="R")]
  drugS = ST[which(GALLOWAY[, drugColumn]=="S")]
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
}

#to change the format of the data sets 
R_Data = reshape2::melt(R_count, id.vars = c("UniqueStrains")) 
S_Data = reshape2::melt(S_count, id.vars = c("UniqueStrains")) 

#to merge resistance and susceptible counts
ST_Data = merge(R_Data, S_Data, by = c("UniqueStrains", "variable"), all.x=T)

#to change column name 
colnames(ST_Data)[1:4] <- c("ST", "Drug", "NumRes", "NumSusc")
ST_Data$Pathogen <- "E. faecium"

#to write into csv file (we can write only once otherwise it will keep modifying the data set in dropbox)
#write.csv(ST_Data, "Data/Galloway_ST_Data.csv") 





