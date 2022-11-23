
getwd()
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
#WST1 <-read.csv
#to read the .csv file
Staph = read.csv(file = "Data/PatientLevel_Data/Manara_S_aureus_Patient_Data.csv", header=TRUE, sep = ",", check.names = F, stringsAsFactors = FALSE)
#added "check.names = F" in previous code since it gave error because of the heading names in columns
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
#MANARA = name of the author
#to extract the cloumns we require; here we only require Sequence Type and Drug Names
#21th colum = Sequence Type
#29-42 = Drug Names
MANARA = select(Staph, 21, 29:42)
MANARA = MANARA[1:184,] #this line removes all the extra "NA" values in the csv file; total number of isolates were 184, but there were NA values even in empty coulmns upto 950
colnames(MANARA)[2:15] <- c("Clindamycin", "Daptomycin", "Fusidic acid", "Gentamicin", "Levofloxacin", "Linezolid", "Rifampicin", "Teicoplanin", "Tigecyclin", "Trimethroprim_Sulfamethoxazole", "Vancomycin", "Oxacillin", "Penicillin G", "Tetracyclin")

#to get Sequence Type (ST) from .csv file
ST = MANARA$ST
#make list of the Sequence Type (each ST is represented only once and sorted)
UniqeStrains = sort(unique(ST))
UniqeStrains

drugnames = names(MANARA)[2:15]#refers two the second through 15th column, this is where the names of drug starts and ends
#drugnames = c("Clindamycin", "Daptomycin", "Fusidic acid", "Gentamicin", "Levofloxacin", "Linezolid", "Rifampicin", "Teicoplanin", "Tigecyclin", "Trimethroprim_Sulfamethoxazole", "Vancomycin", "Oxacillin", "Penicillin G", "Tetracyclin")
#to create new dataframe for Simpson's Diversity Index (SI)
S_count <- data.frame("UniqueStrains" = UniqeStrains, "1" = 0, "2" = 0, "3" = 0, "4" = 0, "5" = 0, "6" = 0, "7" = 0, "8" = 0, "9" = 0, "10" = 0, "11" = 0, "12" = 0, "13" = 0, "14" = 0)
names(S_count)[2:15] <- drugnames #to assign drug names to the column
R_count = S_count
#to check drug name corresponding to the column
#the for loop for all the drugs starts here, therefore make sure you indentate for everything till the first for loop closes
#this for loop fills above created data fame with "R" (resistane) and "S" (susceptible); NOTE:we excluded "I" (intermediate)
for (drugColumn in  2:15) {
  drugName = names(MANARA)[drugColumn]
  drugR = ST[which(MANARA[, drugColumn]=="R")]
  drugS = ST[which(MANARA[, drugColumn]=="S")]
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
ST_Data$Pathogen <- "S. aureus"

#to write into csv file (we can write only once otherwise it will keep modifying the data set in dropbox)
write.csv(ST_Data, "Data/Manara_ST_Data.csv") 






