#This Rscripts divides Kallonen dataset into two datasets: Kallonen_BSAC and Kallonen_CUH
#getwd()
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
#WST1 <-read.csv
#to read the .csv file
Kall = read.csv(file = "Data/PatientLevel_Data/Kallonen_E_coli_Patient_Data.csv", header=TRUE, sep = ",", check.names = F, stringsAsFactors = FALSE)


#install.packages("dplyr")
library(reshape2)
library(dplyr)

#to extract the cloumns we require; here we only require Sequence Type and Drug Names
#7th colum = "MLST" 
#20-42 = Drug Names

DrugNames<-  c("AmoxiClav", "Amoxicillin", "Ceftazidime", "Ciprofloxacin", "Cefotaxime", "Cefuroxime",
               "Gentamicin", "Imipenem", "Tigecycline", "PipTaz","Amikacin", "Ampicillin", "Aztreonam",
               "Cefalotin", "Cefepime", "Cefoxitin", "CefuroximeAxetil", "Ertapenem", "Meropenem",
               "Tobramycin", "Trimethoprim") 
  
#Kall = select(Kall, 3, 7, 20:42) -> it did not work so replaced with new line of code in line 27
#divided the anibiotics column since there are two extra column ("esbl", "esbl_ctxm") in between 
#to slect the columns with  information about name of organization, ST and antibiotic from original csv 
Kall_BSAC = Kall[Kall$Collection=="BSAC", c(3,7, 20:29)] #excluded all the columns which has no any information about antibiotics
Kall_CUH = Kall[Kall$Collection=="CUH", c(3, 7, 20, 22:26, 28:29, 32:42)] ##excluded all the columns which has no any information about antibiotics


#################################################################################################################################################################################################################
#to work on Kallonen datasets collected from BSAC
drugnames = names(Kall_BSAC)[3:12]#refers two the third through 25th column, this is where the names of drug starts and ends
#drugnames = c("Clindamycin", "Daptomycin", "Fusidic acid", "Gentamicin", "Levofloxacin", "Linezolid", "Rifampicin", "Teicoplanin", "Tigecyclin", "Trimethroprim_Sulfamethoxazole", "Vancomycin", "Oxacillin", "Penicillin G", "Tetracyclin")

#to get Sequence Type (ST) from .csv file
ST = Kall_BSAC$MLST

#UniqueCollection = sort(unique(Collection))
#make list of the Sequence Type (each ST is represented only once and sorted)
UniqeStrains = sort(unique(ST))
#UniqeStrains

#to create new dataframe for Simpson's Diversity Index (SI) for Kallonen-BSAC
S_count <- data.frame("UniqueStrains" = UniqeStrains, "1" = 0, "2" = 0, "3" = 0, "4" = 0, "5" = 0, "6" = 0, "7" = 0, "8" = 0, "9" = 0, "10" = 0)
names(S_count)[2:11] <- drugnames #to assign drug names to the column
R_count = S_count
#to check drug name corresponding to the column
#the for loop for all the drugs starts here, therefore make sure you indentate for everything till the first for loop closes
#this for loop fills above created data fame with "R" (resistane) and "S" (susceptible); NOTE:we excluded "I" (intermediate)
for (drugColumn in  2:12) {
  drugName = names(Kall_BSAC)[drugColumn]
  drugR = ST[which(Kall_BSAC[, drugColumn]=="R")]
  drugS = ST[which(Kall_BSAC[, drugColumn]=="S")]
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
#ST_Data$NumRes <- as.integer(as.character(ST_Data$NumRes))
#ST_Data$NumSusc <- as.integer(as.character(ST_Data$NumSusc))
#ST_Data$ST <- as.factor(ST_Data$ST)
ST_Data$Pathogen <- "E. coli"

#to write into csv file (we can write only once otherwise it will keep modifying the data set in dropbox)
write.csv(ST_Data, "Data/Kallonen_BSAC_ST_Data.csv") 


#################################################################################################################################################################################################################
#to work on Kallonen datasets collected from CUH
drugnames = names(Kall_CUH)[3:21]
#to get Sequence Type (ST) from .csv file
ST = Kall_CUH$MLST

#UniqueCollection = sort(unique(Collection))
#make list of the Sequence Type (each ST is represented only once and sorted)
UniqeStrains = sort(unique(ST))
#UniqeStrains

#to create new dataframe for Simpson's Diversity Index (SI) for Kallonen-CUH
S_count <- data.frame("UniqueStrains" = UniqeStrains, "1" = 0, "2" = 0, "3" = 0, "4" = 0, "5" = 0, "6" = 0, "7" = 0, "8" = 0, "9" = 0, "10" = 0, "11" = 0, "12" = 0, "13" = 0, "14" = 0, 
                      "15" = 0, "16" = 0, "17" = 0, "18" = 0, "19" = 0)
names(S_count)[2:20] <- drugnames #to assign drug names to the column
R_count = S_count
#to check drug name corresponding to the column
#the for loop for all the drugs starts here, therefore make sure you indentate for everything till the first for loop closes
#this for loop fills above created data fame with "R" (resistane) and "S" (susceptible); NOTE:we excluded "I" (intermediate)
for (drugColumn in  2:21) {
  drugName = names(Kall_CUH)[drugColumn]
  drugR = ST[which(Kall_CUH[, drugColumn]=="R")]
  drugS = ST[which(Kall_CUH[, drugColumn]=="S")]
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
ST_Data$NumRes <- as.integer(as.character(ST_Data$NumRes))
ST_Data$NumSusc <- as.integer(as.character(ST_Data$NumSusc))
ST_Data$ST <- as.factor(ST_Data$ST)
#as.integer(as.character(ST_Data$NumRes,ST_Data$NumSusc))
ST_Data$Pathogen <- "E. coli"

#to write into csv file (we can write only once otherwise it will keep modifying the data set in dropbox)
write.csv(ST_Data, "Data/Kallonen_CUH_ST_Data.csv") 




