
#######################
####
## Prepare to make calculate diversity indeces 
####
#######################

setwd("~/Documents/GitHub/Antibiotic_Resistance_Data_Analysis")
library(dplyr)

#Write function to calculate SID Simpson Div Index
CalcGiniSimpsonDivIndex<-function(array_of_counts){
  array_of_fractions = array_of_counts/sum(array_of_counts)
  squared_array_of_fractions = array_of_fractions^2
  sum_squared_array_of_fractions = sum(squared_array_of_fractions)
  SIvalue = 1-sum_squared_array_of_fractions
  return(SIvalue)
}
CalcInverseSimpsonIndex<-function(array_of_counts){
  array_of_fractions = array_of_counts/sum(array_of_counts)
  squared_array_of_fractions = array_of_fractions^2
  sum_squared_array_of_fractions = sum(squared_array_of_fractions)
  SIvalue = 1/sum_squared_array_of_fractions
  return(SIvalue)
}
CalcShannonIndex <- function(array_of_counts){
  array_of_fractions = array_of_counts/sum(array_of_counts)
  log_array_of_fractions = array_of_fractions * log(array_of_fractions, base = 2)
  sum_log_array_of_fractions = sum(log_array_of_fractions, na.rm=T)
  Shannonvalue = -sum_log_array_of_fractions
  return(Shannonvalue)
}

#create list of all csv files we will read
ResDocs<-list.files("Data/", pattern ="_ST_Data.csv")
#Create empty dataframe to store data
DivIndices<-data.frame("Dataset"= character(), "Pathogen" = character(), "Drug" = character(), "NumRes"= integer(), "NumSus"= integer(),  "SIvalueR" = double(), "SIvalueS" =  double(), "SIpvalue" = double() )
#I want to remove a bunch of ST that are not really sequence types
ListToRemove <- c("-", "Minor", "ND " ,  "OC39" , "U1"  ,  "U10"  , "U11",   "U12",   "U13" ,  "U14" ,  "U15" ,  "U2"  , 
                  "U3"  ,  "U4" ,   "U5"  ,  "U6"  ,  "U7" ,   "U8"  ,  "U9")

#######################
####
## Loop through each dataset, and then for each dataset loop through all drugs 
####
#######################

for (g in 1:length(ResDocs)){
  print(drugRfiles[g])
  DataSet<-read.csv(paste0("Data/",drugRfiles[g]), stringsAsFactors = F)#, row.names = 1
  ## List of ST to remove (because unclear if they are ST)
  DataSet<-DataSet[which(!DataSet$ST %in% ListToRemove),] #May 2024 Pleuni removed bunch of ST categories. 
  #to change the name of antibiotics that are different in different datasets
  DataSet$Drug[DataSet$Drug == "Trimethoprim-Sulfamethoxazole"]<- "Trim-Sulf"
  DataSet$Drug[DataSet$Drug == "Trimethroprim_Sulfamethoxazole"]<- "Trim-Sulf"
  DataSet$Drug[DataSet$Drug == "AmoxiClav"]<- "Amoxi-Clav"
  DataSet$Drug[DataSet$Drug == "PipTaz"]<- "Pip-Taz"
  DataSet$Drug[DataSet$Drug == "CefuroximeAxetil"]<- "Cefuroxime axetil"
  DataSet$Drug[DataSet$Drug == "Tigecyclin"]<- "Tigecycline"
  DataSet$Drug[DataSet$Drug == "Tetracyclin"]<- "Tetracycline"
  DataSet$Drug[DataSet$Drug == "Penicillin G"]<- "Penicillin"
  
  dname<-gsub("_ST_Data.csv",'',drugRfiles[g]) #Get the study name
  DrugList<-sort(unique(DataSet$Drug)) #get the list of drugs used in this study
  
  sum_numRes<- DataSet %>% group_by(Drug) %>% summarize(sum = sum(NumRes))
  DrugList<-DrugList[which(sum_numRes$sum>0)] #Keep only the drugs for which there is some resistance
  
  pdf(file = paste0("Output/DivIndices_Histograms/Simulated_Histogram_GSI_", dname,".pdf"), width = 8, height = 8) #width = 8. height = 2.5
  par(mfrow=c(2,3), oma=c(0,0,2,0))
  
  for (i in 1:length(DrugList)){ 
    #stepwise process to calculate observed SI value 
    drugName<-DrugList[i]
    #subset the data frame to a specific drug
    df_for_drug<-DataSet[DataSet$Drug==drugName,]
    
    #Calculate SI values of Resistant population
    GSIvalueRealR = CalcGiniSimpsonDivIndex(df_for_drug$NumRes)
    InvSIvalueRealR = CalcInverseSimpsonIndex(df_for_drug$NumRes)
    ShanvalueRealR = CalcShannonIndex(df_for_drug$NumRes)
    
    #Calculate SI values of Susceptible population
    GSIvalueRealS = CalcGiniSimpsonDivIndex(df_for_drug$NumSusc)
    InvSIvalueRealS = CalcInverseSimpsonIndex(df_for_drug$NumSusc)
    ShanvalueRealS = CalcShannonIndex(df_for_drug$NumSusc)
    
    #create a list of all samples to perform random sampling from
    #this list lists one by one the ST for the resistant samples and the ST for the susc samples. 
    ListForSampling <- rep(df_for_drug$ST, times=(df_for_drug$NumRes + df_for_drug$NumSusc))

    GSIvaluesRES = c(); GSIvaluesSUS = c() 
    InvSIvaluesRES = c(); InvSIvaluesSUS = c()
    ShannonvaluesRES = c(); ShannonvaluesSUS = c()
    
    #1000 simulation of random data set
    for(j in 1:1000)  {#added seed for reproducible data 
      set.seed(j)
      #Resistant population
      drugR <-sample(ListForSampling, sum(df_for_drug$NumRes),replace = T) #collecting only randomized data from the oringinal data set
      SimResultsR<-data.frame(table(drugR))
      GSIvaluesRES = c(GSIvaluesRES, CalcGiniSimpsonDivIndex(SimResultsR$Freq))
      InvSIvaluesRES = c(InvSIvaluesRES, CalcInverseSimpsonIndex(SimResultsR$Freq))
      ShannonvaluesRES = c(ShannonvaluesRES, CalcShannonIndex(SimResultsR$Freq))
      #Susceptible population
      drugS <-sample(ListForSampling, sum(df_for_drug$NumSusc),replace = T) #collecting only randomized data from the oringinal data set
      SimResultsS<-data.frame(table(drugS))
      GSIvaluesSUS = c(GSIvaluesSUS, CalcGiniSimpsonDivIndex(SimResultsS$Freq))
      InvSIvaluesSUS = c(InvSIvaluesSUS, CalcInverseSimpsonIndex(SimResultsS$Freq))
      ShannonvaluesSUS = c(ShannonvaluesSUS, CalcShannonIndex(SimResultsS$Freq))
    }
    
    differenceGSIReal = GSIvalueRealS-GSIvalueRealR #SI values from data set 
    differenceGSISim = GSIvaluesSUS-GSIvaluesRES #SI values from 1000 simulation
    differenceInvSIReal = InvSIvalueRealS-InvSIvalueRealR #SI values from data set 
    differenceInvSISim = InvSIvaluesSUS-InvSIvaluesRES #SI values from 1000 simulation
    differenceShannonReal = ShanvalueRealS - ShanvalueRealR #real shannon difference
    differenceShannonSim = ShannonvaluesSUS - ShannonvaluesRES #Shannon from 1000 sims
    
    #to calculate P-value
    pValueGSI = (length(which(differenceGSISim>=differenceGSIReal)))/1000 
    pValueInvSI = (length(which(differenceInvSISim>=differenceInvSIReal)))/1000 
    pValueShannon = (length(which(differenceShannonSim>=differenceShannonReal)))/1000 
    
    #to make histogram
    hist(differenceGSISim, xlim = range(c(differenceGSIReal,differenceGSISim, -0.3,0.3)), ylim=c(0,300),breaks=20,
         main = paste(drugName, "Gini-Simp Ind"),
         xlab = "GSI values")
    mtext("Bootstrapped difference in diversity Sus - Res", line=-0.1, side=3, outer=TRUE, cex=0.9)
    mtext(paste0(DataSet$Pathogen," , ", dname) ,side = 3, line = -0.999, outer=TRUE, cex=0.8)
    lines(x = c(differenceGSIReal,differenceGSIReal), y = c(0,350), col = "blue")
    text(0, 300, labels = paste0("p-value = ", pValueGSI),col = " dark green")
    
    hist(differenceInvSISim, xlim = range(c(differenceInvSIReal,differenceInvSISim, -2,2)), 
         ylim = c(0,350), main = paste(drugName, "InverseSI"))
    lines(x = c(differenceInvSIReal,differenceInvSIReal), y = c(0,300), col = "blue")
    text(0, 300, labels = paste0("p-value = ", pValueInvSI),col = " dark green")
    
    hist(differenceShannonSim, xlim = range(c(differenceShannonReal,differenceShannonSim, -2,2)), 
         ylim = c(0,350), main = paste(drugName, "Shannon"))
    lines(x = c(differenceShannonReal,differenceShannonReal), y = c(0,300), col = "blue")
    text(0, 300, labels = paste0("p-value = ", pValueShannon),col = " dark green")
    
    #add info to DivIndices
    DivIndices<-rbind(DivIndices, 
                        data.frame("Dataset"=dname,"Pathogen" = DataSet$Pathogen[1], 
                                   "Drug" = drugName, "NumRes"= sum(df_for_drug$NumRes), 
                                   "NumSus"= sum(df_for_drug$NumSusc), 
                                   "GSIvalueR" = GSIvalueRealR, "GSIvalueS" = GSIvalueRealS, "GSIpvalue" = pValueGSI, 
                                   "InvSIvalueR" = InvSIvalueRealR, "InvSIvalueS" = InvSIvalueRealS, "InvSIpvalue" = pValueInvSI,
                                   "ShannonvalueR" = ShanvalueRealR, "ShannonvalueS" = ShanvalueRealS, "Shannonpvalue" = pValueShannon
                        ))
  }
  dev.off()
}

#to create a csv file with all the results created above, and stores in "Output" folder
write.csv(x = DivIndices, file = paste("Output/DivIndices",".csv", sep = ""), row.names = FALSE)

