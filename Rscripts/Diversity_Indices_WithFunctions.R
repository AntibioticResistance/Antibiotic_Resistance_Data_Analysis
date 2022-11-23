

#getwd()
#setwd("/Users/Candace Clark/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")
library(dplyr)

############ To calculate random Simpson's Index of Diversity (Difference in SID) ##############
####################### SIMULATION ##########################
#############################################################################

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

drugRfiles<-list.files("Data/", pattern="ST_Data.csv") #read the files with data
DivIndices<-data.frame("Dataset"= character(), "Pathogen" = character(), "Drug" = character(), "NumRes"= integer(), "NumSus"= integer(),  "SIvalueR" = double(), "SIvalueS" =  double(), "SIpvalue" = double() )

#for (g in 1:length(drugRfiles)){
for (g in 1:length(drugRfiles)){
    #g = 1
  print(drugRfiles[g])
  df<-read.csv(paste0("Data/",drugRfiles[g]), stringsAsFactors = F)#, row.names = 1
  #to change the name of antibiotics that are different in different datasets
  df$Drug[df$Drug == "Trimethoprim-Sulfamethoxazole"]<- "Trim-Sulf"
  df$Drug[df$Drug == "Trimethroprim_Sulfamethoxazole"]<- "Trim-Sulf"
  df$Drug[df$Drug == "AmoxiClav"]<- "Amoxi-Clav"
  df$Drug[df$Drug == "PipTaz"]<- "Pip-Taz"
  df$Drug[df$Drug == "CefuroximeAxetil"]<- "Cefuroxime axetil"
  df$Drug[df$Drug == "Tigecyclin"]<- "Tigecycline"
  df$Drug[df$Drug == "Tetracyclin"]<- "Tetracycline"
  df$Drug[df$Drug == "Penicillin G"]<- "Penicillin"
  
  #print(unique(df$ST))
  df<-df[which(df$ST!="Minor" & df$ST != "ND " & df$ST != "-"),] #Oct 2021 Pleuni removed the minor category. 
  #It is not one ST, but many, if we count it as one it messes up our diversity calculations). 
  dname<-gsub("_ST_Data.csv",'',drugRfiles[g])
  DrugList<-sort(unique(df$Drug))
  
  sum_numRes<- df %>% group_by(Drug) %>% summarize(sum = sum(NumRes))
  DrugList<-DrugList[which(sum_numRes$sum>0)]
  
  pdf(file = paste0("Output/DivIndices_Histograms/Simulated_Histogram_GSI_", dname,".pdf"), width = 8, height = 8) #width = 8. height = 2.5
  par(mfrow=c(2,3), oma=c(0,0,2,0))
  
  for (i in 1:length(DrugList)){ 
    #  i=1
    #stepwise process to calculate observed SI value 
    drugName<-DrugList[i]
    #subset the data frame to a specific drug
    df_for_drug<-df[df$Drug==drugName,]
    
    #Calculate SI values of Resistant population
    GSIvalueRealR = CalcGiniSimpsonDivIndex(df_for_drug$NumRes)
    InvSIvalueRealR = CalcInverseSimpsonIndex(df_for_drug$NumRes)
    ShanvalueRealR = CalcShannonIndex(df_for_drug$NumRes)
    
    #Calculate SI values of Susceptible population
    GSIvalueRealS = CalcGiniSimpsonDivIndex(df_for_drug$NumSusc)
    InvSIvalueRealS = CalcInverseSimpsonIndex(df_for_drug$NumSusc)
    ShanvalueRealS = CalcShannonIndex(df_for_drug$NumSusc)
    
    #create the data frame to perform random sampling
    re.ST<- rep(df_for_drug$ST, times=df_for_drug$NumRes)
    su.ST<-rep(df_for_drug$ST, times=df_for_drug$NumSusc)
    
    GSIvaluesRES = c(); GSIvaluesSUS = c() 
    InvSIvaluesRES = c(); InvSIvaluesSUS = c()
    ShannonvaluesRES = c(); ShannonvaluesSUS = c()
    
    #1000 simulation of random data set
    for(j in 1:1000)  {#added seed for reproducible data 
      set.seed(j)
      #Resistant population
      drugR <-sample(re.ST, length(re.ST),replace = T) #collecting only randomized data from the oringinal data set
      SimResultsR<-data.frame(table(drugR))
      GSIvaluesRES = c(GSIvaluesRES, CalcGiniSimpsonDivIndex(SimResultsR$Freq))
      InvSIvaluesRES = c(InvSIvaluesRES, CalcInverseSimpsonIndex(SimResultsR$Freq))
      ShannonvaluesRES = c(ShannonvaluesRES, CalcShannonIndex(SimResultsR$Freq))
      #Susceptible population
      drugS <-sample(su.ST, length(su.ST),replace = T) #collecting only randomized data from the oringinal data set
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
    pValueGSI = (length(which(differenceGSISim<0)))/1000 
    pValueInvSI = (length(which(differenceInvSISim<0)))/1000 
    pValueShannon = (length(which(differenceShannonSim<0)))/1000 
    
    #to make histogram
    hist(differenceGSISim, xlim = range(c(differenceGSISim, -0.3,0.3)), ylim=c(0,300),breaks=20,
         main = paste(drugName, "Gini-Simp Ind"),
         xlab = "GSI values")
    mtext("Bootstrapped difference in diversity Sus - Res", line=-0.1, side=3, outer=TRUE, cex=0.9)
    mtext(paste0(df$Pathogen," , ", dname) ,side = 3, line = -0.999, outer=TRUE, cex=0.8)
    lines(x = c(0,0), y = c(0,300), col = "blue")
    text(0, 300, labels = paste0("p-value = ", pValueGSI),col = " dark green")
    
    hist(differenceInvSISim, xlim = range(c(differenceInvSISim, -2,2)), 
         ylim = c(0,300), main = paste(drugName, "InverseSI"))
    lines(x = c(0,0), y = c(0,300), col = "blue")
    text(0, 300, labels = paste0("p-value = ", pValueInvSI),col = " dark green")
    
    hist(differenceShannonSim, xlim = range(c(differenceShannonSim, -2,2)), 
         ylim = c(0,300), main = paste(drugName, "Shannon"))
    lines(x = c(0,0), y = c(0,300), col = "blue")
    text(0, 300, labels = paste0("p-value = ", pValueShannon),col = " dark green")
    
    #add info to DivIndices
    DivIndices<-rbind(DivIndices, 
                        data.frame("Dataset"=dname,"Pathogen" = df$Pathogen[1], 
                                   "Drug" = drugName, "NumRes"= sum(df_for_drug$NumRes), 
                                   "NumSus"= sum(df_for_drug$NumSusc), 
                                   "GSIvalueR" = GSIvalueRealR, "GSIvalueS" = GSIvalueRealS, "GSIpvalue" = pValueGSI, 
                                   "InvSIvalueR" = InvSIvalueRealR, "InvSIvalueS" = InvSIvalueRealS, "InvSIpvalue" = pValueInvSI,
                                   "ShannonvalueR" = ShanvalueRealR, "ShannonvalueS" = ShanvalueRealS, "Shannonpvalue" = pValueShannon
                        ))
  }
  dev.off()
}

#to create a csv file with all the results created above and date in file name, and stores in "Old" folder 
#file.copy(from=paste("Output/DivIndices",".csv", sep = ""), to=paste("Output/Old/Archived_DivIndices",Sys.Date(),".csv", sep = ""))

#to create a csv file with all the results created above, and stores in "Output" folder
#03/10/2022: Anjani added (sep = "") to remove space in the file name
write.csv(x = DivIndices, file = paste("Output/DivIndices",".csv", sep = ""), row.names = FALSE)

