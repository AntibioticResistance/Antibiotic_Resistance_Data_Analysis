#Comparing normalized vs unnormalized simulated data & effect of fract resistant


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

#create list of all csv files we will read
ResDocs<-list.files("Data/", pattern ="_ST_Data.csv")
#Create empty dataframe to store data
SimulatedData <- data.frame("Dataset"= character(), "Drug" = character(), "NumRes"= integer(), "FractionRes" = double(), "GSI_RESSim" = double(), "GSI_RESSimNormalized" = double())
#I want to remove a bunch of ST that are not really sequence types
ListToRemove <- c("-", "Minor", "ND " ,  "OC39" , "U1"  ,  "U10"  , "U11",   "U12",   "U13" ,  "U14" ,  "U15" ,  "U2"  , 
                  "U3"  ,  "U4" ,   "U5"  ,  "U6"  ,  "U7" ,   "U8"  ,  "U9")

#######################
####
## Loop through each dataset, and then for each dataset loop through all drugs 
####
#######################

for (g in 1:length(ResDocs)){
  print(ResDocs[g])
  DataSet<-read.csv(paste0("Data/",ResDocs[g]), stringsAsFactors = F)#, row.names = 1
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
  
  dname<-gsub("_ST_Data.csv",'',ResDocs[g]) #Get the study name
  DrugList<-sort(unique(DataSet$Drug)) #get the list of drugs used in this study
  
  sum_numRes<- DataSet %>% group_by(Drug) %>% summarize(sum = sum(NumRes))
  DrugList<-DrugList[which(sum_numRes$sum>0)] #Keep only the drugs for which there is some resistance
  
  #pdf(file = paste0("Output/DivIndices_Histograms/Simulated_Histogram_GSI_", dname,".pdf"), width = 8, height = 8) #width = 8. height = 2.5
  par(mfrow=c(1,1), oma=c(1,1,1,1))
  
  for (i in 1:length(DrugList)){ 
    drugName<-DrugList[i]
    #subset the data frame to a specific drug
    df_for_drug<-DataSet[DataSet$Drug==drugName,]
    #create a list of all samples to perform random sampling from
    #this list lists one by one the ST for the resistant samples and the ST for the susc samples. 
    ListForSampling <- rep(df_for_drug$ST, times=(df_for_drug$NumRes + df_for_drug$NumSusc))
    
    #1000 simulation of random data set
    numres = sum_numRes$sum[sum_numRes$Drug==drugName]
    if (numres>1){
      set.seed(numres)
      #Resistant population
      drugR <-sample(ListForSampling, numres,replace = T) #collecting only randomized data from the oringinal data set
      SimResultsR<-data.frame(table(drugR))
      GSI_RESSim = CalcGiniSimpsonDivIndex(SimResultsR$Freq)
      GSI_RESSimNormalized = GSI_RESSim / (1-sum(rep(1/numres,numres)^2))
      SimulatedData<-rbind(SimulatedData, 
                           data.frame("Dataset"= ResDocs[g], "Drug" = drugName, "NumRes"= numres, 
                                      "FractionRes" = numres/length(ListForSampling), "GSI_RESSim" = GSI_RESSim, 
                                      "GSI_RESSimNormalized" = GSI_RESSimNormalized ))
    }
  }
}

NormalizedSimDatalm<-lm(GSI_RESSimNormalized ~ FractionRes + Dataset,  data = SimulatedData) 
summary(NormalizedSimDatalm) 

UnNormalizedSimDatalm<-lm(GSI_RESSim ~ FractionRes + Dataset,  data = SimulatedData) 
summary(UnNormalizedSimDatalm)


sink(file = "Output/Simulated_NormalizedVsUnnormalized_output.txt")
summary(NormalizedSimDatalm)
summary(UnNormalizedSimDatalm)
sink(file = NULL)

write.csv(x = SimulatedData, file = "Output/SimulatedDataNormelizedGSI.csv")

#Make a dataframe with the coefficients from the lm model
dummyUnNormalizedSimDatalm <- data.frame(Dataset = levels(as.factor(SimulatedData$Dataset)), Z = c(0,UnNormalizedSimDatalm$coefficients[3:9]))
dummyUnNormalizedSimDatalm$Dataset <- factor(dummyUnNormalizedSimDatalm$Dataset)
dummyNormalizedSimDatalm <- data.frame(Dataset = levels(as.factor(SimulatedData$Dataset)), Z = c(0,NormalizedSimDatalm$coefficients[3:9]))
dummyNormalizedSimDatalm$Dataset <- factor(dummyNormalizedSimDatalm$Dataset)


pdf("Output/GSINormalizedVsUnNormalizedSimData_Res.pdf",  width=8, height=10)
g <- ggplot(SimulatedData, 
            aes(x = FractionRes, y = GSI_RESSim, color = Dataset)) #Drug should be DrugClass here
g<- g + 
  #labs(title = "Normalized Normalized Gini-Simpson Index vs Fraction Resistant Samples")+
  xlab("Fraction resistant samples") + ylab("Simulated Un_Normalized Gini-Simpson Index")+
  ylim(0.4,1)+
  geom_point()+
  facet_wrap(~Dataset)+
  theme_bw() + 
  theme(panel.border = element_blank(),
        text = element_text(size=11),
        #panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.3, vjust=0.5, size = 7),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
g <- g +   theme(legend.position="none")

g +  geom_abline(data = dummyUnNormalizedSimDatalm, aes(intercept = Z + UnNormalizedSimDatalm$coefficients[1], slope = UnNormalizedSimDatalm$coefficients[2]), color= "darkgrey")

##Normalized simulated datasets

g <- ggplot(SimulatedData, 
            aes(x = FractionRes, y = GSI_RESSimNormalized, color = Dataset)) #Drug should be DrugClass here
g<- g + 
  #labs(title = "Normalized Normalized Gini-Simpson Index vs Fraction Resistant Samples")+
  xlab("Fraction resistant samples") + ylab("Simulated Normalized Gini-Simpson Index")+
  ylim(0.4,1)+
  geom_point()+
  facet_wrap(~Dataset)+
  theme_bw() + 
  theme(panel.border = element_blank(),
        text = element_text(size=11),
        #panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.3, vjust=0.5, size = 7),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
g <- g +   theme(legend.position="none")

g +  geom_abline(data = dummyNormalizedSimDatalm, aes(intercept = Z + NormalizedSimDatalm$coefficients[1], slope = NormalizedSimDatalm$coefficients[2]), color= "darkgrey")

dev.off()
