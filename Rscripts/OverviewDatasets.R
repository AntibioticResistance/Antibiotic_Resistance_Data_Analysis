#getwd()
setwd("~/Dropbox/2017_CodeLab/BacterialResistanceProject_2020/2020_CandaceUTIpopgen/DataAnalysis")

#install.packages("dplyr")
library(dplyr)

#to read csv file
#to read csv with overall datasets with all three Indices
Diversity_Indices <- read.csv(file = "Output/DivIndices.csv")

#to read csv with antibiotics and its classification
AntibioticClass <- read.csv(file = "Data/Antibiotic_Classification.csv")


#to combine Diversity Indices and Antibiotic class
OverviewDatasets <- merge(Diversity_Indices, AntibioticClass, by="Drug",all.x=TRUE) 

#to create a csv file with all the results created above and date in file name, and stores in "Old" folder 
#file.copy(from=paste("Output/OverviewDatasets",".csv", sep = ""), to=paste("Output/Old/Archived_OverviewDatasets",Sys.Date(),".csv", sep = ""))

#to create csv file (without date and stores in Output folder)
write.csv(x = OverviewDatasets,file = paste("Output/OverviewDatasets",".csv", sep = ""), row.names = FALSE)






