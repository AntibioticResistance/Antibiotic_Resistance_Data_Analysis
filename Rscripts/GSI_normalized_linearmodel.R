##With this R script, I will normalize Gini Simpson Index by dividing by the theoretical max value of GSI 
##Use GLM to see effect of fraction resistant 
##Aug 2022
##

#############################################################################

library(ggplot2)

setwd("~/Documents/GitHub/Antibiotic_Resistance_Data_Analysis")
Data<-read.csv("Output/DivIndices.csv")
Data$Dataset<-factor(Data$Dataset, levels = 
                       c("Yamaji_1999", "Yamaji_2016", "Addams-Sapper", "Kallonen_BSAC", "Kallonen_CUH", "Wurster" , "Manara", "Galloway" ))


#I need to add "FracRes"       "TheorMaxGSI"    "GSINormalizedR"
Data$FracRes <- Data$NumRes/(Data$NumRes+Data$NumSus)
Data$TheorMaxGSI <--1 #placeholder, filled up below. 
for (j in 1:nrow(Data)){
  n = Data$NumRes[j]
  Data$TheorMaxGSI[j]<-1-sum(rep(1/n,n)^2)
} #calculating max SI value
Data$GSINormalizedR<-Data$GSIvalueR/Data$TheorMaxGSI

modglm<-glm(GSINormalizedR ~ FracRes + Dataset,  data = Data) #ADD Drug class to the model 
summary(modglm) 

sink(file = "Output/Figure3_lm_output_GSINormalizedR.txt")
summary(modglm)
sink(file = NULL)

#Make a dataframe with the coefficients from the glm model
dummy2 <- data.frame(Dataset = levels(Data$Dataset), Z = c(0,modglm$coefficients[3:9]))
dummy2$Dataset <- factor(dummy2$Dataset)

png("Output/Figure3_GSINormalized_values_Res.png",  width=150, height=150, units = "mm", res = 300)
g <- ggplot(Data, 
            aes(x = FracRes, y = GSIvalueR/TheorMaxGSI, color = Dataset)) #Drug should be DrugClass here
g<- g + 
  #labs(title = "Normalized Normalized Gini-Simpson Index vs Fraction Resistant Samples")+
  xlab("Fraction resistant samples") + ylab("Normalized Gini-Simpson Index")+
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

g +  geom_abline(data = dummy2, aes(intercept = Z + modglm$coefficients[1], slope = modglm$coefficients[2]), color= "darkgrey")

dev.off()

write.csv(x = Data, file = "Output/DataNormelizedGSI.csv")


