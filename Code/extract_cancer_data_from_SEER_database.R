# =================================================================================
# =================================================================================
# Script:"SEER2RColonCancer"
# Date: 2022-04-04
# Implemented by: Sam Palmer
# Description:
# The purpose of this script is to test on which data set that the two candidate models are distinguishable or indistinguishable. The three data sets are corresponds to the increase of incidences of cancer as a function of age in the case of myeloma cancer, colon cancer and chronic myeloid leukemia (CML) respectively. The two candidate models are the power law model (PLM) and the immunological model (IM). To distinguish between these two candidate models the the variance test and Vuong's non-nested likelihood ratio test are conducted. 
library(cowplot)

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

setwd("./")


library("SEER2R")



AllData <- read.SeerStat(DICfileName = "../Data/Incidence-Full-1Y-18.dic", TXTfileName = "../Data/Incidence-Full-1Y-18.txt", UseVarLabelsInData = TRUE)

AllData$`Site_recode_ICDO3/WHO_2008`=trimws(AllData$`Site_recode_ICDO3/WHO_2008`)


Cancer="Chronic Myeloid Leukemia" 
DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                        AllData$Race_recode_W_B_AI_API == "White" &
                        AllData$Sex == "Male and female",]

write.csv(DataSubset,"CML_cancer_all_fields.csv")

