# =================================================================================
# =================================================================================
# Script:"SEER2RColonCancer"
# Date: 2022-04-04
# Implemented by: Sam Palmer
# Description:
# This script is used to extract the data sets from the database.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
library("SEER2R")
setwd("./")
# =================================================================================
# =================================================================================
# Extracting the data from the database
# =================================================================================
# =================================================================================
# Read all data into R
AllData <- read.SeerStat(DICfileName = "../Data/Incidence-Full-1Y-18.dic", TXTfileName = "../Data/Incidence-Full-1Y-18.txt", UseVarLabelsInData = TRUE)
AllData$`Site_recode_ICDO3/WHO_2008`=trimws(AllData$`Site_recode_ICDO3/WHO_2008`)
# Extract the CML data
Cancer="Chronic Myeloid Leukemia" 
DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                        AllData$Race_recode_W_B_AI_API == "White" &
                        AllData$Sex == "Male and female",]
# Write the CML data to a file
write.csv(DataSubset,"CML_cancer_all_fields.csv")

