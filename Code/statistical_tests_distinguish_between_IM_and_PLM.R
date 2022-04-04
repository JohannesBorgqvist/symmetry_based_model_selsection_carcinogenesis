# =================================================================================
# =================================================================================
# Script:"statistical_tests_distinguish_between_models"
# Date: 2022-04-04
# Implemented by: Sam Palmer and Johannes Borgqvist
# Description:
# The purpose of this script is to test on which data set that the two candidate models are distinguishable or indistinguishable. The three data sets are corresponds to the increase of incidences of cancer as a function of age in the case of myeloma cancer, colon cancer and chronic myeloid leukemia (CML) respectively. The two candidate models are the power law model (PLM) and the immunological model (IM). To distinguish between these two candidate models the the variance test and Vuong's non-nested likelihood ratio test are conducted. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
#For the non-nested variance test
if (!require('nonnest2')) install.packages('nonnest2'); library(nonnest2) 
# Analysis of Dose-Response Curves enabling flexible and versatile model fitting and after-fitting functions.
if (!require('drc')) install.packages('drc'); library(drc)
citation('drc')
# =================================================================================
# =================================================================================
# Reading the data and extracting the time series
# =================================================================================
# =================================================================================
# Set the working directory to the current directory
setwd("./")
# Read the myeloma cancer data
myeloma_datacsv = read.csv("../Data/Myeloma_cancer.csv")
myeloma_rates = myeloma_datacsv[,"AgeAdjusted_Rate"][1:86] 
# Myeloma time series: Define the actual time series as a structure. We only consider cancer incidences from the age of 25 years and above
myeloma_data <- structure(list(R = myeloma_rates[25:86], 
                       t = 25:86), 
                  .Names = c("R", "t"), 
                  class = "data.frame", 
                  row.names = c(NA, -86L))
# Read the colon cancer data
colon_datacsv = read.csv("../Data/Colon_cancer.csv")
colon_rates = colon_datacsv[,"AgeAdjusted_Rate"][1:86] 
# Colon time series: Define the actual time series as a structure. We only consider cancer incidences from the age of 12 years and above.
colon_data <- structure(list(R = colon_rates[12:86], 
                       t = 12:86), 
                  .Names = c("R", "t"), 
                  class = "data.frame", 
                  row.names = c(NA, -86L))
# Read the CML data
CML_datacsv = read.csv("../Data/CML_cancer.csv")
CML_rates = CML_datacsv[,"AgeAdjusted_Rate"][1:86] 
# CML time series: Define the actual time series as a structure. We only consider cancer incidences from the age of 10 years and above.
CML_data <- structure(list(R = CML_rates[10:86], 
                       t = 10:86), 
                  .Names = c("R", "t"), 
                  class = "data.frame", 
                  row.names = c(NA, -86L))
# =================================================================================
# =================================================================================
# Define the candidate models
# =================================================================================
# =================================================================================
# MYELOMA CANCER DATA
myeloma.model.IM <- nls(R~A/(exp(exp(-0.044*(t-tau)))-C), data=myeloma_data, start = list(A=100,tau=60.5,C=1)) # IM
myeloma.model.PLM <- nls(R~A*t^B, data=myeloma_data, start = list(A=max(myeloma_data$R)/2, B=4)) # PLM
# COLON CANCER DATA
colon.model.IM <- nls(R~A/(exp(exp(-0.044*(t-tau)))-C), data=colon_data, start = list(A=200,tau=60.5,C=1)) # IM
colon.model.PLM <- nls(R~A*t^B, data=colon_data, start = list(A=max(colon_data$R)/2, B=4)) # PLM
# CML CANCER DATA
CML.model.IM <- nls(R~A/(exp(exp(-0.044*(t-tau)))-C), data=CML_data, start = list(A=1,tau=60.5,C=1)) # IM
CML.model.PLM <- nls(R~A*t^B, data=CML_data, start = list(A=max(CML_data$R)/2, B=4)) # PLM
# =================================================================================
# =================================================================================
# Calculate the summary statistics of all models
# =================================================================================
# =================================================================================
# MYELOMA CANCER DATA
summary(myeloma.model.IM) # IM
summary(myeloma.model.PLM) # PLM
# COLON CANCER DATA
summary(colon.model.IM) # IM
summary(colon.model.PLM) # PLM
# CML CANCER DATA
summary(CML.model.IM) # IM
summary(CML.model.PLM) # PLM
# =================================================================================
# =================================================================================
# Conducts Vuong's nested likelihood ratio test on the candidate models
# =================================================================================
# =================================================================================
# MYELOMA CANCER DATA
vuongtest(myeloma.model.IM,myeloma.model.PLM)
# COLON CANCER DATA
vuongtest(colon.model.IM,colon.model.PLM)
# CML CANCER DATA
vuongtest(CML.model.IM,CML.model.PLM)
# =================================================================================
# =================================================================================
# Conducts Vuong's nested likelihood ratio test on the candidate models
# =================================================================================
# =================================================================================
citation()
citation('drc')
citation('nonnest2')
