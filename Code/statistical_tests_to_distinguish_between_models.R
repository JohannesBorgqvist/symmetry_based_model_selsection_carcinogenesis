# =================================================================================
# =================================================================================
# Script:"statistical_tests_distinguish_between_models"
# Date: 2022-01-04
# Implemented by: Sam Palmer and Johannes Borgqvist
# Description:
# The purpose of this script is to test on which data set that the two candidate models are distinguishable or indistinguishable. The two data sets are corresponds to the increase of incidences of cancer as a function of age in the case of myeloma cancer and colon cancer respectively. The two candidate models are the power law model (PLM) and the immunological model (IM-II). To distinguish between these two candidate models the the variance test and Vuong's non-nested likelihood ratio test are conducted. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
library(nonnest2) # For the non-nested variance test
library(drc) # Analysis of Dose-Response Curves enabling flexible and versatile model fitting and after-fitting functions.
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
# =================================================================================
# =================================================================================
# Define the candidate models
# =================================================================================
# =================================================================================
# MYELOMA CANCER DATA
myeloma.model.IMII <- nls(R~A/(exp(exp(-0.044*(t-tau)))-1), data=myeloma_data, start = list(A=2,tau=60.5)) # IM-II
myeloma.model.PLM <- nls(R~A*t^B, data=myeloma_data, start = list(A=max(myeloma_data$R)/2, B=4)) # PLM
# COLON CANCER DATA
colon.model.IMII <- nls(R~A/(exp(exp(-0.044*(t-tau)))-1), data=colon_data, start = list(A=2,tau=60.5)) # IM-II
colon.model.PLM <- nls(R~A*t^B, data=colon_data, start = list(A=max(colon_data$R)/2, B=4)) # PLM
# =================================================================================
# =================================================================================
# Calculate the summary statistics of all models
# =================================================================================
# =================================================================================
# MYELOMA CANCER DATA
summary(myeloma.model.IMII) # IMII
summary(myeloma.model.PLM) # PLM
# COLON CANCER DATA
summary(colon.model.IMII) # IMII
summary(colon.model.PLM) # PLM
# =================================================================================
# =================================================================================
# Conducts Vuong's nested likelihood ratio test on the candidate models
# =================================================================================
# =================================================================================
# MYELOMA CANCER DATA
vuongtest(myeloma.model.IMII,myeloma.model.PLM)
# COLON CANCER DATA
vuongtest(colon.model.IMII,colon.model.PLM)

