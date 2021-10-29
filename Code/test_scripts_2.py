# =================================================================================
# =================================================================================
# Script:"test_scripts"
# Date: 2021-10-01
# Implemented by: Johannes Borgqvist
# Description:
# The script is made to just test the various scripts and functionalities
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import read_data  # Home-made
import write_output  # Home-made
import fit_to_data # Home-made
import symmetry_based_model_selection # Home-made
import numpy as np
import math
from matplotlib import pyplot as plt
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# Function 1: ""

# =================================================================================
# =================================================================================
# Test the scripts
# =================================================================================
# =================================================================================
#---------------------------------------------------------------------------------
# READ DATA
# Define the name of the file from which we will read the data
file_name = "Colon_cancer"
# Read the data
xlabel_str, ylabel_str, t, R = read_data.time_series_from_csv(file_name)
# Remove first 12 data points
t = t[12:len(t)-1]
R = R[12:len(R)-1]
#t = t[18:len(t)-1]
#R = R[18:len(R)-1]
# Logarithm of the output
for i in range(len(R)):
    R[i] = np.log(R[i])



#---------------------------------------------------------------------------------
# PLOT DATA
# Define the plot
R_plot = [np.exp(R_val) for R_val in R]
#---------------------------------------------------------------------------------
# DEFINE DATA STRING
data_str = "colon"
# FIT EXPONENTIAL MODEL
# Define parameter alpha common to most models
alpha = 0.044
# Save this parameter in a parameter vector
para_exp = [alpha]
para_mixed = [alpha]
para_pow = [alpha]
# Define the model string entailing that we are concentrating on the exponential
# model
model_str = "exponential"
# Set lower and upper bound on our parameter A
# Fit the model to the data
#print("Model 1")
R_hat_exp, opt_para_exp, SS_exp = fit_to_data.PE_risk_profiles(t,R,model_str,para_exp)
# Save the estimated parameters
write_output.save_data_PE(data_str,model_str,opt_para_exp,SS_exp,alpha)
# FIT MIXED MODEL
# Define the model string entailing that we are concentrating on the exponential
# model
model_str = "mixed"
# Fit the model to the data
R_hat_mixed, opt_para_mixed, SS_mixed = fit_to_data.PE_risk_profiles(t,R,model_str,para_mixed)
# Set the parameters for the mixed model instead
A_opt_mixed_Sam = np.exp(4.877637595004758)
#A_opt_mixed_Sam = np.exp(4.878)
tau_opt_mixed_Sam = 58.40023718059162
#tau_opt_mixed_Sam = 58.400
R_hat_mixed_Sam = np.array([fit_to_data.objective_mixed(t_val,A_opt_mixed_Sam,tau_opt_mixed_Sam,alpha) for t_val in list(t)])
R_adj_new = fit_to_data.calculate_R_adj(R,R_hat_mixed_Sam,"mixed")
# FIT POWER LAW
# Define the model string entailing that we are concentrating on the exponential
# model
model_str = "power_law"
# Fit the model to the data
#print("Model 2")
R_hat_pow, opt_para_PLM, SS_PLM = fit_to_data.PE_risk_profiles(t,R,model_str,para_pow)

#---------------------------------------------------------------------------------
# PLOT MODELS
# Define the string were the plot will be stored
file_str = "../Figures/Fig1/Input/exponential.tex"
# Define the string defining the settings for the plot
plot_str = "color=exp_1,line width=2pt,"
# Define the string with the legend
legend_str = "Exponential model, $R(t)$"
# Take the exponential of the fitted models
R_hat_exp_plot = list(R_hat_exp)
R_hat_exp_plot = [np.exp(R_val) for R_val in R_hat_exp_plot]
# Take the exponential of the fitted models
R_hat_pow_plot = list(R_hat_pow)
R_hat_pow_plot = [np.exp(R_val) for R_val in R_hat_pow_plot]
# Take the exponential of the fitted models
R_hat_mixed_plot = list(R_hat_mixed)
R_hat_mixed_plot = [np.exp(R_val) for R_val in R_hat_mixed_plot]
#=====================================================================================================
# Pyplot
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=20)    # fontsize of the x and y label
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
# Subplot 1
#axes[0].plot(t, R_plot, '*', color='black', label='Data colon cancer')
#axes[0].plot(t, R_hat_exp, '-', color='red',label='Exponential model')
#axes[0].legend()
# Subplot 2
axes[0].plot(t, R_plot, '*', color='black', label='Data colon cancer')
axes[0].plot(t, R_hat_pow_plot, '-', color='blue',label='PLM')
axes[0].legend()
# Subplot 3
axes[1].plot(t, R_plot, '*', color='black', label='Data colon cancer')
axes[1].plot(t, R_hat_mixed_plot, '-', color='green',label='IM-II')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
plt.show()

#=====================================================================================================
