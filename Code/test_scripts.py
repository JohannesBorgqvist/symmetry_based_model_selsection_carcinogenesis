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
#file_name = "Colon_cancer"
# Read the data
#xlabel_str, ylabel_str, t, R = read_data.time_series_from_csv(file_name)
# Remove first 12 data points
#t = t[12:len(t)-1]
#R = R[12:len(R)-1]
# Logarithm of the output
#for i in range(len(R)):
    #R[i] = np.log(R[i])
#---------------------------------------------------------------------------------
# PLOT DATA
# Define the string were the plot will be stored
#file_str = "../Figures/Fig1/Input/colon_cancer.tex"
#file_str = "../Figures/FigS4/Input/colon_cancer.tex"
# Define the string defining the settings for the plot
#plot_str = "only marks,scatter,mark=halfcircle*,mark size=2.9pt,color=black"
#plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Define the string with the legend
#legend_str = "Data colon cancer"
# Define the plot
#R_plot = [np.exp(R_val) for R_val in R]
#write_output.plot_LaTeX_2D(t,R_plot,file_str,plot_str,legend_str)
#write_output.plot_LaTeX_2D(t,R,file_str,plot_str,legend_str)
#---------------------------------------------------------------------------------
# DEFINE DATA STRING
#data_str = "colon"
# FIT EXPONENTIAL MODEL
# Define parameter alpha common to most models
#alpha = 0.044
# Save this parameter in a parameter vector
#para_exp = [alpha]
#para_mixed = [alpha]
#para_pow = [alpha]
# Define the model string entailing that we are concentrating on the exponential
# model
#model_str = "exponential"
# Set lower and upper bound on our parameter A
# Fit the model to the data
#print("Model 1")
#R_hat_exp, opt_para_exp, SS_exp = fit_to_data.PE_risk_profiles(t,R,model_str,para_exp)
# Save the estimated parameters
#write_output.save_data_PE(data_str,model_str,opt_para_exp,SS_exp,alpha)
# FIT MIXED MODEL
# Define the model string entailing that we are concentrating on the exponential
# model
#model_str = "mixed"
# Fit the model to the data
#R_hat_mixed, opt_para_mixed, SS_mixed = fit_to_data.PE_risk_profiles(t,R,model_str,para_mixed)
# Save the estimated parameters
#write_output.save_data_PE(data_str,model_str,opt_para_mixed,SS_mixed,alpha)
# FIT POWER LAW
# Define the model string entailing that we are concentrating on the exponential
# model
#model_str = "power_law"
# Fit the model to the data
#print("Model 2")
#R_hat_pow, opt_para_PLM, SS_PLM = fit_to_data.PE_risk_profiles(t,R,model_str,para_pow)
#print("Optimal parameters PLM")
#print(opt_para_PLM)
#print("PLM fit\t=\t%s"%(str(SS_PLM)))
#print("IM-II fit\t=\t%s"%(str(SS_mixed)))

# Save the estimated parameters
#write_output.save_data_PE(data_str,model_str,opt_para_PLM,SS_PLM,alpha)

#---------------------------------------------------------------------------------
# PLOT MODELS
# Define the string were the plot will be stored
#file_str = "../Figures/Fig1/Input/exponential.tex"
# Define the string defining the settings for the plot
#plot_str = "color=exp_1,line width=2pt,"
# Define the string with the legend
#legend_str = "Exponential model, $R(t)$"
# Take the exponential of the fitted models
#R_hat_exp_plot = list(R_hat_exp)
#R_hat_exp_plot = [np.exp(R_val) for R_val in R_hat_exp_plot]
# Define the plot
#write_output.plot_LaTeX_2D(t,R_hat_exp_plot,file_str,plot_str,legend_str)
#write_output.plot_LaTeX_2D(t,R_hat_exp,file_str,plot_str,legend_str)
# Define the string were the plot will be stored
#file_str = "../Figures/Fig1/Input/power_law.tex"
# Define the string defining the settings for the plot
#plot_str = "color=pow_1,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM"
# Take the exponential of the fitted models
#R_hat_pow_plot = list(R_hat_pow)
#R_hat_pow_plot = [np.exp(R_val) for R_val in R_hat_pow_plot]
# Define the plot
#write_output.plot_LaTeX_2D(t,R_hat_pow_plot,file_str,plot_str,legend_str)
#write_output.plot_LaTeX_2D(t,R_hat_pow,file_str,plot_str,legend_str)
# Define the string were the plot will be stored
#file_str = "../Figures/Fig1/Input/mixed_model.tex"
# Define the string defining the settings for the plot
#plot_str = "color=mixed_1,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-II"
# Take the exponential of the fitted models
#R_hat_mixed_plot = list(R_hat_mixed)
#R_hat_mixed_plot = [np.exp(R_val) for R_val in R_hat_mixed_plot]
# Define the plot
#write_output.plot_LaTeX_2D(t,R_hat_mixed_plot,file_str,plot_str,legend_str)
#write_output.plot_LaTeX_2D(t,R_hat_mixed,file_str,plot_str,legend_str)
# Pyplot
# Overall properties
#fig, axes = plt.subplots(1,3,figsize=(15,5))
#plt.rc('axes', labelsize=20)    # fontsize of the x and y label
#plt.rc('legend', fontsize=15)    # legend fontsize
#plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
# Subplot 1
#axes[0].plot(t, R, '*', color='black', label='Data colon cancer')
#axes[0].plot(t, R_hat_exp, '-', color='red',label='Exponential model')
#axes[0].legend()
# Subplot 2
#axes[1].plot(t, R, '*', color='black', label='Data colon cancer')
#axes[1].plot(t, R_hat_pow, '-', color='blue',label='Power law')
#axes[1].legend()
# Subplot 3
#axes[2].plot(t, R, '*', color='black', label='Data colon cancer')
#axes[2].plot(t, R_hat_mixed, '-', color='green',label='Mixed model')
#axes[2].legend()
# add a big axis, hide frame
#fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
#plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
#plt.xlabel("Age, $t$")
#plt.ylabel("Incidence, $R(t)$")
#plt.show()


#---------------------------------------------------------------------------------
# PLOT TRANSFORMATIONS
#-----------------------------------------------------------------
# We will not have a legend for most transformations
#legend_str_empty = []
# Allocate a new more dense time vector
#t_new = np.linspace(12,t[len(t)-1],100)
#-----------------------------------------------------------------
# IM-I SYMMETRY 1
#-----------------------------------------------------------------
# Decide the indicator, meaning we look at the second symmetry
#indicator = 1
# Allocate memory for the output vector
#R_hat_exp = np.zeros(t_new.shape)
# Extract the optimal parameters
#A_opt_exp = opt_para_exp[0]
# Now we calculate the original time series
#for index in range(len(R_hat_exp)):
    #R_hat_exp[index] = np.exp(fit_to_data.objective_exponential(t_new[index],A_opt_exp,alpha))
# Plot the original time series
# Define the string were the plot will be stored
#file_str = "../Figures/FigS1/Input/IM-I_R.tex"
# Define the string defining the settings for the plot
#plot_str = "color=exp_1,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-I original solution, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_new,R_hat_exp,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 0.5
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{3,1}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "}\\right)$"
#legend_str = "Action of symmetry $\\Gamma_{3,2}\\left(\\epsilon=" + str(epsilon) + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.IM_I_transformation(t_new[50],R_hat_exp[50],epsilon,alpha,indicator)
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(52,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_I_transformation(t_new[index],R_hat_exp[index],epsilon,alpha,indicator)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the transformed solution
#t_hat_1 = []
#R_exp_hat_1 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.IM_I_1_symmetry(t_new[index],R_hat_exp[index],epsilon,alpha)
    # Append the transformed values
    #t_hat_1.append(t_hat)
    #R_exp_hat_1.append(R_hat)
# Define the string defining the settings for the plot
#plot_str = "color=exp_2,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-I transformed solution, $\hat{R}_1(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_1,R_exp_hat_1,file_str,plot_str,legend_str)
# Re-define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Loop over and add all other transformations
#for index in range(50,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_I_transformation(t_hat_1[index],R_exp_hat_1[index],epsilon,alpha,indicator)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the second transformed solution
#t_hat_2 = []
#R_exp_hat_2 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.IM_I_1_symmetry(t_hat_1[index],R_exp_hat_1[index],epsilon,alpha)
    # Append the transformed values
    #t_hat_2.append(t_hat)
    #R_exp_hat_2.append(R_hat)
# Define the string defining the settings for the plot
#plot_str = "color=exp_3,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-I transformed solution, $\hat{R}_2(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_2,R_exp_hat_2,file_str,plot_str,legend_str)
#-----------------------------------------------------------------
# IM-I SYMMETRY 2
#-----------------------------------------------------------------
# Decide the indicator meaning we look at the second symmetry
#indicator = 2
# Plot the original time series
# Define the string were the plot will be stored
#file_str = "../Figures/FigS1/Input/IM-I_t.tex"
# Define the string defining the settings for the plot
#plot_str = "color=exp_1,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-I original solution, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_new,R_hat_exp,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 15
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{3,1}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "}\\right)$"
#legend_str = "Action of symmetry $\\Gamma_{3,1}\\left(\\epsilon=" + str(epsilon) + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.IM_I_transformation(t_new[50],R_hat_exp[50],epsilon,alpha,indicator)
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(52,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_I_transformation(t_new[index],R_hat_exp[index],epsilon,alpha,indicator)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the transformed solution
#t_hat_1 = []
#R_exp_hat_1 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.IM_I_2_symmetry(t_new[index],R_hat_exp[index],epsilon)
    # Append the transformed values
    #t_hat_1.append(t_hat)
    #R_exp_hat_1.append(R_hat)
# Define the string defining the settings for the plot
#plot_str = "color=exp_2,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-I transformed solution, $\hat{R}_1(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_1,R_exp_hat_1,file_str,plot_str,legend_str)
# Re-define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Loop over and add all other transformations
#for index in range(50,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_I_transformation(t_hat_1[index],R_exp_hat_1[index],epsilon,alpha,indicator)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the second transformed solution
#t_hat_2 = []
#R_exp_hat_2 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.IM_I_2_symmetry(t_hat_1[index],R_exp_hat_1[index],epsilon)
    # Append the transformed values
    #t_hat_2.append(t_hat)
    #R_exp_hat_2.append(R_hat)
# Define the string defining the settings for the plot
#plot_str = "color=exp_3,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-I transformed solution, $\hat{R}_2(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_2,R_exp_hat_2,file_str,plot_str,legend_str)    
#-----------------------------------------------------------------
# PLM SYMMETRY 1
#-----------------------------------------------------------------
# Decide the indicator meaning we look at the second symmetry
#indicator = 1
# Allocate memory for the output vector
#R_hat_pow = np.zeros(t_new.shape)
# Extract the optimal parameters
#A_opt_PLM = opt_para_PLM[0]
#gamma_opt_PLM = opt_para_PLM[1]
# Now we calculate the original time series
#for index in range(len(R_hat_pow)):
    #R_hat_pow[index] = np.exp(fit_to_data.objective_power_law(t_new[index],A_opt_PLM,gamma_opt_PLM))
    #R_hat_pow[index] = fit_to_data.objective_power_law(t_new[index],A_opt_PLM,gamma_opt_PLM)
# Plot the original time series
# Define the string were the plot will be stored
#file_str = "../Figures/FigS2/Input/PLM_R.tex"
# Define the string defining the settings for the plot
#plot_str = "color=pow_1,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM original solution, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_new,R_hat_pow,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 0.00000025
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{1,2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "}\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[50],R_hat_pow[50],epsilon,gamma_opt_PLM,indicator)
#t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[50],np.exp(R_hat_pow[50]),epsilon,gamma_opt_PLM,indicator)
# Logarithm of transformation
#for index in range(len(R_trans)):
    #R_trans[index] = np.log(R_trans[index])
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(52,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[index],R_hat_pow[index],epsilon,gamma_opt_PLM,indicator)
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[index],np.exp(R_hat_pow[index]),epsilon,gamma_opt_PLM,indicator)
    # Logarithm of transformation
    #for index in range(len(R_trans)):
    #    R_trans[index] = np.log(R_trans[index])
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the transformed solution
#t_hat_1 = []
#R_pow_hat_1 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.PLM_1_symmetry(t_new[index],R_hat_pow[index],epsilon,gamma_opt_PLM)
    #t_hat, R_hat = symmetry_based_model_selection.PLM_1_symmetry(t_new[index],np.exp(R_hat_pow[index]),epsilon,gamma_opt_PLM)
    # Append the transformed values
    #t_hat_1.append(t_hat)
    #R_pow_hat_1.append(R_hat)
    #R_pow_hat_1.append(np.log(R_hat))
# Define the string defining the settings for the plot
#plot_str = "color=pow_2,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM transformed solution, $\hat{R}_1(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_1,R_pow_hat_1,file_str,plot_str,legend_str)
# Re-define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Loop over and add all other transformations
#for index in range(50,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_hat_1[index],R_pow_hat_1[index],epsilon,gamma_opt_PLM,indicator)
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_hat_1[index],np.exp(R_pow_hat_1[index]),epsilon,gamma_opt_PLM,indicator)
    # Logarithm of transformation
    #for index in range(len(R_trans)):
    #    R_trans[index] = np.log(R_trans[index])
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the second transformed solution
#t_hat_2 = []
#R_pow_hat_2 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.PLM_1_symmetry(t_hat_1[index],R_pow_hat_1[index],epsilon,gamma_opt_PLM)
    #t_hat, R_hat = symmetry_based_model_selection.PLM_1_symmetry(t_hat_1[index],np.exp(R_pow_hat_1[index]),epsilon,gamma_opt_PLM)
    # Append the transformed values
    #t_hat_2.append(t_hat)
    #R_pow_hat_2.append(R_hat)
    #R_pow_hat_2.append(np.log(R_hat))
# Define the string defining the settings for the plot
#plot_str = "color=pow_3,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM transformed solution, $\hat{R}_2(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_2,R_pow_hat_2,file_str,plot_str,legend_str)    
#-----------------------------------------------------------------
# PLM SYMMETRY 2
#-----------------------------------------------------------------
# Decide the indicator meaning we look at the second symmetry
#indicator = 2
# Plot the original time series
# Define the string were the plot will be stored
#file_str = "../Figures/Fig2/Input/PLM_2.tex"
#file_str = "../Figures/FigS2/Input/PLM_t.tex"
# Define the string defining the settings for the plot
#plot_str = "color=pow_1,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM original solution, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_new,R_hat_pow,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 0.2
# Create an epsilon delta so that it the arrows end before the new curve.
# This parameter is there purely for esthetical reasons...
#epsilon_delta = 0.001
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{1,2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "}\\right)$"
#legend_str = "Action of symmetry $\\Gamma_{1,1}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[50],R_hat_pow[50],epsilon-epsilon_delta,gamma_opt_PLM,indicator)
#t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[50],np.exp(R_hat_pow[50]),epsilon-epsilon_delta,gamma_opt_PLM,indicator)
# Logarithm of transformation
#for index in range(len(R_trans)):
    #R_trans[index] = np.log(R_trans[index])
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(52,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[index],R_hat_pow[index],epsilon-epsilon_delta,gamma_opt_PLM,indicator)
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_new[index],np.exp(R_hat_pow[index]),epsilon-epsilon_delta,gamma_opt_PLM,indicator)
    # Logarithm of transformation
    #for index in range(len(R_trans)):
    #    R_trans[index] = np.log(R_trans[index])    
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the transformed solution
#t_hat_1 = []
#R_pow_hat_1 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.PLM_2_symmetry(t_new[index],R_hat_pow[index],epsilon)
    # Append the transformed values
    #t_hat_1.append(t_hat)
    #R_pow_hat_1.append(R_hat)
    #R_pow_hat_1.append(np.log(R_hat))
# Define the string defining the settings for the plot
#plot_str = "color=pow_2,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM transformed solution, $\hat{R}_1(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_1,R_pow_hat_1,file_str,plot_str,legend_str)
# Re-define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Loop over and add all other transformations
#for index in range(50,len(t_new),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_hat_1[index],R_pow_hat_1[index],epsilon-epsilon_delta,gamma_opt_PLM,indicator)
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_hat_1[index],np.exp(R_pow_hat_1[index]),epsilon-epsilon_delta,gamma_opt_PLM,indicator)
    # Logarithm of transformation
    #for index in range(len(R_trans)):
        #R_trans[index] = np.log(R_trans[index])    
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the second transformed solution
#t_hat_2 = []
#R_pow_hat_2 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.PLM_2_symmetry(t_hat_1[index],R_pow_hat_1[index],epsilon)
    # Append the transformed values
    #t_hat_2.append(t_hat)
    #R_pow_hat_2.append(R_hat)
    #R_pow_hat_2.append(np.log(R_hat))
# Define the string defining the settings for the plot
#plot_str = "color=pow_3,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM transformed solution, $\hat{R}_2(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_2,R_pow_hat_2,file_str,plot_str,legend_str)    
#-----------------------------------------------------------------
# IM-II
#-----------------------------------------------------------------
# Allocate memory for the output vector
#R_hat_mixed = np.zeros(t_new.shape)
# Extract the optimal parameters
#A_opt_mixed = opt_para_mixed[0]
#tau_opt_mixed = opt_para_mixed[1]
# Now we calculate the original time series
#for index in range(len(R_hat_mixed)):
#    R_hat_mixed[index] = np.exp(fit_to_data.objective_mixed(t_new[index],A_opt_mixed,tau_opt_mixed,alpha))
# This time we take the exponentiation of the time series
#for index in range(len(R_hat_mixed)):
    #R_hat_mixed[index] = np.exp(R_hat_mixed[index])
# Define the string were the plot will be stored
#file_str = "../Figures/FigS3/Input/IM-II_t.tex"
# Define the string defining the settings for the plot
#plot_str = "color=mixed_1,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-II original solution, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_new,R_hat_mixed,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 5
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{1,2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "}\\right)$"
#legend_str = "Action of symmetry $\\Gamma_{2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t_new[80],R_hat_mixed[80],epsilon,tau_opt_mixed,alpha)

# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(81,len(t_new),1):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t_new[index],R_hat_mixed[index],epsilon,tau_opt_mixed,alpha)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the transformed solution
#t_hat_1 = []
#R_mixed_hat_1 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.IM_II_symmetry(t_new[index],R_hat_mixed[index],epsilon,tau_opt_mixed,alpha)
    # Append the transformed values
    #t_hat_1.append(t_hat)
    #R_mixed_hat_1.append(R_hat)
# Define the string defining the settings for the plot
#plot_str = "color=mixed_2,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-II transformed solution, $\hat{R}_1(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_1,R_mixed_hat_1,file_str,plot_str,legend_str)
# Re-define plot string
#plot_str = "color=black,->,>=latex,densely dashed"
# Loop over and add all other transformations
#for index in range(80,len(t_new),1):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t_hat_1[index],R_mixed_hat_1[index],epsilon,tau_opt_mixed,alpha)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Calculate the transformed solution
#t_hat_2 = []
#R_mixed_hat_2 = []
# Loop over values for the original function
#for index in range(len(t_new)):
    # Calculate the transformed values
    #t_hat, R_hat = symmetry_based_model_selection.IM_II_symmetry(t_hat_1[index],R_mixed_hat_1[index],epsilon,tau_opt_mixed,alpha)
    # Append the transformed values
    #t_hat_2.append(t_hat)
    #R_mixed_hat_2.append(R_hat)
# Define the string defining the settings for the plot
#plot_str = "color=mixed_3,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-II transformed solution, $\hat{R}_2(\hat{t})$"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_hat_2,R_mixed_hat_2,file_str,plot_str,legend_str)
#-----------------------------------------------------------------
# NUMERICAL CHECK OF SYMMETRIES
#-----------------------------------------------------------------
#-----------------------------------------------------------------
# PLM 
#-----------------------------------------------------------------
# Add the parameter gamma to our parameters
#para_pow.append(gamma_opt_PLM)
# Allocate a vector of epsilon vectors
#epsilon_vec = np.linspace(0.0, 10.0, num=20)
# Allocate a list of SS
#SS_pow = []
# Define the model string
#model_str = "power_law"
# Loop over epsilon values and calculate the sum of squares
#for eps_ind in range(len(epsilon_vec)):
    # Extract an individual epsilon value
    #epsilon = epsilon_vec[eps_ind]
    # Create a new time series
    #t_pow_dat = []
    #R_pow_dat = []
    # Loop over values for the original function
    #for index in range(len(t_new)):
        # Calculate the transformed values
        #t_hat, R_hat = symmetry_based_model_selection.PLM_2_symmetry(t_new[index],R_hat_pow[index],epsilon)
        # Append the transformed values
        #t_pow_dat.append(t_hat)
        #R_pow_dat.append(R_hat)
    # Convert the data to array
    #t_pow_dat = np.array(t_pow_dat)
    #R_pow_dat = np.array(R_pow_dat)
    # Calculate the fit
    #R_hat_pow_new, opt_para_PLM, SS_PLM = fit_to_data.PE_risk_profiles(t_pow_dat,np.log(R_pow_dat),model_str,para_pow)
    # Append the fit
    #SS_pow.append(SS_PLM)
# Afterwards we convert the list to an array
#SS_pow = np.array(SS_pow)
# Define the string were the plot will be stored
#file_str = "../Figures/FigS5/Input/PLM_2.tex"
# Define the string defining the settings for the plot
#plot_str = "color=pow_1,line width=2pt,"
# Define the string with the legend
#legend_str = "PLM"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(epsilon_vec,SS_pow,file_str,plot_str,legend_str)
#-----------------------------------------------------------------
# IM-II
#-----------------------------------------------------------------
# Add the parameter gamma to our parameters
#para_mixed.append(tau_opt_mixed)
# Allocate a list of SS
#SS_mixed = []
# Define the model string
#model_str = "mixed"
# Extract the value of tau
#tau = opt_para_mixed[1]
# Loop over epsilon values and calculate the sum of squares
#for eps_ind in range(len(epsilon_vec)):
    # Extract an individual epsilon value
    #epsilon = epsilon_vec[eps_ind]
    # Create a new time series
    #t_mixed_dat = []
    #R_mixed_dat = []
    # Loop over values for the original function
    #for index in range(len(t_new)):
        # Calculate the transformed values
        #t_hat, R_hat = symmetry_based_model_selection.IM_II_symmetry(t_new[index],R_hat_mixed[index],epsilon,tau_opt_mixed,alpha)
        # Append the transformed values
        #t_mixed_dat.append(t_hat)
        #R_mixed_dat.append(R_hat)
    # Convert the data to array
    #t_mixed_dat = np.array(t_pow_dat)
    #R_mixed_dat = np.array(R_pow_dat)
    # Calculate the fit
    #R_hat_mixed_new, opt_para_mixed, SS_IM = fit_to_data.PE_risk_profiles(t_mixed_dat,np.log(R_mixed_dat),model_str,para_mixed)
    # Append the fit
    #SS_mixed.append(SS_IM)
# Afterwards we convert the list to an array
#SS_mixed = np.array(SS_mixed)
# Define the string were the plot will be stored
#file_str = "../Figures/FigS5/Input/IM-II.tex"
# Define the string defining the settings for the plot
#plot_str = "color=mixed_1,line width=2pt,"
# Define the string with the legend
#legend_str = "IM-II"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(epsilon_vec,SS_pow,file_str,plot_str,legend_str)
#---------------------------------------------------------------------------------
# ILLUSTRATE MODEL SELECTION POWER LAW
# STEP 1: TRANSFORMED DATA
# Define the string were the plot will be stored
#file_str = "../Figures/FigS6/Input/step_1.tex"
# Define the string defining the settings for the plot
#plot_str = "only marks,scatter,mark=halfcircle*,mark size=2.9pt,color=black"
#plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Define the string with the legend
#legend_str = "Original data"
# Re-define data
#R_ori_data = [np.exp(R[i]) for i in range(len(R))]
# Define the plot
#write_output.plot_LaTeX_2D(t,R_ori_data,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 0.2
# Define a small epsilon_delta for visual purposes
#epsilon_delta = 0.001
# Make sure that we use the t-directional symmetry
#indicator = 2
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{1,2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "}\\right)$"
#legend_str = "Action of symmetry $\\Gamma_{1,1}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t[10],R_ori_data[10],epsilon-epsilon_delta,gamma_opt_PLM,indicator)
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(12,len(t),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t[index],R_ori_data[index],epsilon-epsilon_delta,gamma_opt_PLM,indicator)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Define the transformed data
#R_trans_data = []
#t_trans_data = []
# Loop over the original time series and calculate the transformed time series
#for index in range(len(R_ori_data)):
    # Calculate each transformed coordinate
    #t_hat_temp, R_hat_temp = symmetry_based_model_selection.PLM_2_symmetry(t[index],R_ori_data[index],epsilon)
    # Save the transformed coordinates
    #t_trans_data.append(t_hat_temp)
    #R_trans_data.append(R_hat_temp)
#plot_str = "only marks, mark=diamond*,mark size=1.5pt,color=gray,every mark/.append style={solid, fill=gray}"
# Define the string with the legend
#legend_str = "Transformed data"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_trans_data,file_str,plot_str,legend_str)
# STEP 2: FIT TO TRANSFORMED DATA
# Define the string were the plot will be stored
#file_str = "../Figures/FigS6/Input/step_2.tex"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_trans_data,file_str,plot_str,legend_str)
# Define the logarithm of the data
#R_trans_data_log = [np.log(R_trans) for R_trans in R_trans_data]
# Convert the transformed time series to arrays
#t_trans_data = np.array(t_trans_data)
#R_trans_data_log = np.array(R_trans_data_log)
# FIT POWER LAW
# Define the model string entailing that we are concentrating on the exponential
# model
#model_str = "power_law"
# Fit the model to the data
#print("Model 2")
#R_hat_pow, opt_para_PLM, SS_PLM = fit_to_data.PE_risk_profiles(t_trans_data,R_trans_data_log,model_str,para_pow)
# Convert these to lists
#R_hat_pow = R_hat_pow.tolist()
#t_trans_data = t_trans_data.tolist()
# Take the exponential of R
#R_hat_pow = [np.exp(R_hat) for R_hat in R_hat_pow]
# Define the string defining the settings for the plot
#plot_str = "color=pow_2,line width=2pt,"
# Define the string with the legend
#legend_str = "Fitted PLM, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_hat_pow,file_str,plot_str,legend_str)
# STEP 3: Inversely transform the model back
# Define the string were the plot will be stored
#file_str = "../Figures/FigS6/Input/step_3.tex"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_hat_pow,file_str,plot_str,legend_str)
#legend_str = "Action of inverse symmetry $\\Gamma^{-1}_{1,1}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_trans_data[10],R_hat_pow[10],-(epsilon-epsilon_delta),gamma_opt_PLM,indicator)
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(12,len(t_trans_data),2):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.PLM_transformation(t_trans_data[index],R_hat_pow[index],-(epsilon-epsilon_delta),gamma_opt_PLM,indicator)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Define the inversely transformed model
#R_inv_trans_model = []
#t_inv_trans_model = []
# Loop over the original time series and calculate the inversely transformed model
#for index in range(len(R_hat_pow)):
    # Calculate each transformed coordinate
    #t_hat_temp, R_hat_temp = symmetry_based_model_selection.PLM_2_symmetry(t_trans_data[index],R_hat_pow[index],-epsilon)
    # Save the transformed coordinates
    #t_inv_trans_model.append(t_hat_temp)
    #R_inv_trans_model.append(R_hat_temp)
# Define the string defining the settings for the plot
#plot_str = "color=pow_1,line width=2pt,"
# Define the string with the legend
#legend_str = "Inversely transformed PLM, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_inv_trans_model,R_inv_trans_model,file_str,plot_str,legend_str)    
# STEP 4: Compare inversely fitted model to original data
# Define the string were the plot will be stored
#file_str = "../Figures/FigS6/Input/step_4.tex"
# Define the plot
#write_output.plot_LaTeX_2D(t_inv_trans_model,R_inv_trans_model,file_str,plot_str,legend_str)    
# Define the string defining the settings for the plot
#plot_str = "only marks,scatter,mark=halfcircle*,mark size=2.9pt,color=black"
#plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Define the string with the legend
#legend_str = "Original data"
# Define the plot
#write_output.plot_LaTeX_2D(t,R_ori_data,file_str,plot_str,legend_str)
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# ILLUSTRATE MODEL SELECTION IM-II
# STEP 1: TRANSFORMED DATA
# Define the string were the plot will be stored
#file_str = "../Figures/FigS7/Input/step_1.tex"
# Define the string defining the settings for the plot
#plot_str = "only marks,scatter,mark=halfcircle*,mark size=2.9pt,color=black"
#plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Define the string with the legend
#legend_str = "Original data"
# Re-define data
#R_ori_data = [np.exp(R[i]) for i in range(len(R))]
# Define the plot
#write_output.plot_LaTeX_2D(t,R_ori_data,file_str,plot_str,legend_str)
# Define a transformation parameter
#epsilon = 7
# Define a small epsilon_delta for visual purposes
#epsilon_delta = 0.001
# Define the legend as well
#legend_str = "Action of symmetry $\\Gamma_{2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t[48],R_ori_data[48],epsilon,tau_opt_mixed,alpha)
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(49,len(t),1):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t[index],R_ori_data[index],epsilon,tau_opt_mixed,alpha)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Define the transformed data
#R_trans_data = []
#t_trans_data = []
# Loop over the original time series and calculate the transformed time series
#for index in range(len(R_ori_data)):
    # Calculate each transformed coordinate
    #t_hat_temp, R_hat_temp = symmetry_based_model_selection.IM_II_symmetry(t[index],R_ori_data[index],epsilon,tau_opt_mixed,alpha)
    # Save the transformed coordinates
    #t_trans_data.append(t_hat_temp)
    #R_trans_data.append(R_hat_temp)
#plot_str = "only marks, mark=diamond*,mark size=1.5pt,color=gray,every mark/.append style={solid, fill=gray}"
# Define the string with the legend
#legend_str = "Transformed data"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_trans_data,file_str,plot_str,legend_str)
# STEP 2: FIT TO TRANSFORMED DATA
# Define the string were the plot will be stored
#file_str = "../Figures/FigS7/Input/step_2.tex"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_trans_data,file_str,plot_str,legend_str)
# Define the logarithm of the data
#R_trans_data_log = [np.log(R_trans) for R_trans in R_trans_data]
# Convert the transformed time series to arrays
#t_trans_data = np.array(t_trans_data)
#R_trans_data_log = np.array(R_trans_data_log)
# FIT IM-II
# Define the model string entailing that we are concentrating on the exponential
# model
#model_str = "mixed"
# Fit the model to the data
#print("Model 2")
#R_hat_mixed, opt_para_IM_II, SS_IM_II = fit_to_data.PE_risk_profiles(t_trans_data,R_trans_data_log,model_str,para_mixed)
# Convert these to lists
#R_hat_mixed = R_hat_mixed.tolist()
#t_trans_data = t_trans_data.tolist()
# Take the exponential of R
#R_hat_mixed = [np.exp(R_hat) for R_hat in R_hat_mixed]
# Define the string defining the settings for the plot
#plot_str = "color=mixed_2,line width=2pt,"
# Define the string with the legend
#legend_str = "Fitted IM-II, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_hat_mixed,file_str,plot_str,legend_str)
# STEP 3: Inversely transform the model back
# Define the string were the plot will be stored
#file_str = "../Figures/FigS7/Input/step_3.tex"
# Define the plot
#write_output.plot_LaTeX_2D(t_trans_data,R_hat_mixed,file_str,plot_str,legend_str)
#legend_str = "Action of inverse symmetry $\\Gamma^{-1}_{2}\\left(\\epsilon=" + str(epsilon).replace("e-0","\\cdot 10^{-") + "\\right)$"
# Transform the time series 
#t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t_trans_data[48],R_hat_mixed[48],-epsilon,tau_opt_mixed,alpha)
# Define the plot string, defining the colour
#plot_str = "color=black,->,>=latex,densely dashed"
# Save the transformed stuff
#write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str)
# Loop over and add all other transformations
#for index in range(49,len(t_trans_data),1):
    # Transform the time series 
    #t_trans,R_trans = symmetry_based_model_selection.IM_II_transformation(t_trans_data[index],R_hat_mixed[index],-epsilon,tau_opt_mixed,alpha)
    # Save the transformed stuff
    #write_output.plot_LaTeX_2D(t_trans,R_trans,file_str,plot_str,legend_str_empty)
# Define the inversely transformed model
#R_inv_trans_model = []
#t_inv_trans_model = []
# Loop over the original time series and calculate the inversely transformed model
#for index in range(len(R_hat_mixed)):
    # Calculate each transformed coordinate
    #t_hat_temp, R_hat_temp = symmetry_based_model_selection.IM_II_symmetry(t_trans_data[index],R_hat_mixed[index],-epsilon,tau_opt_mixed,alpha)    
    # Save the transformed coordinates
    #t_inv_trans_model.append(t_hat_temp)
    #R_inv_trans_model.append(R_hat_temp)
# Define the string defining the settings for the plot
#plot_str = "color=mixed_1,line width=2pt,"
# Define the string with the legend
#legend_str = "Inversely transformed IM-II, $R(t)$"
# Define the plot
#write_output.plot_LaTeX_2D(t_inv_trans_model,R_inv_trans_model,file_str,plot_str,legend_str)    
# STEP 4: Compare inversely fitted model to original data
# Define the string were the plot will be stored
#file_str = "../Figures/FigS7/Input/step_4.tex"
# Define the plot
#write_output.plot_LaTeX_2D(t_inv_trans_model,R_inv_trans_model,file_str,plot_str,legend_str)    
# Define the string defining the settings for the plot
#plot_str = "only marks,scatter,mark=halfcircle*,mark size=2.9pt,color=black"
#plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Define the string with the legend
#legend_str = "Original data"
# Define the plot
#write_output.plot_LaTeX_2D(t,R_ori_data,file_str,plot_str,legend_str)
#---------------------------------------------------------------------------------
# Define an epsilon vector
#epsilon_vector = list(np.linspace(0.0,10.0,num=100))
# PLM: Calculate the relative fit
#Delta_fit_pow = symmetry_based_model_selection.sym_model_sel(t,R,epsilon_vector,para_pow,"power_law")
# IM_II: Calculate the relative fit
#Delta_fit_mixed = symmetry_based_model_selection.sym_model_sel(t,R,epsilon_vector,para_mixed,"mixed")
# Save the power law
#write_output.plot_LaTeX_2D(epsilon_vector,Delta_fit_pow,"../Figures/Fig3/Input/PLM.tex","color=pow_2,line width=2pt,","PLM")
# Save the mixed
#write_output.plot_LaTeX_2D(epsilon_vector,Delta_fit_mixed,"../Figures/Fig3/Input/IM_II.tex","color=mixed_2,line width=2pt,","IM-II")



#epsilon_test = [0]
#Delta_fit_test = symmetry_based_model_selection.sym_model_sel(t,R,epsilon_test,para_pow,"power_law")
#Delta_fit_test_2 = symmetry_based_model_selection.sym_model_sel(t,R,epsilon_test,para_mixed,"mixed")
#print("Pow:\t(%0.3f,%0.3f)"%(epsilon_test[0],Delta_fit_test[0]))
#print("Pow:\t(%0.3f,%0.3f)"%(epsilon_test[0],Delta_fit_test_2[0]))
