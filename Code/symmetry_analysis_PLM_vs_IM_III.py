# =================================================================================
# =================================================================================
# Script:"symmetry_based_model_selection.py"
# Date: 2022-01-05
# Implemented by: Johannes Borgqvist
# Description:
# The script re-generates all the results presented in the article.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import read_data  # Home-made
import write_output  # Home-made
import fit_to_data # Home-made
import symmetry_toolbox # Home-made
import numpy as np
import math
from matplotlib import pyplot as plt
# =================================================================================
# =================================================================================
# Read the data
# =================================================================================
# =================================================================================
# MYELOMA
xlabel_myeloma, ylabel_myeloma, t_myeloma, R_myeloma = read_data.time_series_from_csv("Myeloma_cancer")
# Only ages above 25 as we have zero incidences below
t_myeloma = np.array(list(t_myeloma[24:len(t_myeloma)-1]))
R_myeloma = np.array(list(R_myeloma[24:len(R_myeloma)-1]))
# COLON
xlabel_colon, ylabel_colon, t_colon, R_colon = read_data.time_series_from_csv("Colon_cancer")
# Only ages above 12 as we have zero incidences below
t_colon = np.array(list(t_colon[11:len(t_colon)-1]))
R_colon = np.array(list(R_colon[11:len(R_colon)-1]))
# Chronic Myeloid Leukemia (CML)
xlabel_CML, ylabel_CML, t_CML, R_CML = read_data.time_series_from_csv("CML_cancer")
# Only ages above 10 as we have zero incidences below
t_CML = np.array(list(t_CML[9:len(t_CML)-1]))
R_CML = np.array(list(R_CML[9:len(R_CML)-1]))
# =================================================================================
# =================================================================================
# FIT THE MODELS TO THE DATA
# =================================================================================
# =================================================================================
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Orthonal Distance Regression (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_ODR, R_hat_PLM_myeloma_ODR, RMS_PLM_myeloma_ODR  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","ODR",[])
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_ODR, RMS_PLM_myeloma_ODR, "ODR")
# IM-III
IM_III_fitted_to_myeloma_ODR, R_hat_IM_III_myeloma_ODR, RMS_IM_III_myeloma_ODR = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-III","ODR",[])
write_output.save_data_PE("myeloma", "IM-III", IM_III_fitted_to_myeloma_ODR, RMS_IM_III_myeloma_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_ODR, R_hat_PLM_colon_ODR, RMS_PLM_colon_ODR  = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM","ODR",[])
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_ODR, RMS_PLM_colon_ODR, "ODR")
# IM-III
IM_III_fitted_to_colon_ODR, R_hat_IM_III_colon_ODR, RMS_IM_III_colon_ODR = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-III","ODR",[])
write_output.save_data_PE("colon", "IM-III", IM_III_fitted_to_colon_ODR, RMS_IM_III_colon_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_ODR, R_hat_PLM_CML_ODR, RMS_PLM_CML_ODR  = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM","ODR",[])
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_ODR, RMS_PLM_CML_ODR, "ODR")
# IM-III
IM_III_fitted_to_CML_ODR, R_hat_IM_III_CML_ODR, RMS_IM_III_CML_ODR = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-III","ODR",[])
write_output.save_data_PE("CML", "IM-III", IM_III_fitted_to_CML_ODR, RMS_IM_III_CML_ODR, "ODR")

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the data and the fit in Python
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT OF THE FIT OF THE MODEL TO THE DATA
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0].plot(t_myeloma, R_hat_PLM_myeloma_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[0].plot(t_myeloma, R_hat_IM_III_myeloma_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM-III')
axes[0].legend()
# Subplot 2a
axes[1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1].plot(t_colon, R_hat_PLM_colon_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[1].plot(t_colon, R_hat_IM_III_colon_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM-III')
axes[1].legend()
# Subplot 3a
axes[2].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2].plot(t_CML, R_hat_PLM_CML_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[2].plot(t_CML, R_hat_IM_III_CML_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM-III')
axes[2].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("Fit of candidate models to three data sets",fontsize=20, fontweight='bold')
#plt.show()
plt.savefig("../Figures/Fit_of_models_to_cancer_data.png")
#----------------------------------------------------------------------------------
# =================================================================================
# =================================================================================
# ILLUSTRATE THE ACTION OF THE SYMMETRIES
# =================================================================================
# =================================================================================
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE PLM
R_hat_PLM_original = np.array([fit_to_data.objective_PLM(PLM_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Choose an epsilon
epsilon = 0.40
# Allocate memory for a list
R_PLM_trans_1 = []
t_PLM_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(10,len(t_sym)-1,75,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon-0.045,PLM_fitted_to_myeloma_ODR.beta[1])
    # Save the transformed variables
    R_PLM_trans_1.append(R_trans)
    t_PLM_trans_1.append(t_trans)
# Transform the original solution
t_hat_PLM_1,R_hat_PLM_1 = symmetry_toolbox.PLM_transformed_solution(t_sym,R_hat_PLM_original,epsilon,PLM_fitted_to_myeloma_ODR.beta[0],PLM_fitted_to_myeloma_ODR.beta[1])
# Allocate memory for a list
R_PLM_trans_2 = []
t_PLM_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon-0.045,PLM_fitted_to_myeloma_ODR.beta[1])
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2,R_hat_PLM_2 = symmetry_toolbox.PLM_transformed_solution(t_hat_PLM_1,R_hat_PLM_1,epsilon,PLM_fitted_to_myeloma_ODR.beta[0],PLM_fitted_to_myeloma_ODR.beta[1])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE IM-II
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
# Original solution
R_hat_IM_III_original = np.array([fit_to_data.objective_IM_III(IM_III_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_III_trans_1 = []
t_IM_III_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(175,len(t_sym)-1,25,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_sym[index],R_hat_IM_III_original[index],epsilon-0.045,IM_III_fitted_to_myeloma_ODR.beta[1],0.044)
    # Save the transformed variables
    R_IM_III_trans_1.append(R_trans)
    t_IM_III_trans_1.append(t_trans)
# Transform the original solution
t_hat_IM_III_1,R_hat_IM_III_1 = symmetry_toolbox.IM_III_transformed_solution(t_sym,R_hat_IM_III_original,epsilon,IM_III_fitted_to_myeloma_ODR.beta[0],IM_III_fitted_to_myeloma_ODR.beta[1],IM_III_fitted_to_myeloma_ODR.beta[2])
# Allocate memory for a list
R_IM_III_trans_2 = []
t_IM_III_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_hat_IM_III_1[index],R_hat_IM_III_1[index],epsilon-0.045,IM_III_fitted_to_myeloma_ODR.beta[1],0.044)
    # Save the transformed variables
    R_IM_III_trans_2.append(R_trans)
    t_IM_III_trans_2.append(t_trans)
# Transform the second solution      
t_hat_IM_III_2,R_hat_IM_III_2 = symmetry_toolbox.IM_III_transformed_solution(t_hat_IM_III_1,R_hat_IM_III_1,epsilon,IM_III_fitted_to_myeloma_ODR.beta[0],IM_III_fitted_to_myeloma_ODR.beta[1],IM_III_fitted_to_myeloma_ODR.beta[2])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the action of the symmetries in Python
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1
axes[0].plot(t_sym, R_hat_PLM_original, '-', color = (103/256,0/256,31/256), label='Original solution PLM')
axes[0].plot(np.array(t_PLM_trans_1[0]),np.array(R_PLM_trans_1[0]),'--',color='black',label='Symmetry PLM')
for index in range(1,len(t_PLM_trans_1)):
    axes[0].plot(np.array(t_PLM_trans_1[index]),np.array(R_PLM_trans_1[index]),'--',color='black')
axes[0].plot(t_hat_PLM_1, R_hat_PLM_1, '-', color=(206/256,18/256,86/256), label='PLM transformed solution 1')
for index in range(len(t_PLM_trans_2)):
    axes[0].plot(np.array(t_PLM_trans_2[index]),np.array(R_PLM_trans_2[index]),'--',color='black')
axes[0].plot(t_hat_PLM_2, R_hat_PLM_2, '-', color=(223/256,101/256,176/256), label='PLM transformed solution 2')
axes[0].legend()
# Subplot 2
axes[1].plot(t_sym, R_hat_IM_III_original, '-', color = (2/256,56/256,88/256), label='Original solution IM-III')
axes[1].plot(np.array(t_IM_III_trans_1[0]),np.array(R_IM_III_trans_1[0]),'--',color='black',label='Symmetry IM-III')
for index in range(1,len(t_IM_III_trans_1)):
    axes[1].plot(np.array(t_IM_III_trans_1[index]),np.array(R_IM_III_trans_1[index]),'--',color='black')
axes[1].plot(t_hat_IM_III_1, R_hat_IM_III_1, '-', color=(54/256,144/256,192/256), label='IM-II transformed solution 1')
for index in range(len(t_IM_III_trans_2)):
    axes[1].plot(np.array(t_IM_III_trans_2[index]),np.array(R_IM_III_trans_2[index]),'--',color='black')
axes[1].plot(t_hat_IM_III_2, R_hat_IM_III_2, '-', color=(208/256,209/256,230/256), label='IM-III transformed solution 2')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The symmetries of the PLM and the IM-II",fontsize=20, fontweight='bold')
#plt.show()
plt.savefig("../Figures/action_of_symmetries.png")
# =================================================================================
# =================================================================================
# ILLUSTRATE THE MODEL SELECTION FRAMEWORK
# =================================================================================
# =================================================================================
# Pick a value of epsilon
epsilon = 0.60
#------------------------------------------------------------------------------
# PLM
#------------------------------------------------------------------------------
# Allocate memory for a list
R_PLM_trans = []
t_PLM_trans = []
# Allocate an index vector
index_vector = list(np.linspace(10,55,30,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_myeloma[index],R_myeloma[index],epsilon-0.015,PLM_fitted_to_myeloma_ODR.beta[1])
    # Save the transformed variables
    R_PLM_trans.append(R_trans)
    t_PLM_trans.append(t_trans)
# Transform the data
t_myeloma_trans_PLM = np.array([symmetry_toolbox.PLM_symmetry(t_myeloma[index],R_myeloma[index],epsilon)[0] for index in range(len(t_myeloma))])
R_myeloma_trans_PLM = np.array([symmetry_toolbox.PLM_symmetry(t_myeloma[index],R_myeloma[index],epsilon)[1] for index in range(len(t_myeloma))])
# Fit the model to the transformed data
#param_temp, R_hat_PLM_trans, RMS_temp  = fit_to_data.PE_risk_profiles(t_myeloma_trans_PLM,R_myeloma_trans_PLM,"PLM","ODR",[])
# Extract parameters
A_temp = PLM_fitted_to_myeloma_ODR.beta[0]
gamma_temp = PLM_fitted_to_myeloma_ODR.beta[1]
# Update A_temp
A_temp = A_temp*np.exp(-gamma_temp*epsilon)
# Calculate transformed solution
R_hat_PLM_trans = np.array([fit_to_data.objective_PLM((A_temp,gamma_temp),t_i) for t_i in list(t_myeloma_trans_PLM)])
#------------------------------------------------------------------------------
# IM-III
#------------------------------------------------------------------------------
# Allocate memory for a list
R_IM_III_trans = []
t_IM_III_trans = []
# Allocate an index vector
index_vector = list(np.linspace(30,58,30,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_myeloma[index],R_myeloma[index],epsilon-0.015,IM_III_fitted_to_myeloma_ODR.beta[1],0.044)
    # Save the transformed variables
    R_IM_III_trans.append(R_trans)
    t_IM_III_trans.append(t_trans)
# Transform the data
t_myeloma_trans_IM_III = np.array([symmetry_toolbox.IM_III_symmetry(t_myeloma[index],R_myeloma[index],epsilon,IM_III_fitted_to_myeloma_ODR.beta[1],0.044)[0] for index in range(len(t_myeloma))])
R_myeloma_trans_IM_III = np.array([symmetry_toolbox.IM_III_symmetry(t_myeloma[index],R_myeloma[index],epsilon,IM_III_fitted_to_myeloma_ODR.beta[1],0.044)[1] for index in range(len(t_myeloma))])
# IM-III CALCULATE THE TRANSFORMED CURVE
# Extract parameters
A_temp = IM_III_fitted_to_myeloma_ODR.beta[0]
tau_temp = IM_III_fitted_to_myeloma_ODR.beta[1]
C_temp = IM_III_fitted_to_myeloma_ODR.beta[2]
alpha = 0.044
# Update C_temp
C_temp += -alpha*np.exp(alpha*tau_temp)*epsilon
# Calculate transformed solution
R_hat_IM_III_trans = np.array([fit_to_data.objective_IM_III((A_temp,tau_temp,C_temp),t_i) for t_i in list(t_myeloma_trans_IM_III)])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the illustration of the model selection framework
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(2,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0][0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0][0].plot(np.array(t_PLM_trans[0]),np.array(R_PLM_trans[0]),'--',color='black',label='Symmetry PLM')
for index in range(1,len(t_PLM_trans)):
    axes[0][0].plot(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),'--',color='black')
axes[0][0].plot(t_myeloma_trans_PLM, R_myeloma_trans_PLM, '*', color='gray', label='Transformed Data')
axes[0][0].legend()
# Subplot 1b
axes[0][1].plot(t_myeloma_trans_PLM, R_myeloma_trans_PLM, '*', color='gray', label='Transformed Data')
axes[0][1].plot(t_myeloma_trans_PLM, R_hat_PLM_trans, '-', color=(206/256,18/256,86/256), label='Fitted PLM')
axes[0][1].legend()
# Subplot 2a
axes[1][0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[1][0].plot(np.array(t_IM_III_trans[0]),np.array(R_IM_III_trans[0]),'--',color='black',label='Symmetry IM-III')
for index in range(1,len(t_IM_III_trans)):
    axes[1][0].plot(np.array(t_IM_III_trans[index]),np.array(R_IM_III_trans[index]),'--',color='black')
axes[1][0].plot(t_myeloma_trans_IM_III, R_myeloma_trans_IM_III, '*', color='gray', label='Transformed Data')
axes[1][0].legend()
# Subplot 2b
axes[1][1].plot(t_myeloma_trans_IM_III, R_myeloma_trans_IM_III, '*', color='gray', label='Transformed Data')
axes[1][1].plot(t_myeloma_trans_IM_III, R_hat_IM_III_trans, '-', color=(54/256,144/256,192/256), label='Fitted IM-II')
axes[1][1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The steps of the symmetry based model selection",fontsize=20, fontweight='bold')
plt.savefig("../Figures/illustrate_framework.png")
# =================================================================================
# =================================================================================
# CONDUCT THE SYMMETRY BASED MODEL SELECTION
# =================================================================================
# =================================================================================
# Allocate an epsilon vector with transformation parameters
epsilon_vector = np.linspace(0.0,0.60,num=10,endpoint=True)
# MYELOMA CANCER
RMS_transf_PLM_myeloma = symmetry_toolbox.symmetry_based_model_selection(t_myeloma,R_myeloma,epsilon_vector,PLM_fitted_to_myeloma_ODR.beta,"PLM")
RMS_transf_IM_III_myeloma = symmetry_toolbox.symmetry_based_model_selection(t_myeloma,R_myeloma,epsilon_vector,IM_III_fitted_to_myeloma_ODR.beta,"IM-III")
# COLON CANCER
RMS_transf_PLM_colon = symmetry_toolbox.symmetry_based_model_selection(t_colon,R_colon,epsilon_vector,PLM_fitted_to_colon_ODR.beta,"PLM")
RMS_transf_IM_III_colon = symmetry_toolbox.symmetry_based_model_selection(t_colon,R_colon,epsilon_vector,IM_III_fitted_to_colon_ODR.beta,"IM-III")
# CML CANCER
RMS_transf_PLM_CML = symmetry_toolbox.symmetry_based_model_selection(t_CML,R_CML,epsilon_vector,PLM_fitted_to_CML_ODR.beta,"PLM")
RMS_transf_IM_III_CML = symmetry_toolbox.symmetry_based_model_selection(t_CML,R_CML,epsilon_vector,IM_III_fitted_to_CML_ODR.beta,"IM-III")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the symmetry based model selection
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: Myeloma
axes[0].plot(epsilon_vector,RMS_transf_PLM_myeloma,'-', color = (103/256,0/256,31/256),label='PLM Myeloma cancer')
axes[0].plot(epsilon_vector,RMS_transf_IM_III_myeloma,'-', color = (2/256,56/256,88/256),label='IM-III Myeloma cancer')
# Subplot 1: Myeloma
axes[1].plot(epsilon_vector,RMS_transf_PLM_colon,'-', color = (103/256,0/256,31/256),label='PLM Colon cancer')
axes[1].plot(epsilon_vector,RMS_transf_IM_III_colon,'-', color = (2/256,56/256,88/256),label='IM-III Colon cancer')
# Subplot 1: Myeloma
axes[2].plot(epsilon_vector,RMS_transf_PLM_CML,'-', color = (103/256,0/256,31/256),label='PLM CML')
axes[2].plot(epsilon_vector,RMS_transf_IM_III_CML,'-', color = (2/256,56/256,88/256),label='IM-III CML')
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("The symmetry based model selection",fontsize=20, fontweight='bold')
plt.savefig("../Figures/symmetry_based_model_selection.png")
plt.show()


