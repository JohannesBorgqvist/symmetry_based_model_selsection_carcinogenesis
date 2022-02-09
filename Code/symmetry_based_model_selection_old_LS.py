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
#import symmetry_based_model_selection # Home-made
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
# Fit the models to the data
# =================================================================================
# =================================================================================
# Myeloma
opt_para_PLM_myeloma, R_hat_PLM_myeloma, RMS_PLM_myeloma = fit_to_data.PE_risk_profiles(t_myeloma, R_myeloma, "PLM","all")
opt_para_IM_II_myeloma, R_hat_IM_II_myeloma, RMS_IM_II_myeloma = fit_to_data.PE_risk_profiles(t_myeloma, R_myeloma, "IM-II","not alpha")
# Colon
opt_para_PLM_colon, R_hat_PLM_colon, RMS_PLM_colon = fit_to_data.PE_risk_profiles(t_colon, R_colon, "PLM","all")
opt_para_IM_II_colon, R_hat_IM_II_colon, RMS_IM_II_colon = fit_to_data.PE_risk_profiles(t_colon, R_colon, "IM-II","not alpha")
# CML
opt_para_PLM_CML, R_hat_PLM_CML, RMS_PLM_CML = fit_to_data.PE_risk_profiles(t_CML, R_CML, "PLM","all")
opt_para_IM_II_CML, R_hat_IM_II_CML, RMS_IM_II_CML = fit_to_data.PE_risk_profiles(t_CML, R_CML, "IM-II","not alpha")
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Myeloma
axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data myeloma')
axes[0].plot(t_myeloma, R_hat_PLM_myeloma, '-',  color = (103/256,0/256,31/256), label='PML')
axes[0].plot(t_myeloma, R_hat_IM_II_myeloma, '-',  color = (2/256,56/256,88/256), label='IM_II')
axes[0].legend()
# Colon
axes[1].plot(t_colon, R_colon, '*', color='black', label='Data colon')
axes[1].plot(t_colon, R_hat_PLM_colon, '-',  color = (103/256,0/256,31/256), label='PML')
axes[1].plot(t_colon, R_hat_IM_II_colon, '-',  color = (2/256,56/256,88/256), label='IM_II')
axes[1].legend()
# CML
axes[2].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2].plot(t_CML, R_hat_PLM_CML, '-',  color = (103/256,0/256,31/256), label='PML')
axes[2].plot(t_CML, R_hat_IM_II_CML, '-',  color = (2/256,56/256,88/256), label='IM_II')
axes[2].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The steps of the symmetry based model selection",fontsize=20, fontweight='bold')
plt.savefig("../Figures/fitting_curve_fit_CML.png")
#plt.show()

# =================================================================================
# =================================================================================
# Illustrate symmetries
# =================================================================================
# =================================================================================
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE PLM
R_hat_PLM_original = np.array([fit_to_data.objective_PLM(t_sym[index],opt_para_PLM_CML[0],opt_para_PLM_CML[1]) for index in range(len(t_sym))])
# Choose an epsilon
epsilon = 0.31
# Allocate memory for a list
R_PLM_trans_1 = []
t_PLM_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(10,len(t_sym)-1,75,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon-0.01,opt_para_PLM_CML[1],2)
    # Save the transformed variables
    R_PLM_trans_1.append(R_trans)
    t_PLM_trans_1.append(t_trans)
# Transform the original solution
t_hat_PLM_1 = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon)[0] for index in range(len(t_sym))])
R_hat_PLM_1 = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon)[1] for index in range(len(t_sym))])
# Allocate memory for a list
R_PLM_trans_2 = []
t_PLM_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon-0.01,opt_para_PLM_CML[1],2)
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2 = np.array([symmetry_toolbox.PLM_2_symmetry(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon)[0] for index in range(len(t_sym))])
R_hat_PLM_2 = np.array([symmetry_toolbox.PLM_2_symmetry(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon)[1] for index in range(len(t_sym))])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE IM-II
R_hat_IM_II_original = np.array([fit_to_data.objective_IM_II(t_sym[index],opt_para_IM_II_CML[0],opt_para_IM_II_CML[1],0.044) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_II_trans_1 = []
t_IM_II_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(175,len(t_sym)-1,25,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_sym[index],R_hat_IM_II_original[index],epsilon-0.01,opt_para_IM_II_CML[1],0.044)
    # Save the transformed variables
    R_IM_II_trans_1.append(R_trans)
    t_IM_II_trans_1.append(t_trans)
# Transform the original solution
t_hat_IM_II_1 = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon,opt_para_IM_II_CML[1],0.044)[0] for index in range(len(t_sym))])
R_hat_IM_II_1 = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon,opt_para_IM_II_CML[1],0.044)[1] for index in range(len(t_sym))])    
# Allocate memory for a list
R_IM_II_trans_2 = []
t_IM_II_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon-0.01,opt_para_IM_II_CML[1],0.044)
    # Save the transformed variables
    R_IM_II_trans_2.append(R_trans)
    t_IM_II_trans_2.append(t_trans)
# Transform the first solution
t_hat_IM_II_2 = np.array([symmetry_toolbox.IM_II_symmetry(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon,opt_para_IM_II_CML[1],0.044)[0] for index in range(len(t_sym))])
R_hat_IM_II_2 = np.array([symmetry_toolbox.IM_II_symmetry(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon,opt_para_IM_II_CML[1],0.044)[1] for index in range(len(t_sym))])        
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
axes[1].plot(t_sym, R_hat_IM_II_original, '-', color = (2/256,56/256,88/256), label='Original solution IM-II')
axes[1].plot(np.array(t_IM_II_trans_1[0]),np.array(R_IM_II_trans_1[0]),'--',color='black',label='Symmetry IM-II')
for index in range(1,len(t_IM_II_trans_1)):
    axes[1].plot(np.array(t_IM_II_trans_1[index]),np.array(R_IM_II_trans_1[index]),'--',color='black')
axes[1].plot(t_hat_IM_II_1, R_hat_IM_II_1, '-', color=(54/256,144/256,192/256), label='IM-II transformed solution 1')
for index in range(len(t_IM_II_trans_2)):
    axes[1].plot(np.array(t_IM_II_trans_2[index]),np.array(R_IM_II_trans_2[index]),'--',color='black')
axes[1].plot(t_hat_IM_II_2, R_hat_IM_II_2, '-', color=(208/256,209/256,230/256), label='IM-II transformed solution 2')
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


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustrate framework
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Allocate memory for a list
R_PLM_trans = []
t_PLM_trans = []
# Allocate an index vector
index_vector = list(np.linspace(40,len(t_CML)-1,30,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_CML[index],R_CML[index],epsilon-0.01,opt_para_PLM_CML[1],2)
    # Save the transformed variables
    R_PLM_trans.append(R_trans)
    t_PLM_trans.append(t_trans)
# Transform the data
t_CML_trans = np.array([symmetry_toolbox.PLM_2_symmetry(t_CML[index],R_CML[index],epsilon)[0] for index in range(len(t_CML))])
R_CML_trans = np.array([symmetry_toolbox.PLM_2_symmetry(t_CML[index],R_CML[index],epsilon)[1] for index in range(len(t_CML))])
# Fit the model to the transformed data
opt_para_PLM_CML, R_hat_PLM_CML_trans, RMS_PLM_CML = fit_to_data.PE_risk_profiles(t_CML_trans, R_CML_trans, "PLM","all")
# ----------------------------------------------------------------------------------
# Plot the illustration of the model selection framework
# ----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0].plot(t_CML, R_CML, '*', color='black', label='Data Myeloma cancer')
axes[0].plot(np.array(t_PLM_trans[0]),np.array(R_PLM_trans[0]),'--',color='black',label='Symmetry')
for index in range(1,len(t_PLM_trans)):
    axes[0].plot(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),'--',color='black')
axes[0].plot(t_CML_trans, R_CML_trans, '*', color='gray', label='Transformed Data')
axes[0].legend()
# Subplot 1a
axes[1].plot(t_CML_trans, R_CML_trans, '*', color='gray', label='Transformed Data')
axes[1].plot(t_CML_trans, R_hat_PLM_CML_trans, '-', color=(206/256,18/256,86/256), label='Fitted model')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The steps of the symmetry based model selection",fontsize=20, fontweight='bold')
plt.savefig("../Figures/illustrate_framework.png")
#plt.show()


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Validation plot
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Construct a t vector
t_sym = np.linspace(t_CML[0],t_CML[len(t_CML)-1],200)
# Create an epsilon vector
epsilon_vec = np.linspace(0,0.05,50)
# CREATE TWO RMS VECTORS FOR EACH MODEL
RMS_PLM = []
RMS_IM_II = []
# opt_para_PLM_CML
# THE ORIGINAL VECTORS OR ORIGINAL DATA IF YOU WILL
R_hat_PLM_original = np.array([fit_to_data.objective_PLM(t_sym[index],opt_para_PLM_CML[0],opt_para_PLM_CML[1]) for index in range(len(t_sym))])
R_hat_IM_II_original = np.array([fit_to_data.objective_IM_II(t_sym[index],opt_para_IM_II_CML[0],opt_para_IM_II_CML[1],0.044) for index in range(len(t_sym))])
# LOOP OVER THE EPSILON VECTORS AND FIT THE CANDIDATE MODELS TO THE TRANSFORMED TIME SERIES
for epsilon_temp in list(epsilon_vec):
    # TRANSFORM THE TIME ORIGINAL SOLUTION
    # PLM
    t_hat_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[0] for index in range(len(t_sym))])
    R_hat_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[1] for index in range(len(t_sym))])
    # IM-II
    t_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_temp,opt_para_IM_II_CML[1],0.044)[0] for index in range(len(t_sym))])
    R_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_temp,opt_para_IM_II_CML[1],0.044)[1] for index in range(len(t_sym))])
    # FIT MODELS TO TRANSFORMED SOLUTION
    RMS_PLM.append(fit_to_data.PE_risk_profiles(t_hat_PLM, R_hat_PLM, "PLM","all")[2])
    RMS_IM_II.append(fit_to_data.PE_risk_profiles(t_hat_IM_II, R_hat_IM_II, "IM-II","all")[2])
# ----------------------------------------------------------------------------------
# Plot the RMS vs epsilon
# ----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1
axes[0].plot(epsilon_vec,np.array(RMS_PLM),'-',color = (103/256,0/256,31/256),label='PLM')
axes[0].legend()
# Subplot 2
axes[1].plot(epsilon_vec,np.array(RMS_IM_II),'-',color = (2/256,56/256,88/256),label='IM_II')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("",fontsize=20, fontweight='bold')
plt.savefig("../Figures/illustrate_framework.png")

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Validation plot
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
print("opt_para_PLM_CML")
print(opt_para_PLM_CML)
print(np.log(opt_para_PLM_CML))
print("opt_para_IM_II")
print(opt_para_IM_II_CML)
print(np.log(opt_para_PLM_CML))
# CREATE TWO RMS VECTORS FOR EACH MODEL
RMS_PLM_log = []
RMS_IM_II_log = []
# opt_para_PLM_CML
# THE ORIGINAL VECTORS OR ORIGINAL DATA IF YOU WILL
R_hat_PLM_original_log = np.array([np.exp(fit_to_data.objective_PLM_log(t_sym[index],np.log(opt_para_PLM_CML[0]),opt_para_PLM_CML[1])) for index in range(len(t_sym))])
R_hat_IM_II_original_log = np.array([np.exp(fit_to_data.objective_IM_II_log(t_sym[index],np.log(opt_para_IM_II_CML[0]),opt_para_IM_II_CML[1],0.044)) for index in range(len(t_sym))])
# LOOP OVER THE EPSILON VECTORS AND FIT THE CANDIDATE MODELS TO THE TRANSFORMED TIME SERIES
for epsilon_temp in list(epsilon_vec):
    # TRANSFORM THE TIME ORIGINAL SOLUTION
    # PLM
    t_hat_PLM_log = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[0] for index in range(len(t_sym))])
    R_hat_PLM_log = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[1] for index in range(len(t_sym))])
    # IM-II
    t_hat_IM_II_log = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original_log[index],epsilon_temp,opt_para_IM_II_CML[1],0.044)[0] for index in range(len(t_sym))])
    R_hat_IM_II_log = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original_log[index],epsilon_temp,opt_para_IM_II_CML[1],0.044)[1] for index in range(len(t_sym))])
    # FIT MODELS TO TRANSFORMED SOLUTION
    RMS_PLM_log.append(fit_to_data.PE_risk_profiles_log(t_hat_PLM_log, R_hat_PLM_log, "PLM","all")[2])
    RMS_IM_II_log.append(fit_to_data.PE_risk_profiles_log(t_hat_IM_II_log, R_hat_IM_II_log, "IM-II","all")[2])

# ----------------------------------------------------------------------------------
# Plot the RMS vs epsilon
# ----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1
axes[0].plot(epsilon_vec,np.array(RMS_PLM_log),'-',color = (103/256,0/256,31/256),label='PLM')
axes[0].legend()
# Subplot 2
axes[1].plot(epsilon_vec,np.array(RMS_IM_II_log),'-',color = (2/256,56/256,88/256),label='IM_II')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("",fontsize=20, fontweight='bold')
plt.savefig("../Figures/illustrate_framework.png")
plt.show()




