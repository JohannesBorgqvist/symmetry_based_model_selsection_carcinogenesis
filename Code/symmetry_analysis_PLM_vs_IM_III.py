# =================================================================================
# =================================================================================
# Script:"symmetry_analysis_PLM_vs_IM_III.py"
# Date: 2022-02-25
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
PLM_fitted_to_myeloma_ODR, R_hat_PLM_myeloma_ODR, RMS_PLM_myeloma_ODR  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","ODR",[],[])
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_ODR, RMS_PLM_myeloma_ODR, "ODR")
# IM-III
IM_III_fitted_to_myeloma_ODR, R_hat_IM_III_myeloma_ODR, RMS_IM_III_myeloma_ODR = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-III","ODR",[],[])
write_output.save_data_PE("myeloma", "IM-III", IM_III_fitted_to_myeloma_ODR, RMS_IM_III_myeloma_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_ODR, R_hat_PLM_colon_ODR, RMS_PLM_colon_ODR  = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM","ODR",[],[])
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_ODR, RMS_PLM_colon_ODR, "ODR")
# IM-III
IM_III_fitted_to_colon_ODR, R_hat_IM_III_colon_ODR, RMS_IM_III_colon_ODR = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-III","ODR",[],[])
write_output.save_data_PE("colon", "IM-III", IM_III_fitted_to_colon_ODR, RMS_IM_III_colon_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_ODR, R_hat_PLM_CML_ODR, RMS_PLM_CML_ODR  = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM","ODR",[],[])
#write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_ODR, RMS_PLM_CML_ODR, "ODR")
# IM-III
IM_III_fitted_to_CML_ODR, R_hat_IM_III_CML_ODR, RMS_IM_III_CML_ODR = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-III","ODR",[],[])
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
#fig, axes = plt.subplots(1,1,figsize=(15,5))
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
plt.savefig("../Figures/Fit_of_models_to_cancer_data.png")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the data and the fit in latex
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT DATA
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Myeloma cancer")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_colon,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Colon cancer")
# CML
write_output.plot_LaTeX_2D(t_CML,R_CML,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Chronic Myeloid Leukemia (CML)")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Orthogonal Distance Regression (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT PLM
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_PLM_myeloma_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","color=pow_1,line width=2pt,","PLM")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_PLM_colon_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","color=pow_1,line width=2pt,","PLM")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_PLM_CML_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","color=pow_1,line width=2pt,","PLM")
# PLOT IM-III
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_IM_III_myeloma_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","color=mixed_1,line width=2pt,","IM-III")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_IM_III_colon_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","color=mixed_1,line width=2pt,","IM-III")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_IM_III_CML_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","color=mixed_1,line width=2pt,","IM-III")
#----------------------------------------------------------------------------------
# =================================================================================
# =================================================================================
# CALCULATE THE TRANSFORMATION SCALES
# =================================================================================
# What epsilon should we use to transform the age 85 years to 170 years? Here, comes
# the answer in the case of our two models.
# The PLM
epsilon_scale_PLM = symmetry_toolbox.PLM_transformation_scale(2)
# The IM-III
# Myeloma data
epsilon_scale_IM_III_myeloma = symmetry_toolbox.IM_III_transformation_scale(80,2,IM_III_fitted_to_myeloma_ODR.beta[3],IM_III_fitted_to_myeloma_ODR.beta[1])
# Colon data
epsilon_scale_IM_III_colon = symmetry_toolbox.IM_III_transformation_scale(80,2,IM_III_fitted_to_colon_ODR.beta[3],IM_III_fitted_to_colon_ODR.beta[1])
# CML data
epsilon_scale_IM_III_CML = symmetry_toolbox.IM_III_transformation_scale(80,2,IM_III_fitted_to_CML_ODR.beta[3],IM_III_fitted_to_CML_ODR.beta[1])
# Prompt to the user
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tThe transformation scales increasing the age from 80 years to 160 years\n")
print("\n\t--------------------------------------------------------------------------------------\n")
print("\t\tThe PLM:\tepsilon_PLM\t=\t%0.12f"%(epsilon_scale_PLM))
print("\t\tThe IM-III myeloma:\tepsilon_IM_III_myeloma\t=\t%0.12f"%(epsilon_scale_IM_III_myeloma))
print("\t\tThe IM-III colon:\tepsilon_IM_III_colon\t=\t%0.12f"%(epsilon_scale_IM_III_colon))
print("\t\tThe IM-III CML:\tepsilon_IM_III_CML\t=\t%0.12f"%(epsilon_scale_IM_III_CML))
epsilon_scale_IM_III_myeloma = symmetry_toolbox.IM_III_transformation_scale(85,2,IM_III_fitted_to_myeloma_ODR.beta[3],IM_III_fitted_to_myeloma_ODR.beta[1])
print("\t\tThe new IM-III myeloma scale for the plot:\tepsilon_IM_III_myeloma\t=\t%0.12f"%(epsilon_scale_IM_III_myeloma))
# =================================================================================
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
epsilon = epsilon_scale_PLM/2
# Allocate memory for a list
R_PLM_trans_1 = []
t_PLM_trans_1 = []
# Allocate an index vector
index_vector = list(np.arange(10,len(t_sym)-1,2))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon-0.01,PLM_fitted_to_myeloma_ODR.beta[1])
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
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon-0.01,PLM_fitted_to_myeloma_ODR.beta[1])
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2,R_hat_PLM_2 = symmetry_toolbox.PLM_transformed_solution(t_hat_PLM_1,R_hat_PLM_1,epsilon,PLM_fitted_to_myeloma_ODR.beta[0],PLM_fitted_to_myeloma_ODR.beta[1])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
epsilon = epsilon_scale_IM_III_myeloma/2
# THE IM-II
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
# Original solution
R_hat_IM_III_original = np.array([fit_to_data.objective_IM_III(IM_III_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_III_trans_1 = []
t_IM_III_trans_1 = []
# Allocate an index vector
index_vector = list(np.arange(193,len(t_sym)-1,2))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_sym[index],R_hat_IM_III_original[index],epsilon*0.85,IM_III_fitted_to_myeloma_ODR.beta[1],IM_III_fitted_to_myeloma_ODR.beta[3])
    # Save the transformed variables
    R_IM_III_trans_1.append(R_trans)
    t_IM_III_trans_1.append(t_trans)
# Transform the original solution
t_hat_IM_III_1,R_hat_IM_III_1 = symmetry_toolbox.IM_III_transformed_solution(t_sym,R_hat_IM_III_original,epsilon,IM_III_fitted_to_myeloma_ODR.beta[1],IM_III_fitted_to_myeloma_ODR.beta[3])
# Allocate memory for a list
R_IM_III_trans_2 = []
t_IM_III_trans_2 = []
# Update the index vector
index_vector = list(np.arange(190,len(t_sym)-1,2))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_hat_IM_III_1[index],R_hat_IM_III_1[index],epsilon*0.85,IM_III_fitted_to_myeloma_ODR.beta[1],IM_III_fitted_to_myeloma_ODR.beta[3])
    # Save the transformed variables
    R_IM_III_trans_2.append(R_trans)
    t_IM_III_trans_2.append(t_trans)
# Transform the second solution      
t_hat_IM_III_2,R_hat_IM_III_2 = symmetry_toolbox.IM_III_transformed_solution(t_hat_IM_III_1,R_hat_IM_III_1,epsilon,IM_III_fitted_to_myeloma_ODR.beta[1],IM_III_fitted_to_myeloma_ODR.beta[3])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the action of the symmetries in Python
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
#fig, axes = plt.subplots(1,1,figsize=(15,5))
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
plt.title("The symmetries of the PLM and the IM-III",fontsize=20, fontweight='bold')
plt.savefig("../Figures/action_of_symmetries.png")
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ILLUSTRATE THE ACTION OF THE SYMMETRIES IN LATEX
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
# PLM
#-----------------------------------------------------------------------------------------------------
write_output.plot_LaTeX_2D(t_sym,R_hat_PLM_original,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=pow_1,line width=2pt,","PLM sol., $R(t)$")
write_output.plot_LaTeX_2D(t_hat_PLM_1,R_hat_PLM_1,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=pow_2,line width=2pt,","PLM Transf. sol., $\hat{R}_1(t,\epsilon)$")
write_output.plot_LaTeX_2D(t_hat_PLM_2,R_hat_PLM_2,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=pow_3,line width=2pt,","PLM Transf. sol., $\hat{R}_2(t,\epsilon)$")
write_output.plot_LaTeX_2D(np.array(t_PLM_trans_1[0]),np.array(R_PLM_trans_1[0]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed","Symmetry $\Gamma^{\mathrm{PLM}}_{\epsilon}$")
for index in range(1,len(t_PLM_trans_1)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans_1[index]),np.array(R_PLM_trans_1[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])    
for index in range(len(t_PLM_trans_2)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans_2[index]),np.array(R_PLM_trans_2[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])    
#-----------------------------------------------------------------------------------------------------
# IM-II
#-----------------------------------------------------------------------------------------------------
write_output.plot_LaTeX_2D(t_sym,R_hat_IM_III_original,"../Figures/latex_figures/action_of_symmetries/Input/IM_III.tex","color=mixed_1,line width=2pt,","IM-III sol., $R(t)$")
write_output.plot_LaTeX_2D(t_hat_IM_III_1,R_hat_IM_III_1,"../Figures/latex_figures/action_of_symmetries/Input/IM_III.tex","color=mixed_2,line width=2pt,","IM-III Transf. sol., $\hat{R}_1(\hat{t}_1)$")
write_output.plot_LaTeX_2D(t_hat_IM_III_2,R_hat_IM_III_2,"../Figures/latex_figures/action_of_symmetries/Input/IM_III.tex","color=mixed_3,line width=2pt,","IM-III Transf. sol., $\hat{R}_2(\hat{t}_2)$")
write_output.plot_LaTeX_2D(np.array(t_IM_III_trans_1[0]),np.array(R_IM_III_trans_1[0]),"../Figures/latex_figures/action_of_symmetries/Input/IM_III.tex","color=black,->,>=latex,densely dashed","Symmetry $\Gamma^{\mathrm{IM-III}}_{\epsilon}$")
for index in range(1,len(t_IM_III_trans_1)):
    write_output.plot_LaTeX_2D(np.array(t_IM_III_trans_1[index]),np.array(R_IM_III_trans_1[index]),"../Figures/latex_figures/action_of_symmetries/Input/IM_III.tex","color=black,->,>=latex,densely dashed",[])  
for index in range(len(t_IM_III_trans_2)):
    write_output.plot_LaTeX_2D(np.array(t_IM_III_trans_2[index]),np.array(R_IM_III_trans_2[index]),"../Figures/latex_figures/action_of_symmetries/Input/IM_III.tex","color=black,->,>=latex,densely dashed",[])
# =================================================================================
# =================================================================================
# CONDUCT THE SYMMETRY BASED MODEL SELECTION
# =================================================================================
# =================================================================================
# Reset the epsilon scale for the PLM
# Myeloma data
epsilon_scale_IM_III_myeloma = symmetry_toolbox.IM_III_transformation_scale(80,2,IM_III_fitted_to_myeloma_ODR.beta[3],IM_III_fitted_to_myeloma_ODR.beta[1])
# Prompt to the user
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tThe symmetry based framework for model selection\n")
print("\n\t--------------------------------------------------------------------------------------\n")
# Allocate four epsilon vectors with transformation parameters
epsilon_vector_PLM = np.linspace(0.0,2*epsilon_scale_PLM,num=100,endpoint=True)
epsilon_vector_IM_III_myeloma = np.linspace(0.0,epsilon_scale_IM_III_myeloma,num=100,endpoint=True)
epsilon_vector_IM_III_colon = np.linspace(0.0,0.89*epsilon_scale_IM_III_colon,num=100,endpoint=True)
epsilon_vector_IM_III_CML = np.linspace(0.0,0.8286*epsilon_scale_IM_III_CML,num=100,endpoint=True)
# MYELOMA CANCER
print("\t\tModel\t=\t PLM,\tDataset\t=\t Myeloma")
RMS_transf_PLM_myeloma = symmetry_toolbox.symmetry_based_model_selection(t_myeloma,R_myeloma,epsilon_vector_PLM,PLM_fitted_to_myeloma_ODR.beta,"PLM")
print("\t\t\tDone!\n")
print("\t\tModel\t=\t IM-III,\tDataset\t=\t Myeloma")
RMS_transf_IM_III_myeloma = symmetry_toolbox.symmetry_based_model_selection(t_myeloma,R_myeloma,epsilon_vector_IM_III_myeloma,IM_III_fitted_to_myeloma_ODR.beta,"IM-III")
print("\t\t\tDone!\n")
# COLON CANCER
print("\t\tModel\t=\t PLM,\tDataset\t=\t Colon cancer")
RMS_transf_PLM_colon = symmetry_toolbox.symmetry_based_model_selection(t_colon,R_colon,epsilon_vector_PLM,PLM_fitted_to_colon_ODR.beta,"PLM")
print("\t\t\tDone!\n")
print("\t\tModel\t=\t IM-III,\tDataset\t=\t Colon cancer")
RMS_transf_IM_III_colon = symmetry_toolbox.symmetry_based_model_selection(t_colon,R_colon,epsilon_vector_IM_III_colon,IM_III_fitted_to_colon_ODR.beta,"IM-III")
print("\t\t\tDone!\n")
# CML CANCER
print("\t\tModel\t=\t PLM,\tDataset\t=\t CML")
RMS_transf_PLM_CML = symmetry_toolbox.symmetry_based_model_selection(t_CML,R_CML,epsilon_vector_PLM,PLM_fitted_to_CML_ODR.beta,"PLM")
print("\t\t\tDone!\n")
print("\t\tModel\t=\t IM-III,\tDataset\t=\t CML")
RMS_transf_IM_III_CML = symmetry_toolbox.symmetry_based_model_selection(t_CML,R_CML,epsilon_vector_IM_III_CML,IM_III_fitted_to_CML_ODR.beta,"IM-III")
print("\t\t\tDone!\n")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ILLUSTRATE THE SYMMETRY BASED MODEL SELECTION IN LATEX
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA
write_output.plot_LaTeX_2D(epsilon_vector_PLM,RMS_transf_PLM_myeloma,"../Figures/latex_figures/symmetry_based_model_selection/Input/PLM.tex","color=pow_1,line width=2pt,","PLM Myeloma cancer")
write_output.plot_LaTeX_2D(epsilon_vector_IM_III_myeloma,RMS_transf_IM_III_myeloma,"../Figures/latex_figures/symmetry_based_model_selection/Input/IM_III_myeloma.tex","color=mixed_1,line width=2pt,","IM-III Myeloma cancer")
# COLON CANCER
write_output.plot_LaTeX_2D(epsilon_vector_PLM,RMS_transf_PLM_colon,"../Figures/latex_figures/symmetry_based_model_selection/Input/PLM.tex","color=pow_2,line width=2pt,","PLM Colon cancer")
write_output.plot_LaTeX_2D(epsilon_vector_IM_III_colon,RMS_transf_IM_III_colon,"../Figures/latex_figures/symmetry_based_model_selection/Input/IM_III_colon.tex","color=mixed_2,line width=2pt,","IM-III Colon cancer")
# CML
write_output.plot_LaTeX_2D(epsilon_vector_PLM,RMS_transf_PLM_CML,"../Figures/latex_figures/symmetry_based_model_selection/Input/PLM.tex","color=pow_3,line width=2pt,","PLM CML")
write_output.plot_LaTeX_2D(epsilon_vector_IM_III_CML,RMS_transf_IM_III_CML,"../Figures/latex_figures/symmetry_based_model_selection/Input/IM_III_colon.tex","color=mixed_3,line width=2pt,","IM-III CML")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the symmetry based model selection
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLM RESULTS
# Overall properties
fig, axes = plt.subplots(1,1,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
axes.plot(epsilon_vector_PLM,RMS_transf_PLM_myeloma,'-', color = (103/256,0/256,31/256),label='PLM Myeloma cancer')
axes.plot(epsilon_vector_PLM,RMS_transf_PLM_colon,'-', color = (206/256,18/256,86/256),label='PLM Colon cancer')
axes.plot(epsilon_vector_PLM,RMS_transf_PLM_CML,'-', color = (223/256,101/256,176/256),label='PLM CML')
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("The symmetry based model selection for the PLM",fontsize=20, fontweight='bold')
plt.savefig("../Figures/symmetry_based_model_selection_PLM.png")
# IM-III
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: Myeloma
axes[0].plot(epsilon_vector_IM_III_myeloma,RMS_transf_IM_III_myeloma,'-', color = (2/256,56/256,88/256),label='IM-III Myeloma cancer')
# Subplot 2: Colon cancer
axes[1].plot(epsilon_vector_IM_III_colon,RMS_transf_IM_III_colon,'-', color = (54/256,144/256,192/256),label='IM-III Colon cancer')
# Subplot 3: CML
axes[2].plot(epsilon_vector_IM_III_CML,RMS_transf_IM_III_CML,'-', color = (208/256,209/256,230/256),label='IM-III CML')
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("The symmetry based model selection for the IM-III",fontsize=20, fontweight='bold')
plt.savefig("../Figures/symmetry_based_model_selection_IM_III.png")
plt.show()    
# Prompt to the user
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tCalculations are done!\n")
print("\n\t--------------------------------------------------------------------------------------\n")

    
