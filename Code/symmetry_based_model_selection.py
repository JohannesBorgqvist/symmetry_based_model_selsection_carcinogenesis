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
# Fit the models to the data
# =================================================================================
# =================================================================================
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# STANDARD LEAST SQUARE (LS) FITTING
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_LS, R_hat_PLM_myeloma_LS, RMS_PLM_myeloma_LS  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","LS",[])
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_LS, RMS_PLM_myeloma_LS, "LS")
# IM-III
IM_III_fitted_to_myeloma_LS, R_hat_IM_III_myeloma_LS, RMS_IM_III_myeloma_LS = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-III","LS",[])
write_output.save_data_PE("myeloma", "IM-III", IM_III_fitted_to_myeloma_LS, RMS_IM_III_myeloma_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_LS, R_hat_PLM_colon_LS, RMS_PLM_colon_LS  = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM","LS",[])
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_LS, RMS_PLM_colon_LS, "LS")
# IM-III
IM_III_fitted_to_colon_LS, R_hat_IM_III_colon_LS, RMS_IM_III_colon_LS = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-III","LS",[])
write_output.save_data_PE("colon", "IM-III", IM_III_fitted_to_colon_LS, RMS_IM_III_colon_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_LS, R_hat_PLM_CML_LS, RMS_PLM_CML_LS  = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM","LS",[])
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_LS, RMS_PLM_CML_LS, "LS")
# IM-III
IM_III_fitted_to_CML_LS, R_hat_IM_III_CML_LS, RMS_IM_III_CML_LS = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-III","LS",[])
write_output.save_data_PE("CML", "IM-III", IM_III_fitted_to_CML_LS, RMS_IM_III_CML_LS, "LS")

#----------------------------------------------------------------------------------
# =================================================================================
# =================================================================================
# Plot the data and the fit in Python
# =================================================================================
# =================================================================================
# PLOT OF THE FIT OF THE MODEL TO THE DATA
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0].plot(t_myeloma, R_hat_PLM_myeloma_LS, '-', color = (103/256,0/256,31/256),label='LS fit PLM')
axes[0].plot(t_myeloma, R_hat_IM_III_myeloma_LS, '-', color = (2/256,56/256,88/256),label='LS fit IM-III')
axes[0].legend()
# Subplot 2a
axes[1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1].plot(t_colon, R_hat_PLM_colon_LS, '-', color = (103/256,0/256,31/256),label='LS fit PLM')
axes[1].plot(t_colon, R_hat_IM_III_colon_LS, '-', color = (2/256,56/256,88/256),label='LS fit IM-III')
axes[1].legend()
# Subplot 3a
axes[2].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2].plot(t_CML, R_hat_PLM_CML_LS, '-', color = (103/256,0/256,31/256),label='LS fit PLM')
axes[2].plot(t_CML, R_hat_IM_III_CML_LS, '-', color = (2/256,56/256,88/256),label='LS fit IM-III')
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
# =================================================================================
# =================================================================================
# Plot the data and the fit in latex
# =================================================================================
# =================================================================================
# PLOT DATA
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Myeloma cancer")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_colon,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Colon cancer")
# CML
write_output.plot_LaTeX_2D(t_CML,R_CML,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Chronic Myeloid Leukemia (CML)")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# STANDARD LEAST SQUARE (LS) FITTING
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT PLM
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_PLM_myeloma_LS,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","color=pow_1,line width=2pt,","PLM")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_PLM_colon_LS,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","color=pow_1,line width=2pt,","PLM")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_PLM_CML_LS,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","color=pow_1,line width=2pt,","PLM")
# PLOT IM-III
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_IM_III_myeloma_LS,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","color=mixed_1,line width=2pt,","IM-III")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_IM_III_colon_LS,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","color=mixed_1,line width=2pt,","IM-III")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_IM_III_CML_LS,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","color=mixed_1,line width=2pt,","IM-III")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustrate symmetries
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE PLM
R_hat_PLM_original = np.array([fit_to_data.objective_PLM(PLM_fitted_to_myeloma_LS.beta,t_sym[index]) for index in range(len(t_sym))])
# Choose an epsilon
epsilon = 0.45
# Allocate memory for a list
R_PLM_trans_1 = []
t_PLM_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(10,len(t_sym)-1,75,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon-0.045,PLM_fitted_to_myeloma_LS.beta[1])
    # Save the transformed variables
    R_PLM_trans_1.append(R_trans)
    t_PLM_trans_1.append(t_trans)
# Transform the original solution
t_hat_PLM_1,R_hat_PLM_1 = symmetry_toolbox.PLM_transformed_solution(t_sym,R_hat_PLM_original,epsilon,PLM_fitted_to_myeloma_LS.beta[0],PLM_fitted_to_myeloma_LS.beta[1])
# Allocate memory for a list
R_PLM_trans_2 = []
t_PLM_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon-0.045,PLM_fitted_to_myeloma_LS.beta[1])
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2,R_hat_PLM_2 = symmetry_toolbox.PLM_transformed_solution(t_hat_PLM_1,R_hat_PLM_1,epsilon,PLM_fitted_to_myeloma_LS.beta[0],PLM_fitted_to_myeloma_LS.beta[1])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE IM-II
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
# Original solution
R_hat_IM_III_original = np.array([fit_to_data.objective_IM_III(IM_III_fitted_to_myeloma_LS.beta,t_sym[index]) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_III_trans_1 = []
t_IM_III_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(175,len(t_sym)-1,25,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_sym[index],R_hat_IM_III_original[index],epsilon-0.045,IM_III_fitted_to_myeloma_LS.beta[1],0.044)
    # Save the transformed variables
    R_IM_III_trans_1.append(R_trans)
    t_IM_III_trans_1.append(t_trans)
# Transform the original solution
t_hat_IM_III_1,R_hat_IM_III_1 = symmetry_toolbox.IM_III_transformed_solution(t_sym,R_hat_IM_III_original,epsilon,IM_III_fitted_to_myeloma_LS.beta[0],IM_III_fitted_to_myeloma_LS.beta[1],IM_III_fitted_to_myeloma_LS.beta[2])
# Allocate memory for a list
R_IM_III_trans_2 = []
t_IM_III_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_III_transformation(t_hat_IM_III_1[index],R_hat_IM_III_1[index],epsilon-0.045,IM_III_fitted_to_myeloma_LS.beta[1],0.044)
    # Save the transformed variables
    R_IM_III_trans_2.append(R_trans)
    t_IM_III_trans_2.append(t_trans)
# Transform the second solution      
t_hat_IM_III_2,R_hat_IM_III_2 = symmetry_toolbox.IM_III_transformed_solution(t_hat_IM_III_1,R_hat_IM_III_1,epsilon,IM_III_fitted_to_myeloma_LS.beta[0],IM_III_fitted_to_myeloma_LS.beta[1],IM_III_fitted_to_myeloma_LS.beta[2])
# =================================================================================
# =================================================================================
# Plot the illustration of the model selection framework
# =================================================================================
# =================================================================================
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
plt.show()
plt.savefig("../Figures/action_of_symmetries.png")
# =================================================================================
# =================================================================================
# ILLUSTRATE THE ACTION OF THE SYMMETRIES IN LATEX
# =================================================================================
# =================================================================================
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
