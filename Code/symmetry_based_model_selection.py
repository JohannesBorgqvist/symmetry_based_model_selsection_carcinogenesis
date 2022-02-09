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
import fit_to_data_ODR # Home-made
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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# STANDARD LEAST SQUARE (LS) FITTING
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_LS, R_hat_PLM_myeloma_LS, RMS_PLM_myeloma_LS  = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","LS",[])
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_LS, RMS_PLM_myeloma_LS, "LS")
# IM-II
IM_II_fitted_to_myeloma_LS, R_hat_IM_II_myeloma_LS, RMS_IM_II_myeloma_LS = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","LS",[])
write_output.save_data_PE("myeloma", "IM-III", IM_II_fitted_to_myeloma_LS, RMS_IM_II_myeloma_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_LS, R_hat_PLM_colon_LS, RMS_PLM_colon_LS  = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"PLM","LS",[])
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_LS, RMS_PLM_colon_LS, "LS")
# IM-II
IM_II_fitted_to_colon_LS, R_hat_IM_II_colon_LS, RMS_IM_II_colon_LS = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"IM-II","LS",[])
write_output.save_data_PE("colon", "IM-III", IM_II_fitted_to_colon_LS, RMS_IM_II_colon_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_LS, R_hat_PLM_CML_LS, RMS_PLM_CML_LS  = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"PLM","LS",[])
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_LS, RMS_PLM_CML_LS, "LS")
# IM-II
IM_II_fitted_to_CML_LS, R_hat_IM_II_CML_LS, RMS_IM_II_CML_LS = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"IM-II","LS",[])
write_output.save_data_PE("CML", "IM-III", IM_II_fitted_to_CML_LS, RMS_IM_II_CML_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ORTHOGONAL DISTANCE REGRESSION (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_ODR, R_hat_PLM_myeloma_ODR, RMS_PLM_myeloma_ODR  = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","ODR",[])
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_ODR, RMS_PLM_myeloma_ODR, "ODR")
# IM-II
IM_II_fitted_to_myeloma_ODR, R_hat_IM_II_myeloma_ODR, RMS_IM_II_myeloma_ODR = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","ODR",[])
write_output.save_data_PE("myeloma", "IM-III", IM_II_fitted_to_myeloma_ODR, RMS_IM_II_myeloma_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_ODR, R_hat_PLM_colon_ODR, RMS_PLM_colon_ODR  = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"PLM","ODR",[])
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_ODR, RMS_PLM_colon_ODR, "ODR")
# IM-II
IM_II_fitted_to_colon_ODR, R_hat_IM_II_colon_ODR, RMS_IM_II_colon_ODR = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"IM-II","ODR",[])
write_output.save_data_PE("colon", "IM-III", IM_II_fitted_to_colon_ODR, RMS_IM_II_colon_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_ODR, R_hat_PLM_CML_ODR, RMS_PLM_CML_ODR  = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"PLM","ODR",[])
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_ODR, RMS_PLM_CML_ODR, "ODR")
# IM-II
IM_II_fitted_to_CML_ODR, R_hat_IM_II_CML_ODR, RMS_IM_II_CML_ODR = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"IM-II","ODR",[])
write_output.save_data_PE("CML", "IM-III", IM_II_fitted_to_CML_ODR, RMS_IM_II_CML_ODR, "ODR")
#----------------------------------------------------------------------------------
# =================================================================================
# =================================================================================
# Plot the data and the fit in Python
# =================================================================================
# =================================================================================
# PLOT OF THE FIT OF THE MODEL TO THE DATA
# Overall properties
fig, axes = plt.subplots(3,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0][0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0][0].plot(t_myeloma, R_hat_PLM_myeloma_LS, '-', color = (103/256,0/256,31/256),label='LS fit PLM')
axes[0][0].plot(t_myeloma, R_hat_IM_II_myeloma_LS, '-', color = (2/256,56/256,88/256),label='LS fit IM-II')
axes[0][0].legend()
# Subplot 1b
axes[0][1].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0][1].plot(t_myeloma, R_hat_PLM_myeloma_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[0][1].plot(t_myeloma, R_hat_IM_II_myeloma_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM-II')
axes[0][1].legend()
# Subplot 2a
axes[1][0].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1][0].plot(t_colon, R_hat_PLM_colon_LS, '-', color = (103/256,0/256,31/256),label='LS fit PLM')
axes[1][0].plot(t_colon, R_hat_IM_II_colon_LS, '-', color = (2/256,56/256,88/256),label='LS fit IM-II')
axes[1][0].legend()
# Subplot 2b
axes[1][1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1][1].plot(t_colon, R_hat_PLM_colon_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[1][1].plot(t_colon, R_hat_IM_II_colon_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM-II')
axes[1][1].legend()
# Subplot 3a
axes[2][0].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2][0].plot(t_CML, R_hat_PLM_CML_LS, '-', color = (103/256,0/256,31/256),label='LS fit PLM')
axes[2][0].plot(t_CML, R_hat_IM_II_CML_LS, '-', color = (2/256,56/256,88/256),label='LS fit IM-II')
axes[2][0].legend()
# Subplot 3b
axes[2][1].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2][1].plot(t_CML, R_hat_PLM_CML_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[2][1].plot(t_CML, R_hat_IM_II_CML_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM-II')
axes[2][1].legend()
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
write_output.plot_LaTeX_2D(t_myeloma,R_hat_PLM_myeloma_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","color=pow_1,line width=2pt,","PLM")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_PLM_colon_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","color=pow_1,line width=2pt,","PLM")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_PLM_CML_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","color=pow_1,line width=2pt,","PLM")
# PLOT IM-II
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_IM_II_myeloma_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/myeloma.tex","color=mixed_1,line width=2pt,","IM-III")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_IM_II_colon_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/colon.tex","color=mixed_1,line width=2pt,","IM-III")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_IM_II_CML_ODR,"../Figures/latex_figures/fit_of_models_to_cancer_data/Input/CML.tex","color=mixed_1,line width=2pt,","IM-III")

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
R_hat_PLM_original = np.array([fit_to_data_ODR.objective_PLM(PLM_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Choose an epsilon
epsilon = 0.32
# Allocate memory for a list
R_PLM_trans_1 = []
t_PLM_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(10,len(t_sym)-1,75,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon-0.01,PLM_fitted_to_myeloma_ODR.beta[1],2)
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
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon-0.01,PLM_fitted_to_myeloma_ODR.beta[1],2)
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2 = np.array([symmetry_toolbox.PLM_2_symmetry(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon)[0] for index in range(len(t_sym))])
R_hat_PLM_2 = np.array([symmetry_toolbox.PLM_2_symmetry(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon)[1] for index in range(len(t_sym))])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE IM-II
R_hat_IM_II_original = np.array([fit_to_data_ODR.objective_IM_II(IM_II_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_II_trans_1 = []
t_IM_II_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(175,len(t_sym)-1,25,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_sym[index],R_hat_IM_II_original[index],epsilon-0.01,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)
    # Save the transformed variables
    R_IM_II_trans_1.append(R_trans)
    t_IM_II_trans_1.append(t_trans)
# Transform the original solution
t_hat_IM_II_1 = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[0] for index in range(len(t_sym))])
R_hat_IM_II_1 = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[1] for index in range(len(t_sym))])    
# Allocate memory for a list
R_IM_II_trans_2 = []
t_IM_II_trans_2 = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon-0.01,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)
    # Save the transformed variables
    R_IM_II_trans_2.append(R_trans)
    t_IM_II_trans_2.append(t_trans)
# Transform the first solution
t_hat_IM_II_2 = np.array([symmetry_toolbox.IM_II_symmetry(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[0] for index in range(len(t_sym))])
R_hat_IM_II_2 = np.array([symmetry_toolbox.IM_II_symmetry(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[1] for index in range(len(t_sym))])        
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
# =================================================================================
# =================================================================================
# ILLUSTRATE THE ACTION OF THE SYMMETRIES IN LATEX
# =================================================================================
# =================================================================================
#-----------------------------------------------------------------------------------------------------
# PLM
#-----------------------------------------------------------------------------------------------------
write_output.plot_LaTeX_2D(t_sym,R_hat_PLM_original,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=pow_1,line width=2pt,","PLM Original solution, $R(t)$")
write_output.plot_LaTeX_2D(t_hat_PLM_1,R_hat_PLM_1,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=pow_2,line width=2pt,","PLM Transformed solution, $\hat{R}_1(\hat{t}_1)$")
write_output.plot_LaTeX_2D(t_hat_PLM_2,R_hat_PLM_2,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=pow_3,line width=2pt,","PLM Transformed solution, $\hat{R}_2(\hat{t}_2)$")
write_output.plot_LaTeX_2D(np.array(t_PLM_trans_1[0]),np.array(R_PLM_trans_1[0]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed","Action of symmetry $\Gamma_{1,1}(\epsilon)$")
for index in range(1,len(t_PLM_trans_1)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans_1[index]),np.array(R_PLM_trans_1[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])    
for index in range(len(t_PLM_trans_2)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans_2[index]),np.array(R_PLM_trans_2[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])    
#-----------------------------------------------------------------------------------------------------
# IM-II
#-----------------------------------------------------------------------------------------------------
write_output.plot_LaTeX_2D(t_sym,R_hat_IM_II_original,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=mixed_1,line width=2pt,","IM-II Original solution, $R(t)$")
write_output.plot_LaTeX_2D(t_hat_IM_II_1,R_hat_PLM_1,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=mixed_2,line width=2pt,","IM-II Transformed solution, $\hat{R}_1(\hat{t}_1)$")
write_output.plot_LaTeX_2D(t_hat_IM_II_2,R_hat_PLM_2,"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=mixed_3,line width=2pt,","IM-II Transformed solution, $\hat{R}_2(\hat{t}_2)$")
write_output.plot_LaTeX_2D(np.array(t_IM_II_trans_1[0]),np.array(R_IM_II_trans_1[0]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed","Action of symmetry $\Gamma_{2,1}(\epsilon)$")
for index in range(1,len(t_IM_II_trans_1)):
    write_output.plot_LaTeX_2D(np.array(t_IM_II_trans_1[index]),np.array(R_IM_II_trans_1[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])  
for index in range(len(t_IM_II_trans_2)):
    write_output.plot_LaTeX_2D(np.array(t_IM_II_trans_2[index]),np.array(R_IM_II_trans_2[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])  
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustrate framework
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Pick a value of epsilon
epsilon = 0.30
# PLM
# Allocate memory for a list
R_PLM_trans = []
t_PLM_trans = []
# Allocate an index vector
index_vector = list(np.linspace(10,55,30,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_myeloma[index],R_myeloma[index],epsilon-0.01,PLM_fitted_to_myeloma_LS.beta[1],2)
    # Save the transformed variables
    R_PLM_trans.append(R_trans)
    t_PLM_trans.append(t_trans)
# Transform the data
t_myeloma_trans_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon)[0] for index in range(len(t_myeloma))])
R_myeloma_trans_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon)[1] for index in range(len(t_myeloma))])
# Fit the model to the transformed data
param_temp, R_hat_PLM_trans, RMS_temp  = fit_to_data_ODR.PE_risk_profiles(t_myeloma_trans_PLM,R_myeloma_trans_PLM,"PLM","LS",[])
# IM-III
# Allocate memory for a list
R_IM_III_trans = []
t_IM_III_trans = []
# Allocate an index vector
index_vector = list(np.linspace(10,55,30,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_myeloma[index],R_myeloma[index],epsilon-0.01,IM_II_fitted_to_myeloma_LS.beta[1],0.044)
    # Save the transformed variables
    R_IM_III_trans.append(R_trans)
    t_IM_III_trans.append(t_trans)
# Transform the data
t_myeloma_trans_IM_III = np.array([symmetry_toolbox.IM_II_symmetry(t_myeloma[index],R_myeloma[index],epsilon,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[0] for index in range(len(t_myeloma))])
R_myeloma_trans_IM_III = np.array([symmetry_toolbox.IM_II_symmetry(t_myeloma[index],R_myeloma[index],epsilon,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[1] for index in range(len(t_myeloma))])
# Fit the model to the transformed data
param_temp, R_hat_IM_III_trans, RMS_temp  = fit_to_data_ODR.PE_risk_profiles(t_myeloma_trans_IM_III,R_myeloma_trans_IM_III,"IM-II","LS",[])
# =================================================================================
# =================================================================================
# Plot the illustration of the model selection framework
# =================================================================================
# =================================================================================
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
axes[1][0].plot(t_myeloma_trans_PLM, R_myeloma_trans_PLM, '*', color='gray', label='Transformed Data')
axes[1][0].plot(t_myeloma_trans_PLM, R_hat_PLM_trans, '-', color=(206/256,18/256,86/256), label='Fitted PLM')
axes[1][0].legend()
# Subplot 2a
axes[0][1].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0][1].plot(np.array(t_IM_III_trans[0]),np.array(R_IM_III_trans[0]),'--',color='black',label='Symmetry IM-III')
for index in range(1,len(t_IM_III_trans)):
    axes[0][1].plot(np.array(t_IM_III_trans[index]),np.array(R_IM_III_trans[index]),'--',color='black')
axes[0][1].plot(t_myeloma_trans_IM_III, R_myeloma_trans_IM_III, '*', color='gray', label='Transformed Data')
axes[0][1].legend()
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
#plt.show()
# =================================================================================
# =================================================================================
# ILLUSTRATE THE FRAMEWORK IN LATEX
# =================================================================================
# =================================================================================
# Myeloma data
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/illustrate_data/Input/data.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data")
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/illustrate_data/Input/PLM.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data")
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/illustrate_data/Input/IM_II.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data")
#-----------------------------------------------------------------------------------------------------
# PLM
#-----------------------------------------------------------------------------------------------------
write_output.plot_LaTeX_2D(np.array(t_PLM_trans[0]),np.array(R_PLM_trans[0]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed","Action of symmetry $\Gamma_{1,1}(\epsilon)$")
axes[0][0].plot(,'--',color='black',label='Symmetry PLM')
for index in range(1,len(t_PLM_trans)):
    axes[0][0].plot(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),'--',color='black')
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),"../Figures/latex_figures/action_of_symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])

axes[0][0].plot(t_myeloma_trans_PLM, R_myeloma_trans_PLM, '*', color='gray', label='Transformed Data')
axes[1][0].plot(t_myeloma_trans_PLM, R_hat_PLM_trans, '-', color=(206/256,18/256,86/256), label='Fitted PLM')
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Validation plot
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Construct a t vector
t_sym = np.linspace(t_myeloma[0],t_myeloma[len(t_myeloma)-1],200)
# Create an epsilon vector
epsilon_vec = np.linspace(0,0.2,200)
# CREATE TWO RMS VECTORS FOR EACH MODEL
RMS_PLM = []
RMS_IM_II = []
# THE ORIGINAL VECTORS OR ORIGINAL DATA IF YOU WILL
R_hat_PLM_original = np.array([fit_to_data_ODR.objective_PLM(PLM_fitted_to_myeloma_LS.beta,t_sym[index]) for index in range(len(t_sym))])
R_hat_IM_II_original = np.array([fit_to_data_ODR.objective_IM_II(IM_II_fitted_to_myeloma_LS.beta,t_sym[index]) for index in range(len(t_sym))])
# LOOP OVER THE EPSILON VECTORS AND FIT THE CANDIDATE MODELS TO THE TRANSFORMED TIME SERIES
for epsilon_temp in list(epsilon_vec):
    # TRANSFORM THE TIME ORIGINAL SOLUTION
    # PLM
    t_hat_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[0] for index in range(len(t_sym))])
    R_hat_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[1] for index in range(len(t_sym))])
    # IM-II
    t_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_temp,IM_II_fitted_to_myeloma_LS.beta[1],0.044)[0] for index in range(len(t_sym))])
    R_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_temp,IM_II_fitted_to_myeloma_LS.beta[1],0.044)[1] for index in range(len(t_sym))])
    # FIT MODELS TO TRANSFORMED SOLUTION
    RMS_PLM.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM,R_hat_PLM,"PLM","LS",[])[2])
    RMS_IM_II.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II,R_hat_IM_II,"IM-II","LS",[])[2])
# Look at a single transformed solution
epsilon_fixed = 0.12    
# Transform time series
t_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_fixed,IM_II_fitted_to_myeloma_LS.beta[1],0.044)[0] for index in range(len(t_sym))])
R_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_fixed,IM_II_fitted_to_myeloma_LS.beta[1],0.044)[1] for index in range(len(t_sym))])
# Fit to transformed time series
IM_II_hat_fit, R_hat_IM_II_fit, RMS_hat_IM_II_fit = fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II,R_hat_IM_II,"IM-II","LS",[])
# =================================================================================
# =================================================================================
# Plot the validation plot
# =================================================================================
# =================================================================================
# Validation plot
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
plt.title("Numerical validation of the implemented symmetries",fontsize=20, fontweight='bold')
plt.savefig("../Figures/validation_plot.png")
#plt.show()
# Illustration of the framework
# Overall properties
#fig, axes = plt.subplots(1,1,figsize=(15,5))
plt.figure()
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.plot(t_sym,R_hat_IM_II,'*',color = 'black',label='IM-II transformed solution')
plt.plot(t_sym,R_hat_IM_II_fit,'-',color = (103/256,0/256,31/256),label='IM-II fitted curve')
plt.legend()
# add a big axis, hide frame
#fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidences of cancer, $R(t)$")
# displaying the title
plt.title("Fit to transformed solution with $\epsilon=%0.3f$: $\mathrm{RMS}(\epsilon)=%0.3f$"%(epsilon_fixed,RMS_hat_IM_II_fit),fontsize=20, fontweight='bold')
plt.savefig("../Figures/validation_fit_IM_II_individual_curve.png")
#plt.show()
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Symmetry based model selection LS
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Update the epsilon vector
epsilon_vec = np.linspace(0,0.6,200)
#epsilon_vec = np.linspace(0,0.4,200)
# CREATE TWO RMS VECTORS FOR EACH MODEL AND EACH DATA SET WHICH GIVES SIX VECTORS
# IN TOTAL
RMS_PLM_myeloma = []
RMS_IM_II_myeloma = []
RMS_PLM_colon = []
RMS_IM_II_colon = []
RMS_PLM_CML = []
RMS_IM_II_CML = []
# LOOP OVER THE EPSILON VECTORS AND FIT THE CANDIDATE MODELS TO THE TRANSFORMED TIME SERIES
for epsilon_temp in list(epsilon_vec):
    # TRANSFORM THE MYELOMA DATA
    # PLM
    t_hat_PLM_myeloma = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp)[0] for index in range(len(t_myeloma))])
    R_hat_PLM_myeloma = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp)[1] for index in range(len(t_myeloma))])
    # IM-II
    t_hat_IM_II_myeloma = np.array([symmetry_toolbox.IM_II_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp,30,0.044)[0] for index in range(len(t_myeloma))])
    R_hat_IM_II_myeloma = np.array([symmetry_toolbox.IM_II_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp,30,0.044)[1] for index in range(len(t_myeloma))])
    # TRANSFORM THE COLON DATA
    # PLM
    t_hat_PLM_colon = np.array([symmetry_toolbox.PLM_2_symmetry(t_colon[index],R_colon[index],epsilon_temp)[0] for index in range(len(t_colon))])
    R_hat_PLM_colon = np.array([symmetry_toolbox.PLM_2_symmetry(t_colon[index],R_colon[index],epsilon_temp)[1] for index in range(len(t_colon))])
    # IM-II
    t_hat_IM_II_colon = np.array([symmetry_toolbox.IM_II_symmetry(t_colon[index],R_colon[index],epsilon_temp,30,0.044)[0] for index in range(len(t_colon))])
    R_hat_IM_II_colon = np.array([symmetry_toolbox.IM_II_symmetry(t_colon[index],R_colon[index],epsilon_temp,30,0.044)[1] for index in range(len(t_colon))])
    # TRANSFORM THE CML DATA
    # PLM
    t_hat_PLM_CML = np.array([symmetry_toolbox.PLM_2_symmetry(t_CML[index],R_CML[index],epsilon_temp)[0] for index in range(len(t_CML))])
    R_hat_PLM_CML = np.array([symmetry_toolbox.PLM_2_symmetry(t_CML[index],R_CML[index],epsilon_temp)[1] for index in range(len(t_CML))])
    # IM-II
    t_hat_IM_II_CML = np.array([symmetry_toolbox.IM_II_symmetry(t_CML[index],R_CML[index],epsilon_temp,30,0.044)[0] for index in range(len(t_CML))])
    R_hat_IM_II_CML = np.array([symmetry_toolbox.IM_II_symmetry(t_CML[index],R_CML[index],epsilon_temp,30,0.044)[1] for index in range(len(t_CML))])    
    # FIT MODELS TO TRANSFORMED DATA
    RMS_PLM_myeloma.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM_myeloma,R_hat_PLM_myeloma,"PLM","LS",[])[2])
    RMS_IM_II_myeloma.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II_myeloma,R_hat_IM_II_myeloma,"IM-II","LS",[])[2])
    RMS_PLM_colon.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM_colon,R_hat_PLM_colon,"PLM","LS",[])[2])
    RMS_IM_II_colon.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II_colon,R_hat_IM_II_colon,"IM-II","LS",[])[2])    
    RMS_PLM_CML.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM_CML,R_hat_PLM_CML,"PLM","LS",[])[2])
    RMS_IM_II_CML.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II_CML,R_hat_IM_II_CML,"IM-II","LS",[])[2])    
# =================================================================================
# =================================================================================
# Symmetry based model selection LS
# =================================================================================
# =================================================================================
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(3,1,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: Myeloma
axes[0].plot(epsilon_vec,np.array(RMS_PLM_myeloma),'-',color = (103/256,0/256,31/256),label='PLM Myeloma')
axes[0].plot(epsilon_vec,np.array(RMS_IM_II_myeloma),'-',color = (2/256,56/256,88/256),label='IM_II Myeloma')
axes[0].legend()
# Subplot 2: Colon
axes[1].plot(epsilon_vec,np.array(RMS_PLM_colon),'-',color = (103/256,0/256,31/256),label='PLM Colon')
axes[1].plot(epsilon_vec,np.array(RMS_IM_II_colon),'-',color = (2/256,56/256,88/256),label='IM_II Colon')
axes[1].legend()
# Subplot 3: CML
axes[2].plot(epsilon_vec,np.array(RMS_PLM_CML),'-',color = (103/256,0/256,31/256),label='PLM CML')
axes[2].plot(epsilon_vec,np.array(RMS_IM_II_CML),'-',color = (2/256,56/256,88/256),label='IM_II CML')
axes[2].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("Symmetry based model selection using LS",fontsize=20, fontweight='bold')
plt.savefig("../Figures/symmetry_based_model_selection_LS.png")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Symmetry based model selection ODR
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CREATE TWO RMS VECTORS FOR EACH MODEL AND EACH DATA SET WHICH GIVES SIX VECTORS
# IN TOTAL
RMS_PLM_myeloma_ODR = []
RMS_IM_II_myeloma_ODR = []
RMS_PLM_colon_ODR = []
RMS_IM_II_colon_ODR = []
RMS_PLM_CML_ODR = []
RMS_IM_II_CML_ODR = []
# LOOP OVER THE EPSILON VECTORS AND FIT THE CANDIDATE MODELS TO THE TRANSFORMED TIME SERIES
for epsilon_temp in list(epsilon_vec):
    # TRANSFORM THE MYELOMA DATA
    # PLM
    t_hat_PLM_myeloma = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp)[0] for index in range(len(t_myeloma))])
    R_hat_PLM_myeloma = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp)[1] for index in range(len(t_myeloma))])
    # IM-II
    t_hat_IM_II_myeloma = np.array([symmetry_toolbox.IM_II_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp,30,0.044)[0] for index in range(len(t_myeloma))])
    R_hat_IM_II_myeloma = np.array([symmetry_toolbox.IM_II_symmetry(t_myeloma[index],R_myeloma[index],epsilon_temp,30,0.044)[1] for index in range(len(t_myeloma))])
    # TRANSFORM THE COLON DATA
    # PLM
    t_hat_PLM_colon = np.array([symmetry_toolbox.PLM_2_symmetry(t_colon[index],R_colon[index],epsilon_temp)[0] for index in range(len(t_colon))])
    R_hat_PLM_colon = np.array([symmetry_toolbox.PLM_2_symmetry(t_colon[index],R_colon[index],epsilon_temp)[1] for index in range(len(t_colon))])
    # IM-II
    t_hat_IM_II_colon = np.array([symmetry_toolbox.IM_II_symmetry(t_colon[index],R_colon[index],epsilon_temp,30,0.044)[0] for index in range(len(t_colon))])
    R_hat_IM_II_colon = np.array([symmetry_toolbox.IM_II_symmetry(t_colon[index],R_colon[index],epsilon_temp,30,0.044)[1] for index in range(len(t_colon))])
    # TRANSFORM THE CML DATA
    # PLM
    t_hat_PLM_CML = np.array([symmetry_toolbox.PLM_2_symmetry(t_CML[index],R_CML[index],epsilon_temp)[0] for index in range(len(t_CML))])
    R_hat_PLM_CML = np.array([symmetry_toolbox.PLM_2_symmetry(t_CML[index],R_CML[index],epsilon_temp)[1] for index in range(len(t_CML))])
    # IM-II
    t_hat_IM_II_CML = np.array([symmetry_toolbox.IM_II_symmetry(t_CML[index],R_CML[index],epsilon_temp,30,0.044)[0] for index in range(len(t_CML))])
    R_hat_IM_II_CML = np.array([symmetry_toolbox.IM_II_symmetry(t_CML[index],R_CML[index],epsilon_temp,30,0.044)[1] for index in range(len(t_CML))])    
    # FIT MODELS TO TRANSFORMED DATA
    RMS_PLM_myeloma_ODR.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM_myeloma,R_hat_PLM_myeloma,"PLM","ODR",[])[2])
    RMS_IM_II_myeloma_ODR.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II_myeloma,R_hat_IM_II_myeloma,"IM-II","ODR",[])[2])
    RMS_PLM_colon_ODR.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM_colon,R_hat_PLM_colon,"PLM","ODR",[])[2])
    RMS_IM_II_colon_ODR.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II_colon,R_hat_IM_II_colon,"IM-II","ODR",[])[2])    
    RMS_PLM_CML_ODR.append(fit_to_data_ODR.PE_risk_profiles(t_hat_PLM_CML,R_hat_PLM_CML,"PLM","ODR",[])[2])
    RMS_IM_II_CML_ODR.append(fit_to_data_ODR.PE_risk_profiles(t_hat_IM_II_CML,R_hat_IM_II_CML,"IM-II","ODR",[])[2])    
# =================================================================================
# =================================================================================
# Symmetry based model selection
# =================================================================================
# =================================================================================
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(3,1,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: Myeloma
axes[0].plot(epsilon_vec,np.array(RMS_PLM_myeloma_ODR),'-',color = (103/256,0/256,31/256),label='PLM Myeloma')
axes[0].plot(epsilon_vec,np.array(RMS_IM_II_myeloma_ODR),'-',color = (2/256,56/256,88/256),label='IM_II Myeloma')
axes[0].legend()
# Subplot 2: Colon
axes[1].plot(epsilon_vec,np.array(RMS_PLM_colon_ODR),'-',color = (103/256,0/256,31/256),label='PLM Colon')
axes[1].plot(epsilon_vec,np.array(RMS_IM_II_colon_ODR),'-',color = (2/256,56/256,88/256),label='IM_II Colon')
axes[1].legend()
# Subplot 3: CML
axes[2].plot(epsilon_vec,np.array(RMS_PLM_CML_ODR),'-',color = (103/256,0/256,31/256),label='PLM CML')
axes[2].plot(epsilon_vec,np.array(RMS_IM_II_CML_ODR),'-',color = (2/256,56/256,88/256),label='IM_II CML')
axes[2].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("Symmetry based model selection using ODR",fontsize=20, fontweight='bold')
plt.savefig("../Figures/symmetry_based_model_selection_ODR.png")
plt.show()
