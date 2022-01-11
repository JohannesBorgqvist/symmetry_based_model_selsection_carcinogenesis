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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# STANDARD LEAST SQUARE (LS) FITTING
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_LS, R_hat_PLM_myeloma_LS, RMS_PLM_myeloma_LS  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","LS")
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_LS, RMS_PLM_myeloma_LS, "LS")
# IM-II
IM_II_fitted_to_myeloma_LS, R_hat_IM_II_myeloma_LS, RMS_IM_II_myeloma_LS = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","LS")
write_output.save_data_PE("myeloma", "IM-II", IM_II_fitted_to_myeloma_LS, RMS_IM_II_myeloma_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_LS, R_hat_PLM_colon_LS, RMS_PLM_colon_LS  = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM","LS")
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_LS, RMS_PLM_colon_LS, "LS")
# IM-II
IM_II_fitted_to_colon_LS, R_hat_IM_II_colon_LS, RMS_IM_II_colon_LS = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-II","LS")
write_output.save_data_PE("colon", "IM-II", IM_II_fitted_to_colon_LS, RMS_IM_II_colon_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_LS, R_hat_PLM_CML_LS, RMS_PLM_CML_LS  = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM","LS")
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_LS, RMS_PLM_CML_LS, "LS")
# IM-II
IM_II_fitted_to_CML_LS, R_hat_IM_II_CML_LS, RMS_IM_II_CML_LS = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-II","LS")
write_output.save_data_PE("CML", "IM-II", IM_II_fitted_to_CML_LS, RMS_IM_II_CML_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ORTHOGONAL DISTANCE REGRESSION (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_ODR, R_hat_PLM_myeloma_ODR, RMS_PLM_myeloma_ODR  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","ODR")
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_ODR, RMS_PLM_myeloma_ODR, "ODR")
# IM-II
IM_II_fitted_to_myeloma_ODR, R_hat_IM_II_myeloma_ODR, RMS_IM_II_myeloma_ODR = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","ODR")
write_output.save_data_PE("myeloma", "IM-II", IM_II_fitted_to_myeloma_ODR, RMS_IM_II_myeloma_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_ODR, R_hat_PLM_colon_ODR, RMS_PLM_colon_ODR  = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM","ODR")
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_ODR, RMS_PLM_colon_ODR, "ODR")
# IM-II
IM_II_fitted_to_colon_ODR, R_hat_IM_II_colon_ODR, RMS_IM_II_colon_ODR = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-II","ODR")
write_output.save_data_PE("colon", "IM-II", IM_II_fitted_to_colon_ODR, RMS_IM_II_colon_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_ODR, R_hat_PLM_CML_ODR, RMS_PLM_CML_ODR  = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM","ODR")
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_ODR, RMS_PLM_CML_ODR, "ODR")
# IM-II
IM_II_fitted_to_CML_ODR, R_hat_IM_II_CML_ODR, RMS_IM_II_CML_ODR = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-II","ODR")
write_output.save_data_PE("CML", "IM-II", IM_II_fitted_to_CML_ODR, RMS_IM_II_CML_ODR, "ODR")
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
#plt.savefig("../Figures/Fit_of_models_to_cancer_data.png")
# =================================================================================
# =================================================================================
# Plot the data and the fit in latex
# =================================================================================
# =================================================================================
# PLOT DATA
# Define the string defining the settings for the plot
plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/FigS4/myeloma/Input/data.tex",plot_str,"Data Myeloma cancer")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_colon,"../Figures/latex_figures/FigS4/colon/Input/data.tex",plot_str,"Data Colon cancer")
# CML
write_output.plot_LaTeX_2D(t_CML,R_CML,"../Figures/latex_figures/FigS4/CML/Input/data.tex",plot_str,"Data CML")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# STANDARD LEAST SQUARE (LS) FITTING
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT PLM
# Define the string defining the settings for the plot
plot_str = "color=pow_1,line width=2pt,"
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_PLM_myeloma_LS,"../Figures/latex_figures/FigS4/myeloma/Input/PLM_LS.tex",plot_str,"PLM LS fit")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_PLM_colon_LS,"../Figures/latex_figures/FigS4/colon/Input/PLM_LS.tex",plot_str,"PLM LS fit")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_PLM_CML_LS,"../Figures/latex_figures/FigS4/CML/Input/PLM_LS.tex",plot_str,"PLM LS fit")
# PLOT IM-II
# Define the string defining the settings for the plot
plot_str = "color=mixed_1,line width=2pt,"
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_IM_II_myeloma_LS,"../Figures/latex_figures/FigS4/myeloma/Input/IM-II_LS.tex",plot_str,"IM-II LS fit")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_IM_II_colon_LS,"../Figures/latex_figures/FigS4/colon/Input/IM-II_LS.tex",plot_str,"IM-II LS fit")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_IM_II_CML_LS,"../Figures/latex_figures/FigS4/CML/Input/IM-II_LS.tex",plot_str,"IM-II LS fit")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ORTHOGONAL DISTANCE REGRESSION (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT PLM
# Define the string defining the settings for the plot
plot_str = "color=pow_1,line width=2pt,"
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_PLM_myeloma_ODR,"../Figures/latex_figures/FigS4/myeloma/Input/PLM_ODR.tex",plot_str,"PLM ODR fit")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_PLM_colon_ODR,"../Figures/latex_figures/FigS4/colon/Input/PLM_ODR.tex",plot_str,"PLM ODR fit")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_PLM_CML_ODR,"../Figures/latex_figures/FigS4/CML/Input/PLM_ODR.tex",plot_str,"PLM ODR fit")
# PLOT IM-II
# Define the string defining the settings for the plot
plot_str = "color=mixed_1,line width=2pt,"
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_hat_IM_II_myeloma_ODR,"../Figures/latex_figures/FigS4/myeloma/Input/IM-II_ODR.tex",plot_str,"IM-II ODR fit")
# Colon cancer
write_output.plot_LaTeX_2D(t_colon,R_hat_IM_II_colon_ODR,"../Figures/latex_figures/FigS4/colon/Input/IM-II_ODR.tex",plot_str,"IM-II ODR fit")
# CML
write_output.plot_LaTeX_2D(t_CML,R_hat_IM_II_CML_ODR,"../Figures/latex_figures/FigS4/CML/Input/IM-II_ODR.tex",plot_str,"IM-II ODR fit")
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
R_hat_PLM_original = np.array([fit_to_data.objective_PLM(PLM_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Choose an epsilon
epsilon = 0.35
# Allocate memory for a list
R_PLM_trans_1 = []
t_PLM_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(10,len(t_sym)-1,75,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon-0.035,PLM_fitted_to_myeloma_ODR.beta[1],2)
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
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon-0.035,PLM_fitted_to_myeloma_ODR.beta[1],2)
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2 = np.array([symmetry_toolbox.PLM_2_symmetry(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon)[0] for index in range(len(t_sym))])
R_hat_PLM_2 = np.array([symmetry_toolbox.PLM_2_symmetry(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon)[1] for index in range(len(t_sym))])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# THE IM-II
R_hat_IM_II_original = np.array([fit_to_data.objective_IM_II(IM_II_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_II_trans_1 = []
t_IM_II_trans_1 = []
# Allocate an index vector
index_vector = list(np.linspace(175,len(t_sym)-1,25,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_sym[index],R_hat_IM_II_original[index],epsilon-0.035,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)
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
    t_trans,R_trans = symmetry_toolbox.IM_II_transformation(t_hat_IM_II_1[index],R_hat_IM_II_1[index],epsilon-0.035,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)
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
#-------------------------------------------------------------------------------------
# The LaTeX plot as well
#-------------------------------------------------------------------------------------
# PLM
write_output.plot_LaTeX_2D(t_sym, R_hat_PLM_original,"../Figures/latex_figures/symmetries/Input/PLM.tex","color=pow_1,line width=2pt,","PLM original solution, $R(t)$")
# The first symmetry
write_output.plot_LaTeX_2D(np.array(t_PLM_trans_1[0]),np.array(R_PLM_trans_1[0]),"../Figures/latex_figures/symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed","Action of symmetry, $\\Gamma_{1,1}(\\epsilon)$")
# The remaining ones
for index in range(1,len(t_PLM_trans_1)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans_1[index]),np.array(R_PLM_trans_1[index]),"../Figures/latex_figures/symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])
# The remaining ones
for index in range(1,len(t_PLM_trans_2)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans_2[index]),np.array(R_PLM_trans_2[index]),"../Figures/latex_figures/symmetries/Input/PLM.tex","color=black,->,>=latex,densely dashed",[])
# Transformed solution 1
write_output.plot_LaTeX_2D(t_hat_PLM_1, R_hat_PLM_1,"../Figures/latex_figures/symmetries/Input/PLM.tex","color=pow_2,line width=2pt,","PLM transformed solution, $\\hat{R}_1\\left(\\hat{t}\\right)$")
# Transformed solution 2
write_output.plot_LaTeX_2D(t_hat_PLM_2, R_hat_PLM_2,"../Figures/latex_figures/symmetries/Input/PLM.tex","color=pow_3,line width=2pt,","PLM transformed solution, $\\hat{R}_2\\left(\\hat{t}\\right)$")
# IM_II
write_output.plot_LaTeX_2D(t_sym, R_hat_IM_II_original,"../Figures/latex_figures/symmetries/Input/IM_II.tex","color=mixed_1,line width=2pt,","IM-II original solution, $R(t)$")
# The first symmetry
write_output.plot_LaTeX_2D(np.array(t_IM_II_trans_1[0]),np.array(R_IM_II_trans_1[0]),"../Figures/latex_figures/symmetries/Input/IM_II.tex","color=black,->,>=latex,densely dashed","Action of symmetry, $\\Gamma_{2,1}(\\epsilon)$")
# The remaining ones
for index in range(1,len(t_IM_II_trans_1)):
    write_output.plot_LaTeX_2D(np.array(t_IM_II_trans_1[index]),np.array(R_IM_II_trans_1[index]),"../Figures/latex_figures/symmetries/Input/IM_II.tex","color=black,->,>=latex,densely dashed",[])
# The remaining ones
for index in range(1,len(t_IM_II_trans_2)):
    write_output.plot_LaTeX_2D(np.array(t_IM_II_trans_2[index]),np.array(R_IM_II_trans_2[index]),"../Figures/latex_figures/symmetries/Input/IM_II.tex","color=black,->,>=latex,densely dashed",[])
# Transformed solution 1
write_output.plot_LaTeX_2D(t_hat_IM_II_1, R_hat_IM_II_1,"../Figures/latex_figures/symmetries/Input/IM_II.tex","color=mixed_2,line width=2pt,","IM-II transformed solution, $\\hat{R}_1\\left(\\hat{t}\\right)$")
# Transformed solution 2
write_output.plot_LaTeX_2D(t_hat_IM_II_2, R_hat_IM_II_2,"../Figures/latex_figures/symmetries/Input/IM_II.tex","color=mixed_3,line width=2pt,","IM-II transformed solution, $\\hat{R}_2\\left(\\hat{t}\\right)$")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustrate framework
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Pick a value of epsilon
epsilon = 0.3
# Allocate memory for a list
R_PLM_trans = []
t_PLM_trans = []
# Allocate an index vector
index_vector = list(np.linspace(10,55,30,dtype=int))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_myeloma[index],R_myeloma[index],epsilon-0.01,PLM_fitted_to_myeloma_ODR.beta[1],2)
    # Save the transformed variables
    R_PLM_trans.append(R_trans)
    t_PLM_trans.append(t_trans)
# Transform the data
t_myeloma_trans = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon)[0] for index in range(len(t_myeloma))])
R_myeloma_trans = np.array([symmetry_toolbox.PLM_2_symmetry(t_myeloma[index],R_myeloma[index],epsilon)[1] for index in range(len(t_myeloma))])
# Fit the model to the transformed data
param_temp, R_hat_PLM_trans, RMS_temp  = fit_to_data.PE_risk_profiles(t_myeloma_trans,R_myeloma_trans,"PLM","ODR")
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
# Subplot 1a
axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0].plot(np.array(t_PLM_trans[0]),np.array(R_PLM_trans[0]),'--',color='black',label='Symmetry')
for index in range(1,len(t_PLM_trans)):
    axes[0].plot(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),'--',color='black')
axes[0].plot(t_myeloma_trans, R_myeloma_trans, '*', color='gray', label='Transformed Data')
axes[0].legend()
# Subplot 1a
axes[1].plot(t_myeloma_trans, R_myeloma_trans, '*', color='gray', label='Transformed Data')
axes[1].plot(t_myeloma_trans, R_hat_PLM_trans, '-', color=(206/256,18/256,86/256), label='Fitted model')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The steps of the symmetry based model selection",fontsize=20, fontweight='bold')
plt.show()
plt.savefig("../Figures/illustrate_framework.png")
#-------------------------------------------------------------------------------------
# The LaTeX plot as well
#-------------------------------------------------------------------------------------
# Define the string defining the settings for the plot
plot_str = "only marks, mark=halfcircle*,mark size=1.5pt,color=black,"
# Myeloma
write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/illustrate_framework/Input/data.tex",plot_str,"Data Myeloma cancer")
# Plot the symmetry
plot_str = "color=black,->,>=latex,densely dashed"
# The first symmetry
write_output.plot_LaTeX_2D(np.array(t_PLM_trans[0]),np.array(R_PLM_trans[0]),"../Figures/latex_figures/illustrate_framework/Input/symmetry.tex",plot_str,"Transform data with symmetry")
# The remaining ones
for index in range(1,len(t_PLM_trans)):
    write_output.plot_LaTeX_2D(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),"../Figures/latex_figures/illustrate_framework/Input/symmetry.tex",plot_str,[])
# Add the transformed data
plot_str = "only marks, mark=diamond*,mark size=1.5pt,color=gray,every mark/.append style={solid, fill=gray}"
write_output.plot_LaTeX_2D(t_myeloma_trans,R_myeloma_trans,"../Figures/latex_figures/illustrate_framework/Input/transformed_data.tex",plot_str,"Transformed data")
# Add the fitted model to the transformed data
plot_str = "color=pow_2,line width=2pt,"
write_output.plot_LaTeX_2D(t_myeloma_trans,R_hat_PLM_trans,"../Figures/latex_figures/illustrate_framework/Input/model.tex",plot_str,"Fitted model")

