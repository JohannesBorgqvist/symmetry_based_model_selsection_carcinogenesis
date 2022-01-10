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
PLM_fitted_to_myeloma_LS, R_hat_PLM_myeloma_LS, RMS_PLM_myeloma_LS  = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","LS")
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_LS, RMS_PLM_myeloma_LS, "LS")
# IM-II
IM_II_fitted_to_myeloma_LS, R_hat_IM_II_myeloma_LS, RMS_IM_II_myeloma_LS = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","LS")
write_output.save_data_PE("myeloma", "IM-II", IM_II_fitted_to_myeloma_LS, RMS_IM_II_myeloma_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_LS, R_hat_PLM_colon_LS, RMS_PLM_colon_LS  = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"PLM","LS")
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_LS, RMS_PLM_colon_LS, "LS")
# IM-II
IM_II_fitted_to_colon_LS, R_hat_IM_II_colon_LS, RMS_IM_II_colon_LS = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"IM-II","LS")
write_output.save_data_PE("colon", "IM-II", IM_II_fitted_to_colon_LS, RMS_IM_II_colon_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_LS, R_hat_PLM_CML_LS, RMS_PLM_CML_LS  = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"PLM","LS")
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_LS, RMS_PLM_CML_LS, "LS")
# IM-II
IM_II_fitted_to_CML_LS, R_hat_IM_II_CML_LS, RMS_IM_II_CML_LS = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"IM-II","LS")
write_output.save_data_PE("CML", "IM-II", IM_II_fitted_to_CML_LS, RMS_IM_II_CML_LS, "LS")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ORTHOGONAL DISTANCE REGRESSION (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_ODR, R_hat_PLM_myeloma_ODR, RMS_PLM_myeloma_ODR  = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","ODR")
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_ODR, RMS_PLM_myeloma_ODR, "ODR")
# IM-II
IM_II_fitted_to_myeloma_ODR, R_hat_IM_II_myeloma_ODR, RMS_IM_II_myeloma_ODR = fit_to_data_ODR.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","ODR")
write_output.save_data_PE("myeloma", "IM-II", IM_II_fitted_to_myeloma_ODR, RMS_IM_II_myeloma_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_ODR, R_hat_PLM_colon_ODR, RMS_PLM_colon_ODR  = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"PLM","ODR")
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_ODR, RMS_PLM_colon_ODR, "ODR")
# IM-II
IM_II_fitted_to_colon_ODR, R_hat_IM_II_colon_ODR, RMS_IM_II_colon_ODR = fit_to_data_ODR.PE_risk_profiles(t_colon,R_colon,"IM-II","ODR")
write_output.save_data_PE("colon", "IM-II", IM_II_fitted_to_colon_ODR, RMS_IM_II_colon_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM 
PLM_fitted_to_CML_ODR, R_hat_PLM_CML_ODR, RMS_PLM_CML_ODR  = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"PLM","ODR")
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_ODR, RMS_PLM_CML_ODR, "ODR")
# IM-II
IM_II_fitted_to_CML_ODR, R_hat_IM_II_CML_ODR, RMS_IM_II_CML_ODR = fit_to_data_ODR.PE_risk_profiles(t_CML,R_CML,"IM-II","ODR")
write_output.save_data_PE("CML", "IM-II", IM_II_fitted_to_CML_ODR, RMS_IM_II_CML_ODR, "ODR")
#----------------------------------------------------------------------------------
# =================================================================================
# =================================================================================
# Plot the data and the fit in Python
# =================================================================================
# =================================================================================
# Pyplot
# Overall properties
fig, axes = plt.subplots(3,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0][0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0][0].plot(t_myeloma, R_hat_IM_II_myeloma_LS, '-', color='magenta',label='LS fit IM-II')
axes[0][0].plot(t_myeloma, R_hat_PLM_myeloma_LS, '-', color='blue',label='LS fit PLM')
axes[0][0].legend()
# Subplot 1b
axes[0][1].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0][1].plot(t_myeloma, R_hat_IM_II_myeloma_ODR, '-', color='magenta',label='ODR fit IM-II')
axes[0][1].plot(t_myeloma, R_hat_PLM_myeloma_ODR, '-', color='blue',label='ODR fit PLM')
axes[0][1].legend()
# Subplot 2a
axes[1][0].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1][0].plot(t_colon, R_hat_IM_II_colon_LS, '-', color='magenta',label='LS fit IM-II')
axes[1][0].plot(t_colon, R_hat_PLM_colon_LS, '-', color='blue',label='LS fit PLM')
axes[1][0].legend()
# Subplot 2b
axes[1][1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1][1].plot(t_colon, R_hat_IM_II_colon_ODR, '-', color='magenta',label='ODR fit IM-II')
axes[1][1].plot(t_colon, R_hat_PLM_colon_ODR, '-', color='blue',label='ODR fit PLM')
axes[1][1].legend()
# Subplot 3a
axes[2][0].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2][0].plot(t_CML, R_hat_IM_II_CML_LS, '-', color='magenta',label='LS fit IM-II')
axes[2][0].plot(t_CML, R_hat_PLM_CML_LS, '-', color='blue',label='LS fit PLM')
axes[2][0].legend()
# Subplot 3b
axes[2][1].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2][1].plot(t_CML, R_hat_IM_II_CML_ODR, '-', color='magenta',label='ODR fit IM-II')
axes[2][1].plot(t_CML, R_hat_PLM_CML_ODR, '-', color='blue',label='ODR fit PLM')
axes[2][1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("Fit of candidate models to three data sets",fontsize=20, fontweight='bold')
plt.show()
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
