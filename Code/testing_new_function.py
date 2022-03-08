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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define limits for t as well
t_min = 30
t_max = 150
# Define our parameters
#K = 80
#gamma = 1.5
#A = 1600
# Make a time vector
time_vec = np.linspace(t_min,t_max,100,endpoint=True)
# Define the function we want to fit
#R_hat_PLM_II = np.array([fit_to_data.objective_PLM_II((A,gamma,K),t_i) for t_i in list(time_vec)])
#t = np.linspace(t_min,t_max,100,endpoint=True)
#epsilon = np.linspace(t_min,t_max,100,endpoint=True)
epsilon = np.linspace(t_min,t_max,100,endpoint=True)
#R_hat_PLM_II = np.arccosh(0.5*t**(5))
a = 0.12
C = 1
tau = 30
alpha = 0.12
#R_hat_PLM_II = a*((1+C*(t**(2*a*gamma)))/(1-C*(t**(2*a*gamma))))
#R_hat_PLM_II = ((C*t**(a*gamma))/(1+(C*t**(a*gamma))))
#R_hat_PLM_II = ((1)/(1+np.exp(-t)))

epsilon_crit = np.array([C-np.exp(np.exp(-alpha*(t-tau))) for t in list(time_vec)])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the data and our new function PLM-II in Python
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT OF THE FIT OF THE MODEL TO THE DATA
# Overall properties
#fig, axes = plt.subplots(1,1,figsize=(15,5))
fig, axes = plt.subplots(1,1,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
#axes.plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
#axes.plot(t, R_hat_PLM_II, '-', color = (127/256,39/256,4/256),label='ODR fit PLM-II')
#axes.plot(epsilon, t_hat, '-', color = (127/256,39/256,4/256),label='$\hat{t}$')
axes.plot(time_vec,epsilon_crit, '-', color = (127/256,39/256,4/256),label='$\epsilon(t)$')
axes.legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.xlim(t_min,t_max)
plt.title("Plotting PLM-II along side the data",fontsize=20, fontweight='bold')
plt.show()
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
# PLM-II 
#PLM_II_fitted_to_myeloma_ODR, R_hat_PLM_II_myeloma_ODR, RMS_PLM_II_myeloma_ODR  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM-II","ODR",[],[])
#write_output.save_data_PE("myeloma", "PLM-II", PLM_II_fitted_to_myeloma_ODR, RMS_PLM_II_myeloma_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# COLON CANCER DATA
# PLM-II
#PLM_II_fitted_to_colon_ODR, R_hat_PLM_II_colon_ODR, RMS_PLM_II_colon_ODR  = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM-II","ODR",[],[])
#write_output.save_data_PE("colon", "PLM-II", PLM_II_fitted_to_colon_ODR, RMS_PLM_II_colon_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# PLM-II
#PLM_II_fitted_to_CML_ODR, R_hat_PLM_II_CML_ODR, RMS_PLM_II_CML_ODR  = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM-II","ODR",[],[])
#write_output.save_data_PE("CML", "PLM-II", PLM_II_fitted_to_CML_ODR, RMS_PLM_II_CML_ODR, "ODR")

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the data and the fit in Python
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT OF THE FIT OF THE MODEL TO THE DATA
# Overall properties
#fig, axes = plt.subplots(1,3,figsize=(15,5))
#fig, axes = plt.subplots(1,1,figsize=(15,5))
#plt.rc('axes', labelsize=15)    # fontsize of the x and y label
#plt.rc('legend', fontsize=10)    # legend fontsize
#plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
#axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
#axes[0].plot(t_myeloma, R_hat_PLM_II_myeloma_ODR, '-', color = (127/256,39/256,4/256),label='ODR fit PLM-II')
#axes[0].legend()
# Subplot 2a
#axes[1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
#axes[1].plot(t_colon, R_hat_PLM_II_colon_ODR, '-', color = (127/256,39/256,4/256),label='ODR fit PLM-II')
#axes[1].legend()
# Subplot 3a
#axes[2].plot(t_CML, R_CML, '*', color='black', label='Data CML')
#axes[2].plot(t_CML, R_hat_PLM_II_CML_ODR, '-', color = (127/256,39/256,4/256),label='ODR fit PLM-II')
#axes[2].legend()
# add a big axis, hide frame
#fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
#plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
#plt.xlabel("Age, $t$")
#plt.ylabel("Incidence, $R(t)$")
# displaying the title
#plt.title("Fit of candidate models to three data sets",fontsize=20, fontweight='bold')
#plt.show()
#plt.savefig("../Figures/Fit_new_model.png")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the data and the fit in latex
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT DATA
# Myeloma
#write_output.plot_LaTeX_2D(t_myeloma,R_myeloma,"../Figures/latex_figures/fit_new_model/Input/myeloma.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Myeloma cancer")
# Colon cancer
#write_output.plot_LaTeX_2D(t_colon,R_colon,"../Figures/latex_figures/fit_new_model/Input/colon.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Colon cancer")
# CML
#write_output.plot_LaTeX_2D(t_CML,R_CML,"../Figures/latex_figures/fit_new_model/Input/CML.tex","only marks, mark=halfcircle*,mark size=1.5pt,color=black,","Data Chronic Myeloid Leukemia (CML)")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Orthogonal Distance Regression (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT PLM_II
# Myeloma
#write_output.plot_LaTeX_2D(t_myeloma,R_hat_PLM_II_myeloma_ODR,"../Figures/latex_figures/fit_new_model/Input/myeloma.tex","color=pow_2_1,line width=2pt,","PLM-II")
# Colon cancer
#write_output.plot_LaTeX_2D(t_colon,R_hat_PLM_II_colon_ODR,"../Figures/latex_figures/fit_new_model/Input/colon.tex","color=pow_2_1,line width=2pt,","PLM-II")
# CML
#write_output.plot_LaTeX_2D(t_CML,R_hat_PLM_II_CML_ODR,"../Figures/latex_figures/fit_new_model/Input/CML.tex","color=pow_2_1,line width=2pt,","PLM-II")
#----------------------------------------------------------------------------------
