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
import symmetry_based_model_selection # Home-made
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
t_myeloma = t_myeloma[24:len(t_myeloma)-1]
R_myeloma = R_myeloma[24:len(R_myeloma)-1]
# COLON
xlabel_colon, ylabel_colon, t_colon, R_colon = read_data.time_series_from_csv("Colon_cancer")
# Only ages above 12 as we have zero incidences below
t_colon = t_colon[11:len(t_colon)-1]
R_colon = R_colon[11:len(R_colon)-1]
# Chronic Myeloid Leukemia (CML)
xlabel_CML, ylabel_CML, t_CML, R_CML = read_data.time_series_from_csv("CML_cancer")
# Only ages above 10 as we have zero incidences below
t_CML = t_CML[9:len(t_CML)-1]
R_CML = R_CML[9:len(R_CML)-1]
# =================================================================================
# =================================================================================
# Fit the models to the data
# =================================================================================
# =================================================================================
# Define parameter alpha common to most models
alpha = 0.044
# ---------------------------------------------------------------------------------
# Fit IM-II
# ---------------------------------------------------------------------------------
# Myeloma data
R_hat_IM_II_Myeloma, opt_para_IM_II_Myeloma, SS_IM_II_Myeloma = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II",[alpha])
# Save the estimated parameters: Myeloma
write_output.save_data_PE("Myeloma","IM-II",opt_para_IM_II_Myeloma,SS_IM_II_Myeloma,alpha)
# Colon data
R_hat_IM_II_colon, opt_para_IM_II_colon, SS_IM_II_colon = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-II",[alpha])
# Save the estimated parameters: Myeloma
write_output.save_data_PE("Colon_cancer","IM-II",opt_para_IM_II_colon,SS_IM_II_colon,alpha)
# CML data
R_hat_IM_II_CML, opt_para_IM_II_CML, SS_IM_II_CML = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-II",[alpha])
# Save the estimated parameters: Myeloma
write_output.save_data_PE("CML","IM-II",opt_para_IM_II_CML,SS_IM_II_CML,alpha)
# ---------------------------------------------------------------------------------
# Fit PLM
# ---------------------------------------------------------------------------------
# Myeloma data
R_hat_PLM_Myeloma, opt_para_PLM_Myeloma, SS_PLM_Myeloma = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM",[alpha])
# Save the estimated parameters: Myeloma
write_output.save_data_PE("Myeloma","PLM",opt_para_IM_II_Myeloma,SS_IM_II_Myeloma,alpha)
# Colon data
R_hat_PLM_colon, opt_para_PLM_colon, SS_PLM_colon = fit_to_data.PE_risk_profiles(t_colon,R_colon,"PLM",[alpha])
# Save the estimated parameters: Myeloma
write_output.save_data_PE("Colon_cancer","PLM",opt_para_IM_II_colon,SS_IM_II_colon,alpha)
# CML data
R_hat_PLM_CML, opt_para_PLM_CML, SS_PLM_CML = fit_to_data.PE_risk_profiles(t_CML,R_CML,"PLM",[alpha])
# Save the estimated parameters: Myeloma
write_output.save_data_PE("CML","PLM",opt_para_PLM_CML,SS_PLM_CML,alpha)
# =================================================================================
# =================================================================================
# Plot the data and the fit
# =================================================================================
# =================================================================================
# Pyplot
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1
axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0].plot(t_myeloma, R_hat_IM_II_Myeloma, '-', color='magenta',label='IM-II')
axes[0].plot(t_myeloma, R_hat_PLM_Myeloma, '-', color='blue',label='PLM')
axes[0].legend()
# Subplot 2
axes[1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1].plot(t_colon, R_hat_IM_II_colon, '-', color='magenta',label='IM-II')
axes[1].plot(t_colon, R_hat_PLM_colon, '-', color='blue',label='PLM')
axes[1].legend()
# Subplot 3
axes[2].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2].plot(t_CML, R_hat_IM_II_CML, '-', color='magenta',label='IM-II')
axes[2].plot(t_CML, R_hat_PLM_CML, '-', color='blue',label='PLM')
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








