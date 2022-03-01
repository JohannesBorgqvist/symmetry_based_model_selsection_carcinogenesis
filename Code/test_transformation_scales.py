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
# IM-III
IM_III_fitted_to_colon_ODR, R_hat_IM_III_colon_ODR, RMS_IM_III_colon_ODR = fit_to_data.PE_risk_profiles(t_colon,R_colon,"IM-III","ODR",[],[])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML DATA
# IM-III
#IM_III_fitted_to_CML_ODR, R_hat_IM_III_CML_ODR, RMS_IM_III_CML_ODR = fit_to_data.PE_risk_profiles(t_CML,R_CML,"IM-III","ODR",[],[])
# =================================================================================
# =================================================================================
# CALCULATE THE TRANSFORMATION SCALES
# =================================================================================
# Colon data
epsilon_scale_IM_III_colon = symmetry_toolbox.IM_III_transformation_scale(80,2,IM_III_fitted_to_colon_ODR.beta[3],IM_III_fitted_to_colon_ODR.beta[1])
# CML data
#epsilon_scale_IM_III_CML = symmetry_toolbox.IM_III_transformation_scale(80,2,IM_III_fitted_to_CML_ODR.beta[3],IM_III_fitted_to_CML_ODR.beta[1])
# =================================================================================
# =================================================================================
# CONDUCT THE SYMMETRY BASED MODEL SELECTION
# =================================================================================
# =================================================================================
epsilon_vector_IM_III_colon = np.linspace(0.0,0.89*epsilon_scale_IM_III_colon,num=10,endpoint=True)
#epsilon_vector_IM_III_CML = np.linspace(0.0,0.8286*epsilon_scale_IM_III_CML,num=10,endpoint=True)
print("\t\tModel\t=\t IM-III,\tDataset\t=\t Colon cancer")
RMS_transf_IM_III_colon = symmetry_toolbox.symmetry_based_model_selection(t_colon,R_colon,epsilon_vector_IM_III_colon,IM_III_fitted_to_colon_ODR.beta,"IM-III")
print("\t\t\tDone!\n")
#print("\t\tModel\t=\t IM-III,\tDataset\t=\t CML")
#RMS_transf_IM_III_CML = symmetry_toolbox.symmetry_based_model_selection(t_CML,R_CML,epsilon_vector_IM_III_CML,IM_III_fitted_to_CML_ODR.beta,"IM-III")
#print("\t\t\tDone!\n")
# IM-III
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: Colon cancer
axes[0].plot(epsilon_vector_IM_III_colon,RMS_transf_IM_III_colon,'-', color = (54/256,144/256,192/256),label='IM-III Colon cancer')
# Subplot 2: CML
#axes[1].plot(epsilon_vector_IM_III_CML,RMS_transf_IM_III_CML,'-', color = (208/256,209/256,230/256),label='IM-III CML')
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("The symmetry based model selection for the IM-III",fontsize=20, fontweight='bold')
plt.show()    

