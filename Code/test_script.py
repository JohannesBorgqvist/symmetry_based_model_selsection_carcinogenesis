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





# MYELOMA
xlabel_myeloma, ylabel_myeloma, t_myeloma, R_myeloma = read_data.time_series_from_csv("Myeloma_cancer")

# Only ages above 25 as we have zero incidences below
t_myeloma = np.array(list(t_myeloma[24:len(t_myeloma)-1]))
R_myeloma = np.array(list(R_myeloma[24:len(R_myeloma)-1]))

# IM-II
IM_II_fitted_to_myeloma_ODR, R_hat_IM_II_myeloma_ODR, RMS_IM_II_myeloma_ODR = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"IM-II","ODR")



# PLM 
PLM_fitted_to_myeloma_ODR, R_hat_PLM_myeloma_ODR, RMS_PLM_myeloma_ODR  = fit_to_data.PE_risk_profiles(t_myeloma,R_myeloma,"PLM","ODR")


# Construct a t vector
t_sym = np.linspace(t_myeloma[0],t_myeloma[len(t_myeloma)-1],50)

# THE ORIGINAL VECTORS OR ORIGINAL DATA IF YOU WILL
R_hat_PLM_original = np.array([fit_to_data.objective_PLM(PLM_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
R_hat_IM_II_original = np.array([fit_to_data.objective_IM_II(IM_II_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])

#PLM_fit_original, R_hat_PLM_fit, RMS_PLM_original  = fit_to_data.PE_risk_profiles(t_sym,R_hat_PLM_original,"PLM","ODR")
# IM-II
#IM_II_fit_original, R_hat_IM_II_fit, RMS_IM_II_original = fit_to_data.PE_risk_profiles(t_sym,R_hat_IM_II_original,"IM-II","ODR")

#print("RMS PLM ORIGINAL")
#print(RMS_PLM_original)
#print(PLM_fit_original.sum_square_delta)
#print(RMS_PLM_original)
#print(R_hat_PLM_original-R_hat_PLM_fit)
#print("RMS IM-II ORIGINAL")
#print(RMS_IM_II_original)
#print(IM_II_fit_original.sum_square_delta)
#print(RMS_IM_II_original)
#print(R_hat_IM_II_original-R_hat_IM_II_fit)


epsilon_temp = 0.5

t_hat_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[0] for index in range(len(t_sym))])
R_hat_PLM = np.array([symmetry_toolbox.PLM_2_symmetry(t_sym[index],R_hat_PLM_original[index],epsilon_temp)[1] for index in range(len(t_sym))])
# IM-II
t_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_temp,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[0] for index in range(len(t_sym))])
R_hat_IM_II = np.array([symmetry_toolbox.IM_II_symmetry(t_sym[index],R_hat_IM_II_original[index],epsilon_temp,IM_II_fitted_to_myeloma_ODR.beta[1],0.044)[1] for index in range(len(t_sym))])
# Fit to transformed
PLM_fit_trans, R_hat_PLM_fit_trans, RMS_PLM_trans  = fit_to_data.PE_risk_profiles(t_hat_PLM,R_hat_PLM,"PLM","LS")
# IM-II
IM_II_fit_trans, R_hat_IM_II_fit_trans, RMS_IM_II_trans = fit_to_data.PE_risk_profiles(t_hat_IM_II,R_hat_IM_II,"IM-II","LS")


print("RMS PLM TRANSFORMED")
print(RMS_PLM_trans)
print(PLM_fit_trans.sum_square_delta)
print(RMS_PLM_trans)
print(R_hat_PLM-R_hat_PLM_fit_trans)
print("RMS IM-II TRANSFORMED")
print(RMS_IM_II_trans)
print(IM_II_fit_trans.sum_square_delta)
print(RMS_IM_II_trans)
print(R_hat_IM_II-R_hat_IM_II_fit_trans)



# PLOT
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1
#axes[0].plot(t_sym,R_hat_PLM_original,'*',color = 'black',label='Original PLM')
#axes[0].plot(t_sym,R_hat_PLM_fit,'-',color = (103/256,0/256,31/256),label='Fitted PLM Original')
axes[0].plot(t_hat_PLM,R_hat_PLM,'*',color = 'black',label='transformed PLM')
axes[0].plot(t_hat_PLM,R_hat_PLM_fit_trans,'-',color=(206/256,18/256,86/256),label='Fitted PLM Transformed')
axes[0].legend()
# Subplot 2
#axes[1].plot(t_sym,R_hat_IM_II_original,'*',color = 'black',label='Original IM_II')
#axes[1].plot(t_sym,R_hat_IM_II_fit,'-',color = (2/256,56/256,88/256),label='Fitted IM_II Original')
axes[1].plot(t_hat_IM_II,R_hat_IM_II,'*',color = 'black',label='transformed IM_II')
#axes[1].plot(t_hat_IM_II,R_hat_IM_II_fit_trans,'-',color=(54/256,144/256,192/256),label='Fitted IM_II transformed')
axes[1].legend()
plt.show()
