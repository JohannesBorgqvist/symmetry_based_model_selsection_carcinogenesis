# =================================================================================
# =================================================================================
# Script:"symmetry_analysis_PLM_vs_IM.py"
# Date: 2022-03-29
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
import clean_up_output # Home-made
import numpy as np # For the numerical calculations
import math # For mathematics
from matplotlib import pyplot as plt # For plotting
import multiprocessing as mp # For parallelisation
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
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tThe model fitting\n")
print("\n\t--------------------------------------------------------------------------------------\n")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Orthonal Distance Regression (ODR)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Create a list which we iterate over
input_lists = [[t_myeloma,R_myeloma,"PLM"],[t_myeloma,R_myeloma,"IM"],[t_colon,R_colon,"PLM"],[t_colon,R_colon,"IM"],[t_CML,R_CML,"PLM"],[t_CML,R_CML,"IM"]]
# We do the fitting in parallel (to use all cores, use mp.cpu_count())
pool = mp.Pool(mp.cpu_count())
# Do the fitting in parallel
results = pool.starmap(fit_to_data.PE_risk_profiles,[(input_list[0],input_list[1],input_list[2],"ODR",[],[]) for input_list in input_lists])
# Close the pool
pool.close()
# Prompt to the user
print("\t\tFitting is done!")
# Extract and save all the results
# MYELOMA DATA
# PLM 
PLM_fitted_to_myeloma_ODR = results[0][0]
R_hat_PLM_myeloma_ODR = results[0][1]
RMS_PLM_myeloma_ODR = results[0][2]
write_output.save_data_PE("myeloma", "PLM", PLM_fitted_to_myeloma_ODR, RMS_PLM_myeloma_ODR, "ODR")
# IM
IM_fitted_to_myeloma_ODR = results[1][0]
R_hat_IM_myeloma_ODR = results[1][1]
RMS_IM_myeloma_ODR = results[1][2]
write_output.save_data_PE("myeloma", "IM", IM_fitted_to_myeloma_ODR, RMS_IM_myeloma_ODR, "ODR")
# COLON CANCER DATA
# PLM 
PLM_fitted_to_colon_ODR = results[2][0]
R_hat_PLM_colon_ODR = results[2][1]
RMS_PLM_colon_ODR = results[2][2]
write_output.save_data_PE("colon", "PLM", PLM_fitted_to_colon_ODR, RMS_PLM_colon_ODR, "ODR")
# IM
IM_fitted_to_colon_ODR = results[3][0]
R_hat_IM_colon_ODR = results[3][1]
RMS_IM_colon_ODR = results[3][2]
write_output.save_data_PE("colon", "IM", IM_fitted_to_colon_ODR, RMS_IM_colon_ODR, "ODR")
# CML DATA
# PLM 
PLM_fitted_to_CML_ODR = results[4][0]
R_hat_PLM_CML_ODR = results[4][1]
RMS_PLM_CML_ODR = results[4][2]
write_output.save_data_PE("CML", "PLM", PLM_fitted_to_CML_ODR, RMS_PLM_CML_ODR, "ODR")
# IM
IM_fitted_to_CML_ODR = results[5][0]
R_hat_IM_CML_ODR = results[5][1]
RMS_IM_CML_ODR = results[5][2]
write_output.save_data_PE("CML", "IM", IM_fitted_to_CML_ODR, RMS_IM_CML_ODR, "ODR")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the data and the fit in Python
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# PLOT OF THE FIT OF THE MODEL TO THE DATA
# Overall properties
fig, axes = plt.subplots(1,3,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0].plot(t_myeloma, R_myeloma, '*', color='black', label='Data Myeloma cancer')
axes[0].plot(t_myeloma, R_hat_PLM_myeloma_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[0].plot(t_myeloma, R_hat_IM_myeloma_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM')
axes[0].legend()
# Subplot 2a
axes[1].plot(t_colon, R_colon, '*', color='black', label='Data colon cancer')
axes[1].plot(t_colon, R_hat_PLM_colon_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[1].plot(t_colon, R_hat_IM_colon_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM')
axes[1].legend()
# Subplot 3a
axes[2].plot(t_CML, R_CML, '*', color='black', label='Data CML')
axes[2].plot(t_CML, R_hat_PLM_CML_ODR, '-', color = (103/256,0/256,31/256),label='ODR fit PLM')
axes[2].plot(t_CML, R_hat_IM_CML_ODR, '-', color = (2/256,56/256,88/256),label='ODR fit IM')
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
# =================================================================================
# =================================================================================
# CALCULATE THE TRANSFORMATION SCALES
# =================================================================================
# What epsilon should we use to transform the age 85 years to 170 years? Here, comes
# the answer in the case of our two models.
# The PLM
epsilon_scale_PLM = symmetry_toolbox.PLM_transformation_scale(2)
# The IM
# Myeloma data
epsilon_scale_IM_myeloma = symmetry_toolbox.IM_transformation_scale(85,2,IM_fitted_to_myeloma_ODR.beta[3],IM_fitted_to_myeloma_ODR.beta[1])
# Colon data
epsilon_scale_IM_colon = symmetry_toolbox.IM_transformation_scale(85,2,IM_fitted_to_colon_ODR.beta[3],IM_fitted_to_colon_ODR.beta[1])
# CML data
epsilon_scale_IM_CML = symmetry_toolbox.IM_transformation_scale(85,2,IM_fitted_to_CML_ODR.beta[3],IM_fitted_to_CML_ODR.beta[1])
# Prompt to the user
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tThe transformation scales increasing the age from 85 years to 170 years\n")
print("\n\t--------------------------------------------------------------------------------------\n")
print("\t\tThe PLM:\tepsilon_PLM\t=\t%0.12f"%(2*epsilon_scale_PLM))
print("\t\tThe IM myeloma:\tepsilon_IM_myeloma\t=\t%0.12f"%(epsilon_scale_IM_myeloma))
print("\t\tThe IM colon:\tepsilon_IM_colon\t=\t%0.12f"%(epsilon_scale_IM_colon))
print("\t\tThe IM CML:\tepsilon_IM_CML\t=\t%0.12f"%(epsilon_scale_IM_CML))
# Prompt to the user
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tThe symmetry based framework for model selection\n")
print("\n\t--------------------------------------------------------------------------------------\n")
print("\t\tThe scales for the framework are the following:")
print("\t\t\tPLM\tAll datasets:\t epsilon_scale\t=\t%0.12f"%(2*epsilon_scale_PLM))
print("\t\t\tIM\t myeloma:\t epsilon_scale\t=\t%0.12f"%(epsilon_scale_IM_myeloma))
print("\t\t\tIM\t colon:\t epsilon_scale\t=\t%0.12f"%(epsilon_scale_IM_colon))
print("\t\t\tIM\t CML:\t epsilon_scale\t=\t%0.12f\n\n"%(epsilon_scale_IM_CML))
# Allocate four epsilon vectors with transformation parameters
epsilon_vector_PLM = np.linspace(0,2*epsilon_scale_PLM,200,endpoint=True)
epsilon_vector_IM_myeloma = np.linspace(0,epsilon_scale_IM_myeloma,200,endpoint=True)
epsilon_vector_IM_colon = np.linspace(0,epsilon_scale_IM_colon,200,endpoint=True)
epsilon_vector_IM_CML = np.linspace(0,epsilon_scale_IM_CML,200,endpoint=True)
#--------------------------------------------------------------------------------------------------------
# Create an iterable list
input_lists = [[t_myeloma,R_myeloma,epsilon_vector_PLM,PLM_fitted_to_myeloma_ODR,"PLM"],[t_myeloma,R_myeloma,epsilon_vector_IM_myeloma,IM_fitted_to_myeloma_ODR,"IM"],[t_colon,R_colon,epsilon_vector_PLM,PLM_fitted_to_colon_ODR,"PLM"],[t_colon,R_colon,epsilon_vector_IM_colon,IM_fitted_to_colon_ODR,"IM"],[t_CML,R_CML,epsilon_vector_PLM,PLM_fitted_to_CML_ODR,"PLM"],[t_CML,R_CML,epsilon_vector_IM_CML,IM_fitted_to_CML_ODR,"IM"]]
# We do the fitting in parallel (to use all cores, use mp.cpu_count())
pool = mp.Pool(mp.cpu_count())
# Do the model selection in parallel
results = pool.starmap(symmetry_toolbox.symmetry_based_model_selection,[(input_list[0],input_list[1],input_list[2],input_list[3],input_list[4]) for input_list in input_lists])
# Close the pool
pool.close()
# Extract the results
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Myeloma
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
print("\t\tThe output:")
print("\t\t\tMyeloma:")
# PLM MYELOMA
epsilon_transf_PLM_myeloma = results[0][0] # PLM transformation parameters
RMS_transf_PLM_myeloma = results[0][1] # PLM RMS-values
transformed_data_PLM_myeloma = results[0][2] # PLM Transformed data
fitted_parameters_PLM_myeloma = results[0][3] # PLM fitted parameters
inverse_parameters_PLM_myeloma = results[0][4] # PLM inversely transformed parameters
RMS_forward_PLM_myeloma = results[0][5] # RMS forward PLM myeloma
# Remove outliers
epsilon_transf_PLM_myeloma,RMS_transf_PLM_myeloma, fitted_parameters_PLM_myeloma,inverse_parameters_PLM_myeloma,transformed_data_PLM_myeloma = clean_up_output.remove_outliers(epsilon_transf_PLM_myeloma,RMS_transf_PLM_myeloma,fitted_parameters_PLM_myeloma,inverse_parameters_PLM_myeloma,transformed_data_PLM_myeloma)
# Save the data
write_output.save_data_symmetry_based_model_selection("myeloma","PLM",epsilon_transf_PLM_myeloma,RMS_transf_PLM_myeloma, fitted_parameters_PLM_myeloma,inverse_parameters_PLM_myeloma)
# Reduce density
if len(epsilon_transf_PLM_myeloma) > 200:
    epsilon_transf_PLM_myeloma_sparse = np.linspace(epsilon_transf_PLM_myeloma[0],epsilon_transf_PLM_myeloma[-1],200,endpoint=True)
    RMS_transf_PLM_myeloma_sparse = clean_up_output.reduce_density(epsilon_transf_PLM_myeloma,RMS_transf_PLM_myeloma,epsilon_transf_PLM_myeloma_sparse)
else:
    epsilon_transf_PLM_myeloma_sparse = epsilon_transf_PLM_myeloma
    RMS_transf_PLM_myeloma_sparse = RMS_transf_PLM_myeloma
#----------------------------------------------------------------------------------
# IM MYELOMA
epsilon_transf_IM_myeloma = results[1][0] # IM transformation parameters
RMS_transf_IM_myeloma = results[1][1] # IM RMS-values
transformed_data_IM_myeloma = results[1][2] # IM Transformed data
fitted_parameters_IM_myeloma = results[1][3] # IM fitted parameters
inverse_parameters_IM_myeloma = results[1][4] # IM inversely transformed parameters
RMS_forward_IM_myeloma = results[1][5] # RMS forward IM myeloma
# Remove outliers
epsilon_transf_IM_myeloma,RMS_transf_IM_myeloma,fitted_parameters_IM_myeloma,inverse_parameters_IM_myeloma,transformed_data_IM_myeloma = clean_up_output.remove_outliers(epsilon_transf_IM_myeloma,RMS_transf_IM_myeloma,fitted_parameters_IM_myeloma,inverse_parameters_IM_myeloma,transformed_data_IM_myeloma)
# Save the data
write_output.save_data_symmetry_based_model_selection("myeloma","IM",epsilon_transf_IM_myeloma,RMS_transf_IM_myeloma, fitted_parameters_IM_myeloma,inverse_parameters_IM_myeloma)
# Reduce density
if len(epsilon_transf_IM_myeloma) > 200:
    # Reduce the sparsity of these ridiculously dense vectors
    epsilon_transf_IM_myeloma_sparse = np.linspace(epsilon_transf_IM_myeloma[0],epsilon_transf_IM_myeloma[-1],200,endpoint=True)
    RMS_transf_IM_myeloma_sparse = clean_up_output.reduce_density(epsilon_transf_IM_myeloma,RMS_transf_IM_myeloma,epsilon_transf_IM_myeloma_sparse)
else:
    epsilon_transf_IM_myeloma_sparse = epsilon_transf_IM_myeloma
    RMS_transf_IM_myeloma_sparse = RMS_transf_IM_myeloma
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Colon
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
print("\t\t\tColon cancer:")
# PLM COLON
epsilon_transf_PLM_colon = results[2][0] # PLM transformation parameters
RMS_transf_PLM_colon = results[2][1] # PLM RMS-values
transformed_data_PLM_colon = results[2][2] # PLM Transformed data
fitted_parameters_PLM_colon = results[2][3] # PLM fitted parameters
inverse_parameters_PLM_colon = results[2][4] # PLM inversely transformed parameters
RMS_forward_PLM_colon = results[2][5] # RMS forward PLM colon
# Remove outliers
epsilon_transf_PLM_colon,RMS_transf_PLM_colon,fitted_parameters_PLM_colon,inverse_parameters_PLM_colon,transformed_data_PLM_colon = clean_up_output.remove_outliers(epsilon_transf_PLM_colon,RMS_transf_PLM_colon,fitted_parameters_PLM_colon,inverse_parameters_PLM_colon,transformed_data_PLM_colon)
# Save the data
write_output.save_data_symmetry_based_model_selection("colon","PLM",epsilon_transf_PLM_colon,RMS_transf_PLM_colon, fitted_parameters_PLM_colon,inverse_parameters_PLM_colon)
# Reduce density
if len(epsilon_transf_PLM_colon) > 200:
    # Reduce the sparsity of these ridiculously dense vectors
    epsilon_transf_PLM_colon_sparse = np.linspace(epsilon_transf_PLM_colon[0],epsilon_transf_PLM_colon[-1],200,endpoint=True)
    RMS_transf_PLM_colon_sparse = clean_up_output.reduce_density(epsilon_transf_PLM_colon,RMS_transf_PLM_colon,epsilon_transf_PLM_colon_sparse)
else:
    epsilon_transf_PLM_colon_sparse = epsilon_transf_PLM_colon
    RMS_transf_PLM_colon_sparse = RMS_transf_PLM_colon
#----------------------------------------------------------------------------------
# IM COLON
epsilon_transf_IM_colon = results[3][0] # IM transformation parameters
RMS_transf_IM_colon = results[3][1] # IM RMS-values
transformed_data_IM_colon = results[3][2] # IM Transformed data
fitted_parameters_IM_colon = results[3][3] # IM fitted parameters
inverse_parameters_IM_colon = results[3][4] # IM inversely transformed parameters
RMS_forward_IM_colon = results[3][5] # RMS forward IM colon
# Remove outliers
epsilon_transf_IM_colon,RMS_transf_IM_colon,fitted_parameters_IM_colon,inverse_parameters_IM_colon,transformed_data_IM_colon = clean_up_output.remove_outliers(epsilon_transf_IM_colon,RMS_transf_IM_colon,fitted_parameters_IM_colon,inverse_parameters_IM_colon,transformed_data_IM_colon)
# Save the data
write_output.save_data_symmetry_based_model_selection("colon","IM",epsilon_transf_IM_colon,RMS_transf_IM_colon, fitted_parameters_IM_colon,inverse_parameters_IM_colon)
# Reduce density
if len(epsilon_transf_IM_colon) > 200:
    # Reduce the sparsity of these ridiculously dense vectors
    epsilon_transf_IM_colon_sparse = np.linspace(epsilon_transf_IM_colon[0],epsilon_transf_IM_colon[-1],200,endpoint=True)
    RMS_transf_IM_colon_sparse = clean_up_output.reduce_density(epsilon_transf_IM_colon,RMS_transf_IM_colon,epsilon_transf_IM_colon_sparse)
else:
    epsilon_transf_IM_colon_sparse = epsilon_transf_IM_colon
    RMS_transf_IM_colon_sparse = RMS_transf_IM_colon
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# CML
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
print("\t\t\tCML:")
# PLM CML
epsilon_transf_PLM_CML = results[4][0] # PLM transformation parameters
RMS_transf_PLM_CML = results[4][1] # PLM RMS-values
transformed_data_PLM_CML = results[4][2] # PLM Transformed data
fitted_parameters_PLM_CML = results[4][3] # PLM fitted parameters
inverse_parameters_PLM_CML = results[4][4] # PLM inversely transformed parameters
RMS_forward_PLM_CML = results[4][5] # RMS forward PLM CML
# Remove outliers
epsilon_transf_PLM_CML,RMS_transf_PLM_CML,fitted_parameters_PLM_CML,inverse_parameters_PLM_CML,transformed_data_PLM_CML = clean_up_output.remove_outliers(epsilon_transf_PLM_CML,RMS_transf_PLM_CML,fitted_parameters_PLM_CML,inverse_parameters_PLM_CML,transformed_data_PLM_CML)
# Save the data
write_output.save_data_symmetry_based_model_selection("CML","PLM",epsilon_transf_PLM_CML,RMS_transf_PLM_CML, fitted_parameters_PLM_CML,inverse_parameters_PLM_CML)
# Reduce density
if len(epsilon_transf_PLM_CML)>200:
    # Reduce the sparsity of these ridiculously dense vectors
    epsilon_transf_PLM_CML_sparse = np.linspace(epsilon_transf_PLM_CML[0],epsilon_transf_PLM_CML[-1],200,endpoint=True)
    RMS_transf_PLM_CML_sparse = clean_up_output.reduce_density(epsilon_transf_PLM_CML,RMS_transf_PLM_CML,epsilon_transf_PLM_CML_sparse)
else:
    epsilon_transf_PLM_CML_sparse = epsilon_transf_PLM_CML
    RMS_transf_PLM_CML_sparse = RMS_transf_PLM_CML
#----------------------------------------------------------------------------------
# IM CML
epsilon_transf_IM_CML = results[5][0] # IM
RMS_transf_IM_CML = results[5][1] # IM
transformed_data_IM_CML = results[5][2] # IM Transformed data
fitted_parameters_IM_CML = results[5][3] # IM fitted parameters
inverse_parameters_IM_CML = results[5][4] # IM inversely transformed parameters
RMS_forward_IM_CML = results[5][5] # IM CML
# Remove outliers
epsilon_transf_IM_CML,RMS_transf_IM_CML,fitted_parameters_IM_CML,inverse_parameters_IM_CML,transformed_data_IM_CML = clean_up_output.remove_outliers(epsilon_transf_IM_CML,RMS_transf_IM_CML,fitted_parameters_IM_CML,inverse_parameters_IM_CML,transformed_data_IM_CML)
# Save the data
write_output.save_data_symmetry_based_model_selection("CML","IM",epsilon_transf_IM_CML,RMS_transf_IM_CML, fitted_parameters_IM_CML,inverse_parameters_IM_CML)
# Reduce density
if len(epsilon_transf_IM_CML) > 200:
    # Reduce the sparsity of these ridiculously dense vectors
    epsilon_transf_IM_CML_sparse = np.linspace(epsilon_transf_IM_CML[0],epsilon_transf_IM_CML[-1],200,endpoint=True)
    RMS_transf_IM_CML_sparse = clean_up_output.reduce_density(epsilon_transf_IM_CML,RMS_transf_IM_CML,epsilon_transf_IM_CML_sparse)
else:
    epsilon_transf_IM_CML_sparse = epsilon_transf_IM_CML
    RMS_transf_IM_CML_sparse = RMS_transf_IM_CML
#----------------------------------------------------------------------------------    
#----------------------------------------------------------------------------------
# Prompt to the user
print("\t\tSymmetry framework is done!")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the symmetry based model selection
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Proper framework where we do all four steps:
# 1. Transform data,
# 2. Fit to transformed data,
# 3. Inversely transform model back,
# 4. Compare fit to original data.
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: PLM
axes[0].plot(epsilon_transf_PLM_myeloma,RMS_transf_PLM_myeloma_sparse,'-', color = (103/256,0/256,31/256),label='PLM Myeloma cancer')
axes[0].plot(epsilon_transf_PLM_colon_sparse,RMS_transf_PLM_colon_sparse,'-', color = (206/256,18/256,86/256),label='PLM Colon cancer')
axes[0].plot(epsilon_transf_PLM_CML_sparse,RMS_transf_PLM_CML_sparse,'-', color = (223/256,101/256,176/256),label='PLM CML')
axes[0].legend()
# Subplot 2: IM
axes[1].plot(epsilon_transf_IM_myeloma_sparse,RMS_transf_IM_myeloma_sparse,'-', color = (2/256,56/256,88/256),label='IM Myeloma cancer')
axes[1].plot(epsilon_transf_IM_colon_sparse,RMS_transf_IM_colon_sparse,'-', color = (54/256,144/256,192/256),label='IM Colon cancer')
axes[1].plot(epsilon_transf_IM_CML_sparse,RMS_transf_IM_CML_sparse,'-', color = (208/256,209/256,230/256),label='IM CML')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("The symmetry based model selection",fontsize=20, fontweight='bold')
plt.savefig("../Figures/symmetry_based_model_selection.png")
#-----------------------------------------------------------------------------------------------------------------------------
# Cheating where we only do forward fittig. I.e. we only do two steps:
# 1. Transform data,
# 2. Fit to transformed data.
# Overall properties
fig, axes = plt.subplots(1,2,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1: PLM
axes[0].plot(results[0][0],RMS_forward_PLM_myeloma,'-', color = (103/256,0/256,31/256),label='PLM Myeloma cancer')
axes[0].plot(results[2][0],RMS_forward_PLM_colon,'-', color = (206/256,18/256,86/256),label='PLM Colon cancer')
axes[0].plot(results[4][0],RMS_forward_PLM_CML,'-', color = (223/256,101/256,176/256),label='PLM CML')
axes[0].legend()
# Subplot 2: IM
axes[1].plot(results[1][0],RMS_forward_IM_myeloma,'-', color = (2/256,56/256,88/256),label='IM Myeloma cancer')
axes[1].plot(results[3][0],RMS_forward_IM_colon,'-', color = (54/256,144/256,192/256),label='IM Colon cancer')
axes[1].plot(results[5][0],RMS_forward_IM_CML,'-', color = (208/256,209/256,230/256),label='IM CML')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Transformation parameter, $\epsilon$")
plt.ylabel("Root mean square, $\mathrm{RMS}(\epsilon)$")
# displaying the title
plt.title("Fitting to transformed data",fontsize=20, fontweight='bold')
plt.savefig("../Figures/forward_fitting.png")
# =================================================================================
# =================================================================================
# ILLUSTRATE THE MODEL SELECTION FRAMEWORK
# =================================================================================
# =================================================================================
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# PLM
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
print("\t\tIllustrate framework!\n")
#------------------------------------------------------------------------------
# STEP 1 OUT OF 4: TRANSFORM DATA
#------------------------------------------------------------------------------
# Find the index closest to the epsilon scale of the PLM
index_PLM_colon = np.where(epsilon_transf_PLM_colon==epsilon_transf_PLM_colon[np.abs(epsilon_transf_PLM_colon-(epsilon_scale_PLM)).argmin()])[0][0]
# Extract out epsilon value
epsilon = epsilon_transf_PLM_colon[index_PLM_colon]
print("\t\tPLM\n")
print("\t\t\tPLM epsilon scale:\t%0.3f"%(epsilon))
# Allocate memory for a list
R_PLM_trans = []
t_PLM_trans = []
# Allocate an index vector
index_vector = list(np.arange(10,len(t_colon),2))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_colon[index],R_colon[index],0.9*epsilon,fitted_parameters_PLM_colon[0][1])
    # Save the transformed variables
    R_PLM_trans.append(R_trans)
    t_PLM_trans.append(t_trans)
# Transform the data
t_colon_trans_PLM,R_colon_trans_PLM = transformed_data_PLM_colon[index_PLM_colon]
#------------------------------------------------------------------------------
# STEP 2 OUT OF 4: FIT MODEL TO TRANSFORMED DATA
#------------------------------------------------------------------------------
# Fit the model to the transformed data
R_hat_PLM_trans = np.array([fit_to_data.objective_PLM(fitted_parameters_PLM_colon[index_PLM_colon], t_i) for t_i in list(t_colon_trans_PLM)])
#------------------------------------------------------------------------------
# STEP 3 OUT OF 4: INVERSELY TRANSFORM FITTED MODEL BACK
#------------------------------------------------------------------------------
# Allocate memory for a list
R_PLM_trans_inv = []
t_PLM_trans_inv = []
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_colon_trans_PLM[index],R_hat_PLM_trans[index],-0.9*epsilon,inverse_parameters_PLM_colon[index_PLM_colon][1])
    # Save the transformed variables
    R_PLM_trans_inv.append(R_trans)
    t_PLM_trans_inv.append(t_trans)
# Calculate inversely transformed solution
R_hat_PLM_inv_trans = np.array([fit_to_data.objective_PLM(inverse_parameters_PLM_colon[index_PLM_colon], t_i) for t_i in list(t_colon)])
#------------------------------------------------------------------------------
# IM
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# STEP 1 OUT OF 4: TRANSFORM DATA
#------------------------------------------------------------------------------
# Find the index closest to the epsilon scale of the PLM
index_IM_colon = np.where(epsilon_transf_IM_colon==epsilon_transf_IM_colon[np.abs(epsilon_transf_IM_colon-epsilon_scale_IM_colon).argmin()])[0][0]
# Extract an epsilon value
epsilon = epsilon_transf_IM_colon[index_IM_colon]
print("\t\tIM\n")
print("\t\t\tIM epsilon scale:\t%0.3f"%(epsilon))
# Allocate memory for a list
R_IM_trans = []
t_IM_trans = []
# Allocate an index vector
index_vector = list(np.arange(60,len(t_colon),1))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_transformation(t_colon[index],R_colon[index],epsilon*0.95,fitted_parameters_IM_colon[0][1],fitted_parameters_IM_colon[0][3])
    # Save the transformed variables
    R_IM_trans.append(R_trans)
    t_IM_trans.append(t_trans)
# Transform the data
t_colon_trans_IM,R_colon_trans_IM = transformed_data_IM_colon[index_IM_colon]
#------------------------------------------------------------------------------
# STEP 2 OUT OF 4: FIT MODEL TO TRANSFORMED DATA
#------------------------------------------------------------------------------
t_sym_IM_trans = np.linspace(t_colon_trans_IM[0],t_colon_trans_IM[len(t_colon_trans_IM)-1],100, endpoint=True)
R_hat_IM_trans = np.array([fit_to_data.objective_IM(fitted_parameters_IM_colon[index_IM_colon], t_i) for t_i in list(t_sym_IM_trans)])
#------------------------------------------------------------------------------
# STEP 3 OUT OF 4: INVERSELY TRANSFORM FITTED MODEL BACK
#------------------------------------------------------------------------------
# Allocate memory for a list
R_IM_trans_inv = []
t_IM_trans_inv = []
# Allocate an index vector
index_vector_1 = list(np.arange(len(t_sym_IM_trans)-50,len(t_sym_IM_trans)-20,4))
index_vector_2 = list(np.arange(len(t_sym_IM_trans)-20,len(t_sym_IM_trans)-1,8))
index_vector = index_vector_1 + index_vector_2
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_transformation(t_sym_IM_trans[index],R_hat_IM_trans[index],-epsilon*0.95,inverse_parameters_IM_colon[index_IM_colon][1],inverse_parameters_IM_colon[index_IM_colon][3])
    # Save the transformed variables
    R_IM_trans_inv.append(R_trans)
    t_IM_trans_inv.append(t_trans)
# Calculate inversely transformed solution
t_sym_IM_inv_trans = np.linspace(t_colon[0],t_colon[len(t_colon)-1]+2,100, endpoint=True)
R_hat_IM_inv_trans = np.array([fit_to_data.objective_IM(inverse_parameters_IM_colon[index_IM_colon], t_i) for t_i in list(t_sym_IM_inv_trans)])
#----------------------------------------------------------------------------------
# Plot the illustration of the model selection framework
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Illustration of the framework
# Overall properties
fig, axes = plt.subplots(2,4,figsize=(15,5))
plt.rc('axes', labelsize=15)    # fontsize of the x and y label
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
# Subplot 1a
axes[0][0].plot(t_colon, R_colon, '*', color='black', label='Data Colon cancer') # Original data
axes[0][0].plot(t_colon, R_hat_PLM_colon_ODR, '-', color = (103/256,0/256,31/256),label='Original fit PLM')# Original fitted model
axes[0][0].plot(np.array(t_PLM_trans[0]),np.array(R_PLM_trans[0]),'--',color='black',label='Symmetry PLM') # Original symmetry
for index in range(1,len(t_PLM_trans)):
    axes[0][0].plot(np.array(t_PLM_trans[index]),np.array(R_PLM_trans[index]),'--',color='black')
axes[0][0].plot(t_colon_trans_PLM, R_colon_trans_PLM, '*', color='gray', label='Transformed Data') # Transformed data
axes[0][0].legend()
# Subplot 1b
axes[0][1].plot(t_colon_trans_PLM, R_colon_trans_PLM, '*', color='gray', label='Transformed Data') # Transformed data
axes[0][1].plot(t_colon_trans_PLM, R_hat_PLM_trans, '-', color = (206/256,18/256,86/256), label='Fitted PLM') # Fitted model
axes[0][1].legend()
# Subplot 1c
axes[0][2].plot(t_colon_trans_PLM, R_hat_PLM_trans, '-', color=(206/256,18/256,86/256), label='Fitted PLM') # Fitted model
axes[0][2].plot(np.array(t_PLM_trans_inv[0]),np.array(R_PLM_trans_inv[0]),'--',color='black',label='Inverse symmetry PLM') # Inverse transform
for index in range(1,len(t_PLM_trans_inv)):
    axes[0][2].plot(np.array(t_PLM_trans_inv[index]),np.array(R_PLM_trans_inv[index]),'--',color='black')
axes[0][2].plot(t_colon,R_hat_PLM_inv_trans, '-', color = (223/256,101/256,176/256), label='Inversely transformed PLM') # Inverse model
axes[0][2].legend()
# Subplot 1d
axes[0][3].plot(t_colon, R_colon, '*', color='black', label='Data Colon cancer') # Original data
axes[0][3].plot(t_colon, R_hat_PLM_colon_ODR, '-', color = (103/256,0/256,31/256),label='Original PLM')
axes[0][3].plot(t_colon,R_hat_PLM_inv_trans, '-', color = (223/256,101/256,176/256), label='Inversely transformed PLM') # Inverse model
axes[0][3].legend()
# Subplot 2a
axes[1][0].plot(t_colon, R_colon, '*', color='black', label='Data Colon cancer')
axes[1][0].plot(np.array(t_IM_trans[0]),np.array(R_IM_trans[0]),'--',color='black',label='Symmetry PLM')
for index in range(1,len(t_IM_trans)):
    axes[1][0].plot(np.array(t_IM_trans[index]),np.array(R_IM_trans[index]),'--',color='black')
axes[1][0].plot(t_colon_trans_IM, R_colon_trans_IM, '*', color='gray', label='Transformed Data')
axes[1][0].plot(t_colon, R_hat_IM_colon_ODR, '-', color = (2/256,56/256,88/256),label='Original fit IM')
axes[1][0].legend()
# Subplot 1b
axes[1][1].plot(t_colon_trans_IM, R_colon_trans_IM, '*', color='gray', label='Transformed Data')
axes[1][1].plot(t_sym_IM_trans, R_hat_IM_trans, '-', color = (54/256,144/256,192/256), label='Fitted IM')
axes[1][1].legend()
# Subplot 1c
axes[1][2].plot(t_sym_IM_trans, R_hat_IM_trans, '-', color = (54/256,144/256,192/256), label='Fitted IM')
axes[1][2].plot(np.array(t_IM_trans_inv[0]),np.array(R_IM_trans_inv[0]),'--',color='black',label='Inverse symmetry IM')
for index in range(1,len(t_IM_trans_inv)):
    axes[1][2].plot(np.array(t_IM_trans_inv[index]),np.array(R_IM_trans_inv[index]),'--',color='black')
axes[1][2].plot(t_sym_IM_inv_trans,R_hat_IM_inv_trans, '-', color = (208/256,209/256,230/256), label='Inversely transformed IM')
axes[1][2].legend()
# Subplot 1d
axes[1][3].plot(t_colon, R_colon, '*', color='black', label='Data Colon cancer')
axes[1][3].plot(t_sym_IM_inv_trans,R_hat_IM_inv_trans, '-', color = (208/256,209/256,230/256), label='Inversely transformed IM')
axes[1][3].plot(t_colon, R_hat_IM_colon_ODR, '-', color = (2/256,56/256,88/256),label='Original fit IM')
axes[1][3].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The steps of the symmetry based model selection",fontsize=20, fontweight='bold')
plt.savefig("../Figures/step_symmetry_based_model_selection.png")
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
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_sym[index],R_hat_PLM_original[index],epsilon*0.95,PLM_fitted_to_myeloma_ODR.beta[1])
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
    t_trans,R_trans = symmetry_toolbox.PLM_transformation(t_hat_PLM_1[index],R_hat_PLM_1[index],epsilon*0.95,PLM_fitted_to_myeloma_ODR.beta[1])
    # Save the transformed variables
    R_PLM_trans_2.append(R_trans)
    t_PLM_trans_2.append(t_trans)
# Transform the transformed solution
t_hat_PLM_2,R_hat_PLM_2 = symmetry_toolbox.PLM_transformed_solution(t_hat_PLM_1,R_hat_PLM_1,epsilon,PLM_fitted_to_myeloma_ODR.beta[0],PLM_fitted_to_myeloma_ODR.beta[1])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
epsilon = epsilon_scale_IM_myeloma/2
# THE IM-II
# Construct a t vector
t_sym = np.linspace(0,t_myeloma[len(t_myeloma)-1],200)
# Original solution
R_hat_IM_original = np.array([fit_to_data.objective_IM(IM_fitted_to_myeloma_ODR.beta,t_sym[index]) for index in range(len(t_sym))])
# Allocate memory for a list
R_IM_trans_1 = []
t_IM_trans_1 = []
# Allocate an index vector
index_vector = list(np.arange(193,len(t_sym)-1,2))
#index_vector = list(np.arange(190,len(t_sym)-1,2))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_transformation(t_sym[index],R_hat_IM_original[index],epsilon*0.9,IM_fitted_to_myeloma_ODR.beta[1],IM_fitted_to_myeloma_ODR.beta[3])
    # Save the transformed variables
    R_IM_trans_1.append(R_trans)
    t_IM_trans_1.append(t_trans)
# Transform the original solution
t_hat_IM_1,R_hat_IM_1 = symmetry_toolbox.IM_transformed_solution(t_sym,R_hat_IM_original,epsilon,IM_fitted_to_myeloma_ODR.beta[1],IM_fitted_to_myeloma_ODR.beta[3])
# Allocate memory for a list
R_IM_trans_2 = []
t_IM_trans_2 = []
# Update the index vector
index_vector = list(np.arange(190,len(t_sym)-1,2))
#index_vector = list(np.arange(190,len(t_sym)-1,1))
# Save all the transformed stuff
for index in index_vector:
    # Transform stuff
    t_trans,R_trans = symmetry_toolbox.IM_transformation(t_hat_IM_1[index],R_hat_IM_1[index],epsilon*0.9,IM_fitted_to_myeloma_ODR.beta[1],IM_fitted_to_myeloma_ODR.beta[3])
    # Save the transformed variables
    R_IM_trans_2.append(R_trans)
    t_IM_trans_2.append(t_trans)
# Transform the second solution      
t_hat_IM_2,R_hat_IM_2 = symmetry_toolbox.IM_transformed_solution(t_hat_IM_1,R_hat_IM_1,epsilon,IM_fitted_to_myeloma_ODR.beta[1],IM_fitted_to_myeloma_ODR.beta[3])
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tAction of symmetries!\n")
print("\n\t--------------------------------------------------------------------------------------\n")
print("\t\tEpsilon scales for illustrations (i.e. the value that the transformations are pushed with twice)")
print("\t\tPLM,\t%0.7f"%(epsilon_scale_PLM/2))
print("\t\tIM,\t%0.7f"%(epsilon_scale_IM_myeloma/2))
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the action of the symmetries in Python
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
axes[1].plot(t_sym, R_hat_IM_original, '-', color = (2/256,56/256,88/256), label='Original solution IM')
axes[1].plot(np.array(t_IM_trans_1[0]),np.array(R_IM_trans_1[0]),'--',color='black',label='Symmetry IM')
for index in range(1,len(t_IM_trans_1)):
    axes[1].plot(np.array(t_IM_trans_1[index]),np.array(R_IM_trans_1[index]),'--',color='black')
axes[1].plot(t_hat_IM_1, R_hat_IM_1, '-', color=(54/256,144/256,192/256), label='IM-II transformed solution 1')
for index in range(len(t_IM_trans_2)):
    axes[1].plot(np.array(t_IM_trans_2[index]),np.array(R_IM_trans_2[index]),'--',color='black')
axes[1].plot(t_hat_IM_2, R_hat_IM_2, '-', color=(208/256,209/256,230/256), label='IM transformed solution 2')
axes[1].legend()
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("The symmetries of the PLM and the IM",fontsize=20, fontweight='bold')
plt.savefig("../Figures/action_of_symmetries.png")
plt.show()
# Final prompt to the user
print("\n\t--------------------------------------------------------------------------------------\n")
print("\n\t\tCalculations are done!\n")
print("\n\t--------------------------------------------------------------------------------------\n")    

    
