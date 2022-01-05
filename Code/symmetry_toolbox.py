# =================================================================================
# =================================================================================
# Script:"symmetry_based_model_selection"
# Date: 2021-09-20
# Implemented by: Johannes Borgqvist
# Description:
# The script involves all the symmetry functionalities
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import read_data  # Home-made
import write_output  # Home-made
import fit_to_data # Home-made
import numpy as np 
from matplotlib import pyplot as plt
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# FUNCTION 1: PLM_1_symmetry
# The first symmetry of the power law model (PLM)
def PLM_1_symmetry(t,R,epsilon,gamma):
    t_hat = t
    R_hat = R + epsilon*(t**gamma)
    return t_hat,R_hat
# FUNCTION 2: PLM_2_symmetry
# The second symmetry of the power law model (PLM)
def PLM_2_symmetry(t,R,epsilon):
    t_hat = t * np.exp(epsilon)
    R_hat = R
    return t_hat,R_hat
# FUNCTION 3: PLM_transformation
# The function plots the action of
# the symmetries of the PLM from
# a given starting point and an end point
def PLM_transformation(t_from,R_from,epsilon,gamma,indicator):
    # In order to plot the action of the transformation
    # continuously we will add small plots along the
    # way.
    epsilon_increment = epsilon/10 # Increment rotation angle
    epsilon_temp = epsilon_increment # The current angle we plot with
    # Two lists as output corresponding to the transformation
    t_trans_list = [] # Transformed variable
    R_trans_list = [] # Transformed state
    # Append the starting point
    t_trans_list.append(t_from)
    R_trans_list.append(R_from)
    # Calculate the transformed point of interest
    if indicator == 1:
        # Calculate the next point for the transformation
        t_to, R_to = PLM_1_symmetry(t_from,R_from,epsilon_increment,gamma)
    elif indicator == 2:
        # Calculate the next point for the transformation
        t_to, R_to = PLM_2_symmetry(t_from,R_from,epsilon_increment)        
    # Save the new rotated point
    t_trans_list.append(t_to)
    R_trans_list.append(R_to)
    # Re-do this procedure until we have rotated the original point
    # far enough
    while abs(epsilon_temp)<abs(epsilon):
        # Update point we want to rotate
        t_from = t_to
        R_from = R_to
        # Calculate the transformed point of interest
        if indicator == 1:
            # Next point for the transformation
            t_to, R_to = PLM_1_symmetry(t_from,R_from,epsilon_increment,gamma)
        elif indicator == 2:
            # Calculate the next point for the transformation
            t_to, R_to = PLM_2_symmetry(t_from,R_from,epsilon_increment)                    
        # Save the rotated point
        t_trans_list.append(t_to)
        R_trans_list.append(R_to)
        # Change the angle we rotate with
        epsilon_temp += epsilon_increment
    # Convert the lists to arrays
    t_trans = np.asarray(t_trans_list)
    R_trans = np.asarray(R_trans_list)
    # Return the action of the transformation
    return t_trans,R_trans
# FUNCTION 4: IM_II_symmetry
# The symmetry of the immunological model IM-II
def IM_II_symmetry(t,R,epsilon,tau,alpha):
    # Define t_hat recursively
    t_hat =  np.log(np.log(np.abs((alpha*np.exp(alpha*tau)*epsilon) - np.exp(np.exp(-alpha*(t-tau))))))
    t_hat = tau - ((t_hat)/(alpha))
    R_hat = R
    return t_hat,R_hat
# FUNCTION 5: PLM_transformation
# The function plots the action of
# the symmetries of the PLM from
# a given starting point and an end point
def IM_II_transformation(t_from,R_from,epsilon,tau,alpha):
    # In order to plot the action of the transformation
    # continuously we will add small plots along the
    # way.
    epsilon_increment = epsilon/10 # Increment rotation angle
    epsilon_temp = epsilon_increment # The current angle we plot with
    # Two lists as output corresponding to the transformation
    t_trans_list = [] # Transformed variable
    R_trans_list = [] # Transformed state
    # Append the starting point
    t_trans_list.append(t_from)
    R_trans_list.append(R_from)
    # Calculate the transformed point of interest
    t_to, R_to = IM_II_symmetry(t_from,R_from,epsilon_increment,tau,alpha)
    # Save the new rotated point
    t_trans_list.append(t_to)
    R_trans_list.append(R_to)
    # Re-do this procedure until we have rotated the original point
    # far enough
    while abs(epsilon_temp)<abs(epsilon):
        # Update point we want to rotate
        t_from = t_to
        R_from = R_to
        # Calculate the next transformed point
        t_to, R_to = IM_II_symmetry(t_from,R_from,epsilon_increment,tau,alpha)
        # Save the rotated point
        t_trans_list.append(t_to)
        R_trans_list.append(R_to)
        # Change the angle we rotate with
        epsilon_temp += epsilon_increment
    # Convert the lists to arrays
    t_trans = np.asarray(t_trans_list)
    R_trans = np.asarray(R_trans_list)
    # Return the action of the transformation
    return t_trans,R_trans
# FUNCTION 6: IM_I_1_symmetry
# The first symmetry of the exponential model (IM-I)
def IM_I_1_symmetry(t,R,epsilon,alpha):
    t_hat = t
    R_hat = R + epsilon*(np.exp(alpha*t))
    return t_hat,R_hat
# FUNCTION 7: IM_I_2_symmetry
# The second symmetry of the exponential model (IM-I)
def IM_I_2_symmetry(t,R,epsilon):
    t_hat = t + epsilon
    R_hat = R
    return t_hat,R_hat
# FUNCTION 8: PLM_transformation
# The function plots the action of
# the symmetries of the PLM from
# a given starting point and an end point
def IM_I_transformation(t_from,R_from,epsilon,alpha,indicator):
    # In order to plot the action of the transformation
    # continuously we will add small plots along the
    # way.
    epsilon_increment = epsilon/10 # Increment rotation angle
    epsilon_temp = epsilon_increment # The current angle we plot with
    # Two lists as output corresponding to the transformation
    t_trans_list = [] # Transformed variable
    R_trans_list = [] # Transformed state
    # Append the starting point
    t_trans_list.append(t_from)
    R_trans_list.append(R_from)
    # Calculate the transformed point of interest
    if indicator == 1:
        # Calculate the next point for the transformation
        t_to, R_to = IM_I_1_symmetry(t_from,R_from,epsilon_increment,alpha)
    elif indicator == 2:
        # Calculate the next point for the transformation
        t_to, R_to = IM_I_2_symmetry(t_from,R_from,epsilon_increment)        
    # Save the new rotated point
    t_trans_list.append(t_to)
    R_trans_list.append(R_to)
    # Re-do this procedure until we have rotated the original point
    # far enough
    while abs(epsilon_temp)<abs(epsilon):
        # Update point we want to rotate
        t_from = t_to
        R_from = R_to
        # Calculate the transformed point of interest
        if indicator == 1:
            # Next point for the transformation
            t_to, R_to = IM_I_1_symmetry(t_from,R_from,epsilon_increment,alpha)
        elif indicator == 2:
            # Calculate the next point for the transformation
            t_to, R_to = IM_I_2_symmetry(t_from,R_from,epsilon_increment)                    
        # Save the rotated point
        t_trans_list.append(t_to)
        R_trans_list.append(R_to)
        # Change the angle we rotate with
        epsilon_temp += epsilon_increment
    # Convert the lists to arrays
    t_trans = np.asarray(t_trans_list)
    R_trans = np.asarray(R_trans_list)
    # Return the action of the transformation
    return t_trans,R_trans
# FUNCTION 9: sym_model_sel
# The function conducts the symmetry based model selection.
# It takes the following input:
# 1. The time in the array t,
# 2. The incidences in the array R, 
# 3. The epsilon_range defining for how many parameters we will do the fitting,
# 4. The optimal parameters in the vector opt_para,
# 5. The model_str defining which model we work with. 
# It returns the following output:
# 1. "Delta_fit" being the relative R^2_adj which has the same
# length as the "epsilon_range" value. 
def sym_model_sel(t,R,epsilon_vector,opt_para,model_str):
    # SETTING UP THE MODEL SELECTION
    # Allocate memory for the output
    Delta_fit = []
    # Calculate the original fit
    R_hat, opt_para_original, R_squared_adj_original = fit_to_data.PE_risk_profiles(t,R,model_str,opt_para)
    # Re-define data by taking the exponential of it
    R_ori_data = [np.exp(R[i]) for i in range(len(R))]
    # Extract the value of alpha
    alpha = opt_para[0]
    # Extract the other parameters
    if model_str == "PLM":
        gamma = opt_para[1]
    elif model_str == "IM-II":
        tau = opt_para[1]
    #====================================================================================
    # Loop over our epsilon values and conduct the model selection procedure
    for epsilon in epsilon_vector:
        # ----------------------------------------------------
        # STEP 1: Transform the data
        # ----------------------------------------------------        
        # Define the transformed data
        R_trans_data = []
        t_trans_data = []
        # Loop over the original time series and calculate the transformed time series
        for index in range(len(R_ori_data)):
            # Calculate each transformed coordinate
            if model_str == "PLM":
                t_hat_temp, R_hat_temp = PLM_2_symmetry(t[index],R_ori_data[index],epsilon)
            elif model_str == "IM-II":
                t_hat_temp, R_hat_temp = IM_II_symmetry(t[index],R_ori_data[index],epsilon,tau,alpha)
            # Save the transformed coordinates
            t_trans_data.append(t_hat_temp)
            R_trans_data.append(R_hat_temp)
        # ----------------------------------------------------
        # STEP 2: Fit model to the transformed data
        # ----------------------------------------------------
        # Define the logarithm of the data 
        R_trans_data_log = [np.log(R_trans) for R_trans in R_trans_data]
        # Convert the transformed time series to arrays
        t_trans_data = np.array(t_trans_data)
        R_trans_data_log = np.array(R_trans_data_log)
        # Calculate the fit
        R_hat_trans, opt_para_temp, R_squared_adj_after = fit_to_data.PE_risk_profiles(t_trans_data,R_trans_data_log,model_str,opt_para)
        # Append the relative fit
        #Delta_fit.append(((R_squared_adj_after)/(R_squared_adj_original)))
        Delta_fit.append(R_squared_adj_after)
    #====================================================================================        
    # Return the relative R^2_adj value as a function of epsilon    
    return Delta_fit

