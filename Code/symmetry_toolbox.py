# =================================================================================
# =================================================================================
# Script:"symmetry_based_model_selection"
# Date: 2022-02-11
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
#----------------------------------------------------------------------------------
# FUNCTION 1: PLM_symmetry
# The second symmetry of the power law model (PLM)
def PLM_symmetry(t,R,epsilon):
    t_hat = t * np.exp(epsilon)
    R_hat = R
    return t_hat,R_hat
#----------------------------------------------------------------------------------
# FUNCTION 2: PLM_transformation
# The function plots the action of
# the symmetries of the PLM from
# a given starting point and an end point
def PLM_transformation(t_from,R_from,epsilon,gamma):
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
    t_to, R_to = PLM_symmetry(t_from,R_from,epsilon_increment)        
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
        t_to, R_to = PLM_symmetry(t_from,R_from,epsilon_increment)                    
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
#----------------------------------------------------------------------------------
# FUNCTION 3: PLM_transformed_solution
# This function returns the transformed solution of the PLM for a given value
# of the parameters A and gamma as well as the transformation parameter epsilon
def PLM_transformed_solution(t,R,epsilon,A,gamma):
    # Allocate memory for the transformed solution
    #t_trans = t 
    t_trans = np.asarray([PLM_symmetry(t[index],R[index],epsilon)[0] for index in range(len(t))])
    #R_trans = np.asarray([fit_to_data.objective_PLM((A*np.exp(-gamma*epsilon),gamma),t_i) for t_i in list(t_trans)])
    R_trans = R
    # Return the transformed solution
    return t_trans,R_trans
#----------------------------------------------------------------------------------
# FUNCTION 4: IM_III_symmetry
# The symmetry of the immunological model IM-III
def IM_III_symmetry(t,R,epsilon,tau,alpha):
    # Define t_hat recursively
    t_hat =  np.log(np.log(np.exp(np.exp(-alpha*(t-tau))) - (alpha*np.exp(alpha*tau)*epsilon) ))
    t_hat = tau - ((t_hat)/(alpha))
    R_hat = R
    return t_hat,R_hat
#----------------------------------------------------------------------------------
# FUNCTION 5: IM_III_transformation
# The function plots the action of
# the symmetries of the IM-III from
# a given starting point and an end point
def IM_III_transformation(t_from,R_from,epsilon,tau,alpha):
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
    t_to, R_to = IM_III_symmetry(t_from,R_from,epsilon_increment,tau,alpha)
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
        t_to, R_to = IM_III_symmetry(t_from,R_from,epsilon_increment,tau,alpha)
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
#----------------------------------------------------------------------------------
# FUNCTION 8: IM_III_transformed_solution
# This function returns the transformed solution of the IM-III for a given value
# of the parameters A, tau and C as well as the transformation parameter epsilon
def IM_III_transformed_solution(t,R,epsilon,A,tau,C):
    # Define alpha to its fixed value
    alpha = 0.044
    # Allocate memory for the transformed solution
    #t_trans = t
    t_trans = np.asarray([IM_III_symmetry(t[index],R[index],epsilon,tau,alpha)[0] for index in range(len(t))])
    #R_trans = np.asarray([fit_to_data.objective_IM_III((A,tau,C-alpha*np.exp(alpha*tau)*(epsilon)),t_i) for t_i in list(t_trans)])
    R_trans = R
    # Return the transformed solution
    return t_trans,R_trans
#----------------------------------------------------------------------------------
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

