# =================================================================================
# =================================================================================
# Script:"symmetry_based_model_selection"
# Date: 2022-02-15
# Implemented by: Johannes Borgqvist
# Description:
# The script involves all the symmetry functionalities used for transforming the
# functions as well as conducting the symmetry based model selection. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import read_data  # Home-made
import write_output  # Home-made
import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import fmin_cobyla # To find the orthogonal distance
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
    t_trans = np.asarray([PLM_symmetry(t[index],R[index],epsilon)[0] for index in range(len(t))])
    R_trans = R
    # Return the transformed solution
    return t_trans,R_trans
#----------------------------------------------------------------------------------
# FUNCTION 4: PLM_objective
# This function returns the objective value of the PLM which is used to calculate
# the orthogonal distance between the data point (t_data,R_data) and the solution
# the solution curve (t,R(t)) of the power law model.
def PLM_objective(t,args):
    # Extract the parameters
    A = args[0]
    gamma = args[1]
    # Return the objective value of the PLM
    return A*(t**gamma)
#----------------------------------------------------------------------------------
# FUNCTION 5: PLM_constraint
# This function is a help function that is necessary for finding the orthogonal distance
# between a data point and a point on a solution curve of the PLM.
def PLM_constraint(Model_Point,*args):
    # Extract the data point
    t,R = Model_Point
    # Pass the constraint to the optimisation
    return PLM_objective(t,args) - R
#----------------------------------------------------------------------------------
# FUNCTION 6: IM_III_symmetry
# The symmetry of the immunological model IM-III
def IM_III_symmetry(t,R,epsilon,tau,alpha):
    # Define t_hat recursively
    t_hat =  np.log(np.log(np.exp(np.exp(-alpha*(t-tau))) - (alpha*np.exp(alpha*tau)*epsilon) ))
    t_hat = tau - ((t_hat)/(alpha))
    R_hat = R
    return t_hat,R_hat
#----------------------------------------------------------------------------------
# FUNCTION 7: IM_III_transformation
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
    t_trans = np.asarray([IM_III_symmetry(t[index],R[index],epsilon,tau,alpha)[0] for index in range(len(t))])
    R_trans = R
    # Return the transformed solution
    return t_trans,R_trans
#----------------------------------------------------------------------------------
# FUNCTION 9: IM_III_objective
# This function returns the objective value of the IM_III which is used to calculate
# the orthogonal distance between the data point (t_data,R_data) and the solution
# the solution curve (t,R(t)) of the immunological model.
def IM_III_objective(t,args):
    # Extract the parameters
    A = args[0]
    tau = args[1]
    C = args[2]
    # Define the value of the fixed parameter alpha
    alpha = 0.044
    # Return the objective value of the IM-III
    return ((A)/(np.exp(np.exp(-alpha*(t-tau)))-C))
    #return np.log(A)-np.log(np.exp(np.exp(-alpha*(t-tau)))-C)
#----------------------------------------------------------------------------------
# FUNCTION 10: IM_III_constraint
# This function is a help function that is necessary for finding the orthogonal distance
# between a data point and a point on a solution curve of the IM-III.
def IM_III_constraint(Model_Point,*args):
    # Extract the data point
    t,R = Model_Point
    # Pass the constraint to the optimisation
    return IM_III_objective(t,args) - R
    #return IM_III_objective(t,args) - np.log(R)
#----------------------------------------------------------------------------------
# FUNCTION 11: SS_res_model_data
# This function calculates the Sum of Squares (SS) of the residuals between the
# data and the model
def SS_res_model_data(Model_Point,Data_Point):
    # Extract the model points
    t_model,R_model = Model_Point
    # Extract the model points
    t_data,R_data = Data_Point
    # Return the objective value
    return ((t_model - t_data)**2 + (R_model - R_data)**2)
#----------------------------------------------------------------------------------
# FUNCTION 12: symmetry_based_model_selection
# This function performs the symmetry based model selection. It does so by
# taking the following input:
# 1. DATA: t_data, the ages of the patients,
# 2. DATA: R_data, the number of incidences of cancer,
# 3. epsilon_vector, the transformation parameters we transform the data with,
# 4. model_str, a string indicating which model we consider (i.e. the PLM or the IM-III)
def symmetry_based_model_selection(t_data,R_data,epsilon_vector,parameters,model_str):
    # SETTING UP THE MODEL SELECTION
    # Allocate memory for the output which is the RMS(epsilon)
    # indicating how the fit to transformed data changes as a
    # function of the transformation parameter epsilon
    RMS_transf = []
    # Extract the parameters depending on the model
    if model_str == "PLM":
        # Extract the parameter A
        A = parameters[0]
        # Extract the parameter gamma
        gamma = parameters[1]
    elif model_str == "IM-III":
        # Extract the parameter A
        A = parameters[0]
        # Extract the parameter tau
        tau = parameters[1]
        # Extract the parameter C
        C = parameters[2]
        # Fix the parameter alpha
        alpha = 0.044
    # Loop over the epsilon vectors, transform the data and
    # calculate the fit of the candidate models to the transformed data
    for epsilon in list(epsilon_vector):
        # Transform the data
        if model_str == "PLM":
            t_trans = np.array([PLM_symmetry(t_data[index],R_data[index],epsilon)[0] for index in range(len(t_data))])
            R_trans = np.array([PLM_symmetry(t_data[index],R_data[index],epsilon)[1] for index in range(len(t_data))])
        elif model_str == "IM-III":
            # Transform the data
            t_trans = np.array([IM_III_symmetry(t_data[index],R_data[index],epsilon,tau,alpha)[0] for index in range(len(t_data))])
            R_trans = np.array([IM_III_symmetry(t_data[index],R_data[index],epsilon,tau,alpha)[1] for index in range(len(t_data))])
        # Allocate memory for the sum of squares (SS)
        SS = 0            
        # Loop over the transformed time series
        for time_series_index in range(len(t_data)):
            # Extract a data point
            Data_point = (t_trans[time_series_index],R_trans[time_series_index])
            # Update the curve specific parameters that are transformed corresponding to
            # the parameter A in the case of the PLM and the parameter C in the case of
            # the IM-III. Also, we find the orthogonal point on the solution curve using
            # fmin_cobyla
            if model_str == "PLM":
                # Update the parameter A
                A_transf = A*np.exp(-gamma*epsilon)
                #A_transf = A*np.exp(-gamma*epsilon)                
                # Find the orthogonal point on the solution curve (t,R(t)) of the PLM
                Model_point = fmin_cobyla(SS_res_model_data, x0=list(Data_point), cons=[PLM_constraint], args=(Data_point,),consargs=(A_transf,gamma))          
            elif model_str == "IM-III":
                # Update the parameter C of the IM-III
                C_transf = C - alpha*np.exp(alpha*tau)*epsilon
                #C_transf = C + alpha*np.exp(alpha*tau)*epsilon                
                # Find the orthogonal point on the solution curve (t.R(t)) of the IM-III
                Model_point = fmin_cobyla(SS_res_model_data, x0=list(Data_point), cons=[IM_III_constraint], args=(Data_point,),consargs=(A,tau,C_transf))
            # Add the squared distances to our growing sum of squares (SS)
            SS += SS_res_model_data(Model_point,Data_point)
        # Lastly, append the root mean squared calculated based on the SS-value
        RMS_transf.append(np.sqrt(SS/len(t_data)))
    # Lastly, cast the RMS-values as an array before we return
    return np.array(RMS_transf)    
