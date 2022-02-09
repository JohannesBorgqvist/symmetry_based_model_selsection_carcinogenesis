# =================================================================================
# =================================================================================
# Script:"fit_to_data"
# Date: 2021-09-23
# Implemented by: Johannes Borgqvist
# Description:
# The program conducts a simple curve fitting to a time series and returns the
# fitted model.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
from scipy.odr import * # For calculating the total sum of squares
import numpy as np  # For the exponential function
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# ---------------------------------------------------------------------------------------
# Function 1: "objective_PLM"
# The function returns the objective value of the power law model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "gamma" being the power in the power law model.
def objective_PLM(parameters,t):
    # Extract the parameters to be fitted
    A, gamma = parameters
    # Return the output
    return A*(t**gamma)
# ---------------------------------------------------------------------------------------
# Function 2: "objective_IM_II"
# The function returns the objective value of the exponential model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "tau" being the delay parameter.
def objective_IM_II(parameters, t):
    # Define the constant alpha
    alpha = 0.044
    # Extract the parameters to be fitted
    A, tau = parameters
    # Return the logarithm of the output
    return ((A)/(np.exp(np.exp(-alpha*(t-tau)))-1))
# ---------------------------------------------------------------------------------------
# Function 3: "PE_risk_profiles"
# The function fits one of the two candidate models (the PLM or the IM-II) to experimental
# data. It takes four inputs:
# 1. The vector t being the ages of the patients,
# 2. The vector R being the number of incidences of cancer for a given age,
# 3. The string model_str defining whether it is the PLM or the IM-II that is fitted to the data,
# 4. The string fit_type defining whether it is the standard least square or the Orthogonal Distance Regression
# (ODR) measure that is used.
# The function returns three outputs:
# 1. The structure fitted_model containing the fits, the optimal parameters, the parameter variances etc.
# 2. The vector R_hat containing the simulated incidences with the optimal parameters where each age is given by the data vector t,
# 3. The float number R_squared which quantifies the R_squared fit of the candidate model to the data at hand. 
def PE_risk_profiles(t, R, model_str,fit_string):
    # Define the data at hand as a structure on which
    # the scipy odr package will do the estimation
    #data = RealData(t, R, 0.03, 0.01)
    data = RealData(t, R)
    # We have two models to consider, namely the PLM and the IM-II. In order to do the ODR
    # based model fitting, we need to construct a model object and a start guess for the
    # parameters.
    if model_str == "PLM":# The PLM
        # Define the model
        model = Model(objective_PLM)
        # Define the start guess [A, gamma]
        parameter_guess = [1, 3]
    elif model_str == "IM-II": # The IM-II
        # Define the model
        model = Model(objective_IM_II)
        # Define the start guess [A, tau]
        parameter_guess = [30, 30]
    # Set up ODR with the model and data.
    odr = ODR(data, model, parameter_guess)
    # Define whether we should use standard least square or the fancy ODR fitting
    if fit_string == "LS": # Least square fitting
        odr.set_job(fit_type=2)
    elif fit_string == "ODR": # Orthogonal distance regression
        odr.set_job(fit_type=0)        
    # Run the regression.
    fitted_model = odr.run()
    # Generate the fitted time series as well in the two cases
    if model_str == "PLM":# The PLM    
        R_hat = np.array([objective_PLM(fitted_model.beta,t_i) for t_i in list(t)])
    elif model_str == "IM-II": # The IM-II
        R_hat = np.array([objective_IM_II(fitted_model.beta,t_i) for t_i in list(t)])
    # Calculate the sum of squares
    #SS = np.sum(np.array([(R[index]-R_hat[index])**2 for index in range(len(R))]))
    # Calculate the RMS
    RMS = np.sqrt(((fitted_model.sum_square)/(len(R))))
    #RMS = np.sqrt(((SS)/(len(R))))
    # Return the fitted_model
    return fitted_model, R_hat, RMS
