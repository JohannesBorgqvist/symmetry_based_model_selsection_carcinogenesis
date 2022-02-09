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
    #alpha = 0.044
    # Extract the parameters to be fitted
    A, tau, C, alpha = parameters
    #A, tau, C = parameters
    # Return the logarithm of the output
    #return ((A)/(np.exp(np.exp(-alpha*(t-tau)))-1))
    return ((A)/(np.exp(np.exp(-alpha*(t-tau)))-C*A))
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
def PE_risk_profiles(t, R, model_str,fit_string,fixed_parameters):
    # Define the data at hand as a structure on which
    # the scipy odr package will do the estimation
    #data = RealData(t, R, 0.03, 0.01)
    #data = RealData(t, R)
    #data = Data(t, R, 0.01, 0.01)
    data = Data(t, R)
    # Define the number of start guesses we will try. Note that the optimisation for finding
    # the optimal parameters is local, and the cost function is non-convex meaning that there
    # are multiple local optima which can be found for different start guesses. Therefore, 
    # we use a multiple shooting technique were we test multiple start guesses and then we
    # save the optimal parameters that resulted in the best fit.
    num_of_start_guesses = 7
    # We have two models to consider, namely the PLM and the IM-II. In order to do the ODR
    # based model fitting, we need to construct a model object and a start guess for the
    # parameters.
    if model_str == "PLM":# The PLM
        # Define the model
        model = Model(objective_PLM)
        # Define a set of parameter values for the multishooting technique
        # where we do many local optimisations starting from these startguesses.
        # In the end, we pick the parameters resulting in the lowest minima. This
        # is because we have a non-convex optimisation problem with many local minima.
        A_vec = np.linspace(0.001,10,num_of_start_guesses,endpoint=True)
        gamma_vec = np.linspace(1,10,num_of_start_guesses,endpoint=True)        
        # Extract the first start guess
        parameter_guesses = [[A, gamma] for A in A_vec for gamma in gamma_vec]
        #parameter_guess = [1 , 5]
        #parameter_guess = [0.02, 3]
    elif model_str == "IM-II": # The IM-II
        # Define the model
        model = Model(objective_IM_II)
        # Define the start guess [A, tau, C, alpha]
        #parameter_guess = [1.5, 20, 0.1, 0.044] # Excellent start guess
        # Define a set of parameter values for the multishooting technique
        # where we do many local optimisations starting from these startguesses.
        # In the end, we pick the parameters resulting in the lowest minima. This
        # is because we have a non-convex optimisation problem with many local minima.
        A_vec = np.linspace(0.01,50,num_of_start_guesses,endpoint=True)
        C_vec = np.linspace(0.01,10,num_of_start_guesses,endpoint=True)
        tau_vec = np.linspace(1,200,num_of_start_guesses,endpoint=True)
        alpha_vec = np.linspace(0.001,0.1,num_of_start_guesses,endpoint=True)                
        # Extract the first start guess
        parameter_guesses = [[A, tau, C, alpha] for A in A_vec for tau in tau_vec for C in C_vec for alpha in alpha_vec]        
    # Set an initial value of the RMS so that it will be updated    
    RMS = 50000
    # Also, initiate the other two outputs that we return
    fitted_model = 0
    R_hat = np.array([1, 2, 3])
    #------------------------------------------------------------------------------------------
    # Now we re-do this analysis for the remaining start guesses of the parameters
    for parameter_guess in parameter_guesses:
        # Define the model fitting object
        if len(fixed_parameters)==0:
            # Set up ODR with the model and data.
            odr = ODR(data, model, parameter_guess)        
        else:
            # Set up ODR with the model and data and fixed parameters.
            odr = ODR(data, model, parameter_guess,ifixb=fixed_parameters)        
        # Define whether we should use standard least square or the fancy ODR fitting
        if fit_string == "LS": # Least square fitting
            odr.set_job(fit_type=2)
        elif fit_string == "ODR": # Orthogonal distance regression
            odr.set_job(fit_type=0)
        # Run the regression.
        fitted_model_temp = odr.run()
        # Generate the fitted time series as well in the two cases
        if model_str == "PLM":# The PLM    
            R_hat_temp = np.array([objective_PLM(fitted_model_temp.beta,t_i) for t_i in list(t)])
        elif model_str == "IM-II": # The IM-II
            R_hat_temp = np.array([objective_IM_II(fitted_model_temp.beta,t_i) for t_i in list(t)])
        # Calculate the RMS
        RMS_temp = np.sqrt(((fitted_model_temp.sum_square)/(len(R))))
        # Lastly, if we obtained a better fit than the current minimal fit, we save that fit instead.
        if RMS_temp < RMS:
            RMS = RMS_temp
            R_hat = R_hat_temp
            fitted_model = fitted_model_temp
    #------------------------------------------------------------------------------------------        
    # Return the fitted_model
    return fitted_model, R_hat, RMS
