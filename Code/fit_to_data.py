# =================================================================================
# =================================================================================
# Script:"fit_to_data"
# Date: 2022-02-24
# Implemented by: Johannes Borgqvist
# Description:
# The program conducts a simple curve fitting to a time series and returns the
# fitted model.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import symmetry_toolbox # Home-made
import fit_to_data  # Home-made
from scipy.odr import * # For calculating the total sum of squares
from scipy.optimize import fmin_cobyla # To find the orthogonal distance
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
    #return A + gamma*np.log(t)
# ---------------------------------------------------------------------------------------
# Function 2: "objective_IM_III"
# The function returns the objective value of the exponential model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "tau" being the delay parameter.
def objective_IM_III(parameters, t):
    # Extract the parameters to be fitted
    A, tau, C, alpha = parameters
    # Return the logarithm of the output
    return ((A)/(np.exp(np.exp(-alpha*(t-tau)))-C))
    #return A - np.log(np.exp(np.exp(-alpha*(t-tau)))-C)
# ---------------------------------------------------------------------------------------
# Function 3: "PE_risk_profiles"
# The function fits one of the two candidate models (the PLM or the IM-III) to experimental
# data. It takes four inputs:
# 1. The vector t being the ages of the patients,
# 2. The vector R being the number of incidences of cancer for a given age,
# 3. The string model_str defining whether it is the PLM or the IM-III that is fitted to the data,
# 4. The string fit_type defining whether it is the standard least square or the Orthogonal Distance Regression
# (ODR) measure that is used.
# The function returns three outputs:
# 1. The structure fitted_model containing the fits, the optimal parameters, the parameter variances etc.
# 2. The vector R_hat containing the simulated incidences with the optimal parameters where each age is given by the data vector t,
# 3. The float number R_squared which quantifies the R_squared fit of the candidate model to the data at hand. 
def PE_risk_profiles(t, R, model_str,fit_string,fixed_parameters):
    # Take the logarithm of the data
    #R_log = np.array([np.log(R_temp) for R_temp in list(R)])
    # Define the data at hand as a structure on which
    # the scipy odr package will do the estimation
    #if model_str == "PLM":
    data = Data(t,R)
    #elif model_str == "IM-III":
        #data = Data(t, R_log)
    # Define the number of start guesses we will try. Note that the optimisation for finding
    # the optimal parameters is local, and the cost function is non-convex meaning that there
    # are multiple local optima which can be found for different start guesses. Therefore, 
    # we use a multiple shooting technique were we test multiple start guesses and then we
    # save the optimal parameters that resulted in the best fit.
    num_of_start_guesses = 10
    #num_of_start_guesses = 5
    # We have two models to consider, namely the PLM and the IM-III. In order to do the ODR
    # based model fitting, we need to construct a model object and a start guess for the
    # parameters.
    if model_str == "PLM":# The PLM
        # Define the model
        model = Model(objective_PLM)
        # Define a set of parameter values for the multishooting technique
        # where we do many local optimisations starting from these startguesses.
        # In the end, we pick the parameters resulting in the lowest minima. This
        # is because we have a non-convex optimisation problem with many local minima.
        A_vec = np.linspace(0.00001,1,num_of_start_guesses,endpoint=True)
        gamma_vec = np.linspace(1,10,num_of_start_guesses,endpoint=True)        
        #A_vec = np.linspace(0.0001,0.1,num_of_start_guesses,endpoint=True)
        #gamma_vec = np.linspace(3,5,num_of_start_guesses,endpoint=True)        
        # Extract the first start guess
        parameter_guesses = [[A, gamma] for A in A_vec for gamma in gamma_vec]
    elif model_str == "IM-III": # The IM-III
        # Define the model
        model = Model(objective_IM_III)
        # Define the start guess [A, tau, C, ]
        # Define a set of parameter values for the multishooting technique
        # where we do many local optimisations starting from these startguesses.
        # In the end, we pick the parameters resulting in the lowest minima. This
        # is because we have a non-convex optimisation problem with many local minima.
        #A_vec = np.linspace(1,10,num_of_start_guesses,endpoint=True)
        #C_vec = np.linspace(-2,2,num_of_start_guesses,endpoint=True)
        #tau_vec = np.linspace(1,50,num_of_start_guesses,endpoint=True)
        A_vec = np.linspace(0.01,50,num_of_start_guesses,endpoint=True)
        C_vec = np.linspace(0.01,10,num_of_start_guesses,endpoint=True)
        tau_vec = np.linspace(1,200,num_of_start_guesses,endpoint=True)        
        #alpha_vec = np.linspace(0.0001, 0.01,num_of_start_guesses,endpoint=True)
        alpha_vec = np.array([0.040, 0.045])
        # Extract the first start guess
        parameter_guesses = [[A, tau, C, alpha] for A in A_vec for tau in tau_vec for C in C_vec for alpha in alpha_vec]        
    # Set an initial value of the RMS so that it will be updated    
    RMS = 50000
    # Also, initiate the other two outputs that we return
    fitted_model = 0
    R_hat = np.array([1, 2, 3])
    #------------------------------------------------------------------------------------------
    # Now we do this analysis for the remaining start guesses of the parameters
    for parameter_guess in parameter_guesses:
        # Define the model fitting object
        if len(fixed_parameters)==0:
            # Set up ODR with the model and data.
            odr = ODR(data, model, parameter_guess,sstol=1e-15,partol=1e-15)        
        else:
            # Set up ODR with the model and data and fixed parameters.
            odr = ODR(data, model, parameter_guess,ifixb=fixed_parameters,sstol=1e-15,partol=1e-15)        
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
        elif model_str == "IM-III": # The IM-III
            R_hat_temp = np.array([objective_IM_III(fitted_model_temp.beta,t_i) for t_i in list(t)])
        # Calculate the RMS
        #RMS_temp = symmetry_toolbox.symmetry_based_model_selection(t,R,np.array([0]),fitted_model_temp.beta,model_str)[0]
        # Allocate memory for the sum of squares (SS)
        SS = 0            
        # Loop over the transformed time series
        for time_series_index in range(len(t)):
            # Extract a data point
            Data_point = (t[time_series_index],R[time_series_index])
            # Update the curve specific parameters that are transformed corresponding to
            # the parameter A in the case of the PLM and the parameter C in the case of
            # the IM-III. Also, we find the orthogonal point on the solution curve using
            # fmin_cobyla
            if model_str == "PLM":
                # Find the orthogonal point on the solution curve (t,R(t)) of the PLM
                Model_point = fmin_cobyla(symmetry_toolbox.SS_res_model_data, x0=list(Data_point), cons=[symmetry_toolbox.PLM_constraint], args=(Data_point,),consargs=(fitted_model_temp.beta[0],fitted_model_temp.beta[1]))          
            elif model_str == "IM-III":    
                # Find the orthogonal point on the solution curve (t,R(t)) of the IM-III
                Model_point = fmin_cobyla(symmetry_toolbox.SS_res_model_data, x0=list(Data_point), cons=[symmetry_toolbox.IM_III_constraint], args=(Data_point,),consargs=(fitted_model_temp.beta[0],fitted_model_temp.beta[1],fitted_model_temp.beta[2],fitted_model_temp.beta[3]))
            # Add the squared distances to our growing sum of squares (SS)
            SS += symmetry_toolbox.SS_res_model_data(Model_point,Data_point)
        # Lastly, append the root mean squared calculated based on the SS-value
        RMS_temp = np.sqrt(SS/len(t))        
        # Lastly, if we obtained a better fit than the current minimal fit, we save that fit instead.
        if RMS_temp < RMS:
            RMS = RMS_temp
            R_hat = R_hat_temp
            fitted_model = fitted_model_temp
    #------------------------------------------------------------------------------------------       
    # Return the fitted_model
    return fitted_model, R_hat, RMS
