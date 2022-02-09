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
import numpy as np  # For the exponential function
from scipy.optimize import curve_fit # Curve fit hey :)
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
def objective_PLM(t,A,gamma):
    # Return the output
    return A*(t**gamma)
# ---------------------------------------------------------------------------------------
# Function 2: "objective_IM_II"
# The function returns the objective value of the exponential model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "tau" being the delay parameter.
def objective_IM_II(t,A,tau,alpha):
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
def PE_risk_profiles(t, R, model_str,fit_alternative):
    # We have two models to consider, namely the PLM and the IM-II. In order to do the ODR
    # based model fitting, we need to construct a model object and a start guess for the
    # parameters.
    if model_str == "PLM":# The PLM
        # Guess parameters
        guess = np.array([1, 3])
        # Here, we estimate both A and gamma
        opt_para, _ = curve_fit(lambda t, A, gamma: objective_PLM(t, A, gamma), t, R, guess)
        # Extract parameters
        A, gamma = opt_para
        # Return the simulated output
        R_hat = np.array([objective_PLM(t[i],A,gamma) for i in range(len(t))])        
    elif model_str == "IM-II": # The IM-II
        # Fit alpha or not
        if fit_alternative == "all":
            # Make a start guess
            guess = np.array([5, 50, 0.044])
            # Fit parameters
            opt_para, _ = curve_fit(lambda t, A, tau, alpha: objective_IM_II(t, A, tau, alpha), t, R,guess)
            # Extract optimal parameters
            A, tau, alpha = opt_para
            # Return the simulated output
            R_hat = np.array([objective_IM_II(t[i],A,tau,alpha) for i in range(len(t))])
        elif fit_alternative == "not alpha":
            # Make a start guess
            guess = np.array([5, 50])
            # Fit parameters
            opt_para, _ = curve_fit(lambda t, A, tau: objective_IM_II(t, A, tau, 0.044), t, R,guess)            
            # Extract optimal parameters
            A, tau = opt_para
            # Return the simulated output
            R_hat = np.array([objective_IM_II(t[i],A,tau,0.044) for i in range(len(t))])
    # Calculate the root mean square as well
    RMS = np.sum(np.array([(R[i]-R_hat[i])**2 for i in range(len(R))]))
    # Return the fitted_model
    return opt_para, R_hat, RMS


# ---------------------------------------------------------------------------------------
# Function 4: "objective_PLM_log"
# The function returns the objective value of the power law model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "gamma" being the power in the power law model.
def objective_PLM_log(t,A,gamma):
    # Return the output
    return A+gamma*np.log(t)
# ---------------------------------------------------------------------------------------
# Function 2: "objective_IM_II"
# The function returns the objective value of the exponential model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "tau" being the delay parameter.
def objective_IM_II_log(t,A,tau,alpha):
    # Return the logarithm of the output
    return A - np.log(np.exp(np.exp(-alpha*(t-tau)))-1)
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
def PE_risk_profiles_log(t, R, model_str,fit_alternative):
    # Convert the Concentrations to logarithms
    R_log = np.array([np.log(R[i]) for i in range(len(R))])
    # We have two models to consider, namely the PLM and the IM-II. In order to do the ODR
    # based model fitting, we need to construct a model object and a start guess for the
    # parameters.
    if model_str == "PLM":# The PLM
        # Guess parameters
        guess = np.array([np.log(1e-7), 3])
        # Here, we estimate both A and gamma
        opt_para, _ = curve_fit(lambda t, A, gamma: objective_PLM_log(t, A, gamma), t, R_log, guess)
        # Extract parameters
        A, gamma = opt_para
        # Return the simulated output
        R_hat = np.array([np.exp(objective_PLM_log(t[i],A,gamma)) for i in range(len(t))])        
    elif model_str == "IM-II": # The IM-II
        # Fit alpha or not
        if fit_alternative == "all":
            # Make a start guess
            guess = np.array([np.log(27), 50, 0.044])
            # Fit parameters
            opt_para, _ = curve_fit(lambda t, A, tau, alpha: objective_IM_II_log(t, A, tau, alpha), t, R_log,guess)
            # Extract optimal parameters
            A, tau, alpha = opt_para
            # Return the simulated output
            R_hat = np.array([np.exp(objective_IM_II_log(t[i],A,tau,alpha)) for i in range(len(t))])
        elif fit_alternative == "not alpha":
            # Make a start guess
            guess = np.array([np.log(27), 50])
            # Fit parameters
            opt_para, _ = curve_fit(lambda t, A, tau: objective_IM_II_log(t, A, tau, 0.044), t, R_log,guess)            
            # Extract optimal parameters
            A, tau = opt_para
            # Return the simulated output
            R_hat = np.array([np.exp(objective_IM_II_log(t[i],A,tau,0.044)) for i in range(len(t))])
    # Calculate the root mean square as well
    RMS = np.sum(np.array([(R[i]-R_hat[i])**2 for i in range(len(R))]))
    # Return the fitted_model
    return opt_para, R_hat, RMS
