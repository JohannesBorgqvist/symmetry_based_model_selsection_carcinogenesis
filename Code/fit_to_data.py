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
from scipy.optimize import curve_fit  # For the curve fitting
import numpy as np  # For the exponential function
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# ---------------------------------------------------------------------------------------
# Function 1: "objective_exponential"
# The function returns the objective value of the exponential model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,


def objective_exponential(t, A, alpha):
    # Return the output
    # return A*np.exp(alpha*t)
    # return np.log(A) + (alpha*t)
    return A + (alpha*t)
# ---------------------------------------------------------------------------------------
# Function 2: "objective_power_law"
# The function returns the objective value of the power law model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "gamma" being the power in the power law model.
def objective_power_law(t, A, gamma):
    # Return the output
    # return A*(t**gamma)
    # return np.log(A)+gamma*np.log(t)
    return A+gamma*np.log(t)
# ---------------------------------------------------------------------------------------
# Function 3: "objective_mixed"
# The function returns the objective value of the exponential model and it takes
# the following three inputs
# 1."t" being the ages in the time series (or the independent variable if you will),
# 2. "A" being the scaling parameter in all models,
# 3. "tau" being the delay parameter.


def objective_mixed(t, A, tau, alpha):
    # Return the logarithm of the output
    # return ((A)/(np.exp(np.exp(-alpha*(t-tau)))-1))
    # return np.log(A) - np.log(np.exp(np.exp(-alpha*(t-tau)))-1)
    return A - np.log(np.exp(np.exp(-alpha*(t-tau)))-1)
# ---------------------------------------------------------------------------------------
# Function 4: "calculate_R_adj"
# The function calculates the adjusted R square value.
# It takes the following input:
# 1. R which corresponds to the data,
# 2. R_hat corresponding to the model,
# 3. model_str defining how many parameters that is involved.
# It returns the adjusted R square value


def calculate_R_adj(R, R_hat, model_str):
    # Lastly, we calculate the sum of squares as well.
    # Define the residuals:
    residuals = [np.exp(R[i])-np.exp(R_hat[i]) for i in range(len(R))]
    # Calculate the mean value of the data
    R_exp = np.array([np.exp(R[i]) for i in range(len(R))])
    mean_R = np.mean(R_exp)
    # Calculate the total variation
    total_variation = [np.exp(R[i])-mean_R for i in range(len(R))]
    # Allocate memory for the sum of squares
    res_sum = 0
    tot_sum = 0
    # Loop over the residuals and add them to the sum of squares
    for res_index in range(len(residuals)):
        res_sum += residuals[res_index]**2
        tot_sum += (total_variation[res_index])**2
    # Take 1 minus this value to get the R^2
    R_squared = 1 - ((res_sum)/(tot_sum))
    #correlation_matrix = np.corrcoef(x_values, y_values)
    #correlation_xy = correlation_matrix[0,1]
    #r_squared = correlation_xy**2
    # Calculate the adjusted R^2
    if model_str == "mixed":
        k = 2  # Number of parameters
    else:
        k = 2  # Number of parameters
    # Number of samples
    n = len(R)
    # Calculate the R adjusted
    num = (1-R_squared)*(n-1)
    denom = n-k-1
    R_adj = 1-((num)/(denom))
    # Return the output
    return R_adj
# Function 5: "PE_risk_profiles"
# The function takes in five inputs:
# 1."t" being the years in the time series (the variable),
# 2. "R" being the risk incidence in the time series (the output),
# 3. "model_str" the string defining the model we work with and there are
# three options: "exponential", "power_law", "mixed",
# 4. "lb" being a list of the lower bounds,
# 5. "ub" being a list of the upper bounds.
# The function returns three outputs:
# "R_hat" being the output of the model at hand,
# "opt_para" being the estimated parameters,
# "SS" being the sum of squares.


def PE_risk_profiles(t, R, model_str, parameters):
    # Define the model parameter for transmission
    alpha = parameters[0]
    # Find the optimal parameters to the model at hand
    if model_str == "exponential":  # Exponential model
        # Fit the exponential model to the data at hand
        if len(parameters) == 1:
            opt_para, _ = curve_fit(
                lambda t, A: objective_exponential(t, A, alpha), t, R)
            # Extract the optimal parameter
            A_opt = opt_para[0]
        # Allocate an numpy array for the simulated output
        # R_hat
        R_hat = np.zeros(t.shape)
        # Loop through this array and fill its value
        for index in range(len(R_hat)):
            R_hat[index] = objective_exponential(t[index], A_opt, alpha)
    elif model_str == "power_law":  # Power law
        # Fit the power law model to the data at hand
        if len(parameters) == 1:
            # Here, we estimate both A and gamma
            opt_para, _ = curve_fit(
                lambda t, A, gamma: objective_power_law(t, A, gamma), t, R)
            # Extract the optimal parameter
            A_opt, gamma_opt = opt_para
        else:
            # Here, we estimate only A
            gamma_opt = parameters[1]
            # Fit to the data
            opt_para, _ = curve_fit(
                lambda t, A: objective_power_law(t, A, gamma_opt), t, R)
            # Calculate the output of the model
            A_opt = opt_para
        # Allocate an numpy array for the simulated output
        # R_hat
        R_hat = np.zeros(t.shape)
        # Loop through this array and fill its value
        for index in range(len(R_hat)):
            R_hat[index] = objective_power_law(t[index], A_opt, gamma_opt)
    elif model_str == "mixed":  # Mixed model
        if len(parameters) == 1:
            # Fit the IM-II to the data at hand
            opt_para, _ = curve_fit(
                lambda t, A, tau: objective_mixed(t, A, tau, alpha), t, R)
            # Extract the optimal parameter
            A_opt, tau_opt = opt_para
            # Note that in the case above we estimate both A and tau
        else:
            # Here, we estimate only A
            tau_opt = parameters[1]
            # Fit the IM-II to the data at hand
            opt_para, _ = curve_fit(
                lambda t, A: objective_mixed(t, A, tau_opt, alpha), t, R)
            # Extract the optimal parameter
            A_opt = opt_para
        # Allocate an numpy array for the simulated output
        # R_hat
        R_hat = np.zeros(t.shape)
        # Loop through this array and fill its value
        for index in range(len(R_hat)):
            R_hat[index] = objective_mixed(t[index], A_opt, tau_opt, alpha)
    else:  # We have no other candidates
        return [], (0, 0), 0
    # Calculate the adjusted R square value
    R_adj = calculate_R_adj(R, R_hat, model_str)
    # Return the output
    return R_hat, opt_para, R_adj
# ---------------------------------------------------------------------------------------
