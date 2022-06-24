# =================================================================================
# =================================================================================
# Script:"write_output"
# Date: 2022-02-15
# Implemented by: Johannes Borgqvist
# Description:
# We present one function which saves the data from the model fitting into a
# csv-file that is stored in the Data-folder. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import csv # In order to write to csv format
import numpy as np # To take the exponential of data
import pandas as pd # Pandas for saving numpy arrays efficiently
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# Function 1: "save_data_PE"
# The function saves the data from the parameter estimation
# or curve fitting procedure. It takes in the following input:
# 1. "data_str": Defining the data set,
# 2. "model_str": Defining the model that was fitted to the data,
# 3. "optimal_fits": Defines the structure with all parameters, uncertainties, least squares etc.,
# 4. "R_squared": Defines the R_squared value of the fit,
# 5. "fit_str": Defines the methodology by which the models were fitted to the data. 
def save_data_PE(data_str,model_str,optimal_fits,RMS,fit_str):
    # Define the file name
    file_name = "../Output/data_" + data_str + "_model_" + model_str + "_fit_type_" + fit_str + ".csv" 
    # Define the headers and the data for the different models
    if model_str == "PLM":
        header = ["Data", "Model", "Fitting method","RMS", "Parameter, A", "Std error, A", "Parameter, gamma", "Std error, gamma", "Sum of squares", "Sum of squares variable", "Sum of squares response variable"]
        data = [data_str, model_str, fit_str, RMS, optimal_fits.beta[0], optimal_fits.sd_beta[0], optimal_fits.beta[1], optimal_fits.sd_beta[1], optimal_fits.sum_square, optimal_fits.sum_square_delta, optimal_fits.sum_square_eps]
    elif model_str == "IM":
        header = ["Data", "Model", "Fitting method","RMS", "Parameter, A", "Std error, A", "Parameter, tau", "Std error, tau","Parameter, C", "Std error, C", "Parameter, alpha", "Std error, alpha", "Sum of squares", "Sum of squares variable", "Sum of squares response variable"]
        data = [data_str, model_str, fit_str, RMS, optimal_fits.beta[0], optimal_fits.sd_beta[0], optimal_fits.beta[1], optimal_fits.sd_beta[1],optimal_fits.beta[2], optimal_fits.sd_beta[2], optimal_fits.beta[3], optimal_fits.sd_beta[3], optimal_fits.sum_square, optimal_fits.sum_square_delta, optimal_fits.sum_square_eps]        
    # Write the data to the defined file    
    with open(file_name, "w", encoding="UTF8") as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
        # write the data
        writer.writerow(data)
# Function 2: "save_data_symmetry_based_model_selection"
# The function saves the data from the symmetry based model selection procedure.
# It takes the following inputs:
# 1. data_str: a string defining which data we are working with,
# 2. model_str: a string defining the model we are working with,
# 3. epsilon: a numpy array of the transformation parameters epsilon,
# 4. RMS: a numpy array of the RMS-values,
# 5. fitted_parameters: a list of the fitted parameters,
# 6. inverse_parameters: a list of the inverse parameters.
def save_data_symmetry_based_model_selection(data_str,model_str,epsilon,RMS,fitted_parameters,inverse_parameters):
    # Define the file name
    file_name = "../Output/data_" + data_str + "_model_" + model_str + "_symmetry_based_model_selection.csv"
    # Define list for all parameters
    if model_str == "PLM":
        # Optimal parameters when model is fitted to transformed data
        A_list = [fitted_parameters[index][0] for index in range(len(fitted_parameters))]
        gamma_list = [fitted_parameters[index][1] for index in range(len(fitted_parameters))]
        # Parameters of the inversely transformed model
        A_inv_list = [inverse_parameters[index][0] for index in range(len(inverse_parameters))]
        gamma_inv_list = [inverse_parameters[index][1] for index in range(len(inverse_parameters))]
        # Create our massive array
        temp_data = np.array([list(epsilon),list(RMS),A_list,gamma_list,A_inv_list,gamma_inv_list])
        # Create a list of the names
        col_names = ["Transformation parameter epsilon", "RMS " + data_str + " " + model_str, "Parameter A", "Parameter gamma", "Inverse parameter A", "Inverse parameter gamma"]
    elif model_str == "IM":
        # Optimal parameters when model is fitted to transformed data
        A_list = [fitted_parameters[index][0] for index in range(len(fitted_parameters))]
        tau_list = [fitted_parameters[index][1] for index in range(len(fitted_parameters))]
        C_list = [fitted_parameters[index][2] for index in range(len(fitted_parameters))]
        alpha_list = [fitted_parameters[index][3] for index in range(len(fitted_parameters))]        
        # Parameters of the inversely transformed model
        A_inv_list = [inverse_parameters[index][0] for index in range(len(inverse_parameters))]
        tau_inv_list = [inverse_parameters[index][1] for index in range(len(inverse_parameters))]
        C_inv_list = [inverse_parameters[index][2] for index in range(len(inverse_parameters))]
        alpha_inv_list = [inverse_parameters[index][3] for index in range(len(inverse_parameters))]
        # Create our massive array
        temp_data = np.array([list(epsilon),list(RMS),A_list,tau_list,C_list,alpha_list,A_inv_list,tau_inv_list,C_inv_list,alpha_inv_list])
        # Create a list of the names
        col_names = ["Transformation parameter epsilon", "RMS " + data_str + " " + model_str, "Parameter A", "Parameter tau", "Parameter C", "Parameter alpha", "Inverse parameter A", "Inverse parameter gamma", "Inverse parameter C", "Inverse parameter alpha"]
    # Transpose the data so we get a fewer columns than rows
    data = temp_data.T
    # Create a pandas data frame
    df = pd.DataFrame(data,columns=col_names)
    # Finally, save the data frame to the given file name
    df.to_csv(file_name)

