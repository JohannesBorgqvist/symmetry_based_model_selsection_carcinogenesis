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
    elif model_str == "IM-III":
        header = ["Data", "Model", "Fitting method","RMS", "Parameter, A", "Std error, A", "Parameter, tau", "Std error, tau","Parameter, C", "Std error, C", "Sum of squares", "Sum of squares variable", "Sum of squares response variable"]
        data = [data_str, model_str, fit_str, RMS, optimal_fits.beta[0], optimal_fits.sd_beta[0], optimal_fits.beta[1], optimal_fits.sd_beta[1],optimal_fits.beta[2], optimal_fits.sd_beta[2], optimal_fits.sum_square, optimal_fits.sum_square_delta, optimal_fits.sum_square_eps]        
    # Write the data to the defined file    
    with open(file_name, "w", encoding="UTF8") as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
        # write the data
        writer.writerow(data)
