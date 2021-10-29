# =================================================================================
# =================================================================================
# Script:"write_output"
# Date: 2021-09-20
# Implemented by: Johannes Borgqvist
# Description:
# The function write input data to output files which are mostly tex-files that
# are plotted in LaTeX. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import csv # In order to write to csv format
import numpy as np
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# Function 1: "plot_LaTeX_2D"
# The function takes the following input:
# 1. "t" being the list of t-values or input values,
# 2. "y" being the list of y-values or output values,
# 3. "file_str" being the string defining where the output
# file with format tex is stored,
# 4. "plot_str" is the string defining properties such as
# colour and linewidth,
# 5. "legend_str" is a string containing the legend of the plot.
# The function then creates a tex file which can be plotted
# in LaTeX using pgfplots. 
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for i in range(len(t)):
        temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
#---------------------------------------------------------------------------------------    
# Function 2: "save_data_PE"
# The function saves the data from the parameter estimation
# or curve fitting procedure. It takes in the following input:
# 1. "data_str": Defining the data set,
# 2. "model_str": Defining the model that was fitted to the data,
# 3. "opt_para": Defining the optimal parameters,
# 4. "SS": Defining the sum_of_squares
def save_data_PE(data_str,model_str,opt_para,R_adj,alpha):
    # Define the file name
    file_name = "../Output/data_" + data_str + "_model_" + model_str + ".csv" 
    # Define the headers and the data for the different models
    if model_str == "exponential":
        header = ["Data", "Model", "Adjusted R^2, R_{adj}^2", "Parameter, log(A)", "Parameter alpha"]
        data = [data_str, model_str, R_adj, np.log(opt_para[0]), alpha]
    elif model_str == "power_law":
        header = ["Data", "Model", "Adjusted R^2, R_{adj}^2", "Parameter, log(A)", "Parameter, gamma"]
        data = [data_str, model_str, R_adj, np.log(opt_para[0]), opt_para[1]]
    elif model_str == "mixed":
        header = ["Data", "Model", "Adjusted R^2, R_{adj}^2", "Parameter, log(A)", "Parameter, tau", "Parameter, alpha"]
        data = [data_str, model_str, R_adj, np.log(opt_para[0]), opt_para[1], alpha]        
    # Write the data to the defined file    
    with open(file_name, "w", encoding="UTF8") as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
        # write the data
        writer.writerow(data)
