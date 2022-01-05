# Symmetry based model selection in the context of carcinogenesis
*Written by:* Johannes Borgqvist<br>
*Date:* 2021-12-15<br>




## Installation of the required packages using anaconda
The easiest way to install all relevant libraries and packages is to use [anaconda](https://docs.anaconda.com/anaconda/install/index.html). When anaconda has been installed properly, create an environment from the provided yml file as follows:<br>

*conda env create -f symmetry\_based\_model\_selection.yml*<br>

After this, the environment is activated with the command:<br>
*conda activate symmetry\_based\_model\_selection*,<br>
or if this does not work you can activate the environment using the command 
*source activate symmetry\_based\_model\_selection*,<br>
and this environment is deactivated using<br>
*conda deactivate*.<br> 

Alternatively, you can create the anaconda environment and install all relevant libraries individually. This approach seems to not work especially when it comes to installing R using anaconda. Anyhow, here follows some instruction for creating the same environment (again this approach might not work on some systems).<br> 


Create an anaconda environment with python, numpy, scipy, matplotlib and pandas:<br>
*conda create -n symmetry_based_model_selection python numpy matplotlib pandas scipy*<br>
Add the R essentials bundle:<br>
*conda install -n symmetry_based_model_selection -c r r-essentials*<br>
Add nonnest-2 in R<br>
*conda install -n symmetry_based_model_selection -c r r-nonnest2*<br>
Add drc in R<br>
*conda install -n symmetry_based_model_selection -c conda-forge r-drc*<br>
To activate your new environment, use the command<br>
*conda activate symmetry_based_model_selection*<br>
or if this does not work use<br>
*source activate symmetry_based_model_selection*<br>
and to deactivate the environment type
*conda deactivate*<br>
I had some issues with installing R using anaconda on my personal laptop, but it worked on my work computer. A piece of advice that I saw somewhere was that if there are any issues with conflicting 
versions of various packages, then it is always good to update anaconda. This is done in the following way:<br>
*conda update --all*<br>


## Running the scripts
In an ubuntu terminal, the python scripts are executed as follows:<br>
*python [scriptname].py*<br>
and the R-scripts are executed as follows:<br>
*Rscript [scriptname].R*<br>



