# TiO2_inhalation
This repository contains the data, fitting and visualisation R codes developed under Deliverable 6.3 of the NanoCommons project. The data of the xlsx files were obtained from https://www.doi.org/10.1186/s12989-019-0303-7 

The biodist_data.xlsx contains the mean values in ng in sheet 1, the standrard deviations in sheet 2 and the input parameters in sheet 3. The urine.xlsx and feces.xlsx contain the correpsonging excreta data. Results.RData contains the results of the fitting process under the 'fit' variable. The fitting_R.R and fitting_stan.stan files contain the fitting R and stan code, respectively. Finally, the group_plots.R and VPC_plots.R include the code for generating the group plots and VPC plots. 

To use the code, the user is asked to put all files in a common file and then in the first lines insert the path of the file in a string format.
