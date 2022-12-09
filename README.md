Explanation of scripts for:  

Probabilistic inversion of circular phase spectra: Application to two-station phase-velocity dispersion estimation in western Canada 

See publication in GJI for details
This software is provided ‘as-is’. 

Authors: Jeremy M. Gosselin, Pascal Audet, Clement Esteve, and Andrew Schaeffer



 Disclaimer: 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.




 Basic Instructions: 

You will require a list of seismic stations to analyze for the two-station method. There are several python scripts that are used for seismic data downloading and preparation. These are structured in a way to be optimal for running on a large set of seismic station pairs, as there is a lot of internal book-keeping required. 

Data preparation:  

- ‘stations_list_test.txt’  
This file contains all station information. See this file as an example for formatting and creating a station list.  

- ‘teleseismic_eq_download_iris.py’ 
As the name suggests, this script will download all seismic waveform data from IRIS.  Within this script, there are sever user defined variables. See comments in this script for details. 
Waveform data will be stored in the ‘Data’ directory.
Earthquake catalog for the entire project will be stored in ‘CAT/master.pkl’. Earthquake catalogs for each station pair will be stored in the ‘CAT’ directory.    

- ‘teleseismic_eq_noise_spectra.py’ 
This script will process waveform data and compute quality control (QC) criteria for every waveform.  
Within this script, there are sever user defined variables. See comments in this script for details.
Results will be stored in the ‘Keep_Data_lists’ directory. Plots of waveforms, signal spectra, and noise spectra will be stored in the ‘Noise_plot’ directory.

- ‘teleseismic_eq_cc_process.py’
This script will process cross-correlations of seismograms for the same event recorded by a station pair (for data that pass QC criteria). Within this script, there are sever user defined variables. See comments in this script for details. Furthermore, this script utilizes python joblib for parallel computations. The number of CPUs used can be modified according to your system in order to improve performance of this script.  Results will be stored in the ‘CC_process’ directory.   

- ‘teleseismic_eq_cc_stack_von_mises.py’  
This script will compute properties of von Mises distributions for every station pair and every desired period.  
Within this script, there are sever user defined variables. See comments in this script for details. Furthermore, this script utilizes python joblib for parallel computations. The number of CPUs used can be modified according to your system in order to improve performance of this script. 
Results will be stored in the ‘Phase_VM’ directory. Plots of histograms of phase delay measurements will be stored in the ‘Phase_VM_plots’ directory. 




Probabilistic inversion:   

The codes for performing the probabilistic inversion of phase spectra are written in Fortran, and are contained within the ‘1D_Inversion’ directory. There, you will find a makefile necessary to compile the Fortran code.  You may have to modify this makefile to suit your needs and system requirements.  In your terminal, in the directory with the codes and makefile, type:

make clean   
make

 This will produce two executables:   

- ‘transd_inversion_1d’ 
This is the main inversion code. 
To call this code, type in your terminal: 

./transd_inversion_1d

Input parameters for the inversion are contained within the ‘inversion_parameters.in’ file. See comments in this file for details. Furthermore, this Fortran code uses OpenMP for parallel computations. The number of CPUs used can be modified according to your system in order to improve performance of this script. This is specified in the ‘inversion_parameters.in’ file. 
Inversion results for all station pairs will be stored in the ‘results_1d’ directory.   


- ‘grid_models_1d’ 
This code will use a subset of the ensemble of rjMcMC samples to calculate marginal histograms. Importantly, this code computes the depth marginal that can be used for plotting. Input parameters for the inversion are contained within the ‘inversion_parameters.in’ file. Again, importantly, the ‘cutoff’ value defined the number of initial samples that were SAVED to disc that you wish to exclude from the marginal (as burn-in; see GJI paper for details).  Results for all station pairs will be stored in the ‘results_1d’ directory.   


- ‘1D_plot_profile.py’ 
This script provides an example for how to read the inversion results into python, and recreates the 1D inversion results (example) figure shown in the GJI paper. This script is provided as an example, ‘as-is’.  
