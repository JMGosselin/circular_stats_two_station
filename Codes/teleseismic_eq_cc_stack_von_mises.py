

# Copyright 2022 by Jeremy M. Gosselin
# Department of Geoscience
# University of Calgary

# This program is free software: 
# you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) 
# any later version.

# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
# for more details.

# You should have received a copy of the GNU General Public License along 
# with this program. 
# If not, see <https://www.gnu.org/licenses/>.


# Import modules and functions
#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import obspy.core
import os.path
import pickle
from math import pi, sin

# from scipy.fftpack import fft, fftfreq, fftshift
import scipy.stats as stats
from joblib import Parallel, delayed
from scipy.special import i0  
from scipy.stats import vonmises




def process_phase_VM(istation1,stations,networks,cc_record_path,good_data_record_path,phase_stack_path,phase_plot_path,keep_periods,plot_flag):

    # Extract pre-computed SNR information for station1
    filename_snr = good_data_record_path + '/' + networks[istation1] + '_' + stations[istation1] + '.pkl'
    infile = open(filename_snr,'rb')
    structure = pickle.load(infile)
    infile.close()
    evs1 = structure[0]
    ev_dists1 = structure[1]
    data_source1 = structure[2]
    ev_times1 = structure[3]
    snrs1 = structure[4]
    snrs_isolate1 = structure[5]

    # Keeping track of number of stacked phase spectra
    count = 0

    # Loop over stations to find CC between station1 and station2
    for istation2 in range(istation1,len(stations)):

        if (istation1 == istation2):
            continue

        # Relic for testing specific station pairs
        # if (stations[istation2] != "T35M") :
        #     continue
        # if (stations[istation2] != "R33M") :
        #     continue

        # Check input file
        filename = cc_record_path + '/' + networks[istation1] + '_' + stations[istation1] + '_' + \
                                          networks[istation2] + '_' + stations[istation2] + '_CC.pkl'

        if not os.path.isfile(filename):
            print("No catalog for station pair: " + \
                  networks[istation1] + '_' + stations[istation1] + '<-->' + \
                  networks[istation2] + '_' + stations[istation2] ) 
            continue


        # Load previously saved CC processing results for this station pair
        infile = open(filename,'rb')
        structure = pickle.load(infile)
        infile.close()
        epi_dist3           = structure[0]
        ev_list             = structure[1]
        station_dist_diffs  = structure[2]
        phases_list         = structure[3]
        phases_list_iso     = structure[4]
        amplitudes_list     = structure[5]
        amplitudes_list_iso = structure[6]

        # You can put in extra QC criteria here, to avoid re-running the CC code.

        # if (len(phases_list_iso) < 5):
        #     print('Too few events for station pair: '+ \
        #         networks[istation1] + '_' + stations[istation1] + '<-->' + \
        #           networks[istation2] + '_' + stations[istation2] )
        #     continue

        # Extract pre=computed SNR information for station2
        filename_snr = good_data_record_path + '/' + networks[istation2] + '_' + stations[istation2] + '.pkl'
        infile = open(filename_snr,'rb')
        structure = pickle.load(infile)
        infile.close()
        evs2 = structure[0]
        ev_dists2 = structure[1]
        data_source2 = structure[2]
        ev_times2 = structure[3]
        snrs2 = structure[4]
        snrs_isolate2 = structure[5]


        # Output file
        filename_out = phase_stack_path + \
                   networks[istation1] + '_' + stations[istation1] + '__' + \
                   networks[istation2] + '_' + stations[istation2] + '_phase.txt'
        outfile = open(filename_out,'w') 
        outfile.write('{} \n'.format(len(keep_periods)))
        outfile.write('{} \n'.format(epi_dist3))
        
        # Lists of individual vM properties (locations and concentrations) for each period
        kappas = []
        mus = []

        period_plot_list = keep_periods.copy()



        ## Phase diagram plot
        if plot_flag:
            fig = plt.figure(figsize=(8,14))


        # Loop over desired periods
        for ipp, pp in enumerate(period_plot_list):
            
            if plot_flag:
                ####################################################################################################
                # You will need to manually change this depending on the number of periods you decided to use 
                # In order to make the number of subplots match your keep periods list 
                ####################################################################################################
                ax = fig.add_subplot(5,2,len(period_plot_list) - ipp)


            # Keep track of isolated phases at desired period
            temp_phases = []

            # Loop over events in station pair
            for iphase, phases in enumerate(phases_list_iso):

                # Event that passed initial QC for station pair
                ev = ev_list[iphase]

                # Event index in station 1 list
                iev1 = evs1.index(ev)

                # Event index in station 2 list
                iev2 = evs2.index(ev)
                
                ###################################################################################################
                ###################################################################################################
                ###################################################################################################

                # Uncomment if you want to apply additional QC (period-dependent SNR QC) 
                # We can modify this later
                # if ( (snrs_isolate1[iev2][ipp] > snr_threshold) and (snrs_isolate1[iev1][ipp] > snr_threshold) ):
                #     temp_phases.append(phases[ipp])
                
                # Comment if you uncomment the above
                temp_phases.append(phases[ipp])

                ###################################################################################################
                ###################################################################################################
                ###################################################################################################



            # Check if list is empty (only applies when using additional QC)
            if temp_phases == []:
                kappa = -1.0
                mu = 0.0
            else:
                kappa, mu, scale = vonmises.fit(temp_phases, fscale=1)
                
                # Plotting stuff
                if plot_flag:
                    counts, bins, patches = ax.hist(temp_phases, 20, density=True)
                    ax.set_yticks([])
                    ax.text(1.0, 1.005, str(pp) +' s',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes) 
                    ax.set_ylim([0.,1.5*np.max(counts)])
                    ax.set_xlim([-np.pi, np.pi])
                    ax.grid()

                    ####################################################################################################
                    # If you modify the number of subplots, then you will have to edit the labelling below
                    ####################################################################################################
                    if ipp in [1,3,5,7,9]:
                        ax.set_ylabel('Probability')
                    if ipp in [0,1]:
                        ax.set_xlabel('Phase (rad)')
                        ax.set_xticks([-np.pi,-np.pi/2, 0, np.pi/2, np.pi])
                        ax.set_xticklabels(['$-\pi$','$-\pi/2$', '0', '$\pi/2$', '$\pi$'])
                    else:
                        ax.set_xticks([-np.pi,-np.pi/2, 0, np.pi/2, np.pi])
                        ax.set_xticklabels([])                        


            kappas.append(kappa)
            mus.append(mu)

            outfile.write('{} {} {} \n'.format(pp, mu, kappa))

        outfile.close() 

        # Plotting stuff
        if plot_flag:
            fig_file = phase_plot_path + \
                  networks[istation1] + '_' + stations[istation1] + '__' + \
                  networks[istation2] + '_' + stations[istation2] + '_phase_example.png'

            plt.savefig(fig_file,bbox_inches='tight',dpi=300)
            # plt.show()
            plt.close('all')


        print('Stacked station pair: '+ \
            networks[istation1] + '_' + stations[istation1] + '<-->' + \
              networks[istation2] + '_' + stations[istation2] )

        count = count + 1
    

    return istation1, istation2, count


#####################################################################################
#####################################################################################
#####################################################################################


########################
# User defined variables
########################

snr_threshold = 2

# Periods to extract CC (same as CC processing code)
keep_periods_SNR = [120.,100.,80.,60.,50.,40.,30.,25.,20.,15.]

# vel_axis = np.linspace(3.005,4.995,200)
# slow_axis = np.linspace(0.2005,0.3995,200)

# Results directory for noise spectra code
good_data_record_path = "Keep_Data_lists"

# Results directory for CC processing
cc_record_path = "CC_process"

# Results directory for this script
# Where to save stacked noise spectra and figures
phase_stack_path = 'PHASE_VM/'
phase_plot_path = 'PHASE_VM_plots/'



#####################################################################################
#####################################################################################
#####################################################################################

# Create path to VM results
if not os.path.isdir(phase_stack_path): 
    print('Path to '+phase_stack_path+' doesn`t exist - creating it')
    os.makedirs(phase_stack_path)

# Create path to VM plots
if not os.path.isdir(phase_plot_path): 
    print('Path to '+phase_plot_path+' doesn`t exist - creating it')
    os.makedirs(phase_plot_path)


# Create empty lists
networks = []
stations = []
stlats   = []
stlons   = []
start_times = []
end_times = []  


# Read in station information
with open('station_list_test.txt', 'r') as f:
    # do things with your file
    line = f.readline()

    while line:
        data = line.split(",")
        
        networks.append(data[0])
        stations.append(data[1])
        stlats.append(float(data[2]))
        stlons.append(float(data[3]))
        #stelevs.append(float(data[4]))
        start_times.append(data[5])
        end_times.append(data[6])
        #sources.append(data[7]) 

        line = f.readline()


# #####################################################################################
# #####################################################################################


# Do you want to generate plots of phase histograms for every station pair?
plot_flag = True



# Parallel code that loops over stations, and searches for CC results to stack
results = Parallel(n_jobs=4)(delayed(process_phase_VM)\
    (istation1,stations,networks,cc_record_path,\
                                 good_data_record_path,\
                                 phase_stack_path,\
                                 phase_plot_path,\
                                 keep_periods_SNR,plot_flag) for istation1 in range(0,len(stations)) )
##
# If the call to the parallel function does not work for you, then you can easily run it in serial by uncommenting below  
##

# results = []
# for istation1 in range(0,len(stations)):
#     result = process_phase_VM(istation1,stations,networks,cc_record_path,\
#                                  good_data_record_path,\
#                                  phase_stack_path,\
#                                  phase_plot_path,\
#                                  keep_periods_SNR,plot_flag)
#     results.append(result)



# Save stacking output (just keeps track of final number of events per station pair) 
outfile = open('VM_output.pkl','wb')
pickle.dump(results,outfile)
outfile.close()


