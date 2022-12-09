

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
import numpy.fft as ft

import obspy.core
import os.path
import pickle
from math import pi, sin
from obspy.core import read, Stream, Trace, AttribDict

# from obspy.fdsn import Client
from obspy.clients.fdsn import Client

from obspy.geodetics import gps2dist_azimuth as epi
from obspy.geodetics import kilometers2degrees as k2d

from scipy.fftpack import fft, fftfreq, fftshift
from scipy import signal
from scipy.signal import hilbert

from scipy.stats import norm

from joblib import Parallel, delayed


# Wrapper function that can be called for parallel processing
# #####################################################################################
def parallel_process_cc(iev, ev, nev, station1, station2, network1, network2, \
                        datapath_iris1, datapath_iris2, epi_dist1, epi_dist2, \
                        vr_min, vr_max):

    print('Event: ' + str(iev+1) + ' out of ' + str(nev) + \
           '  between ' + station1 + '<-->' + station2 )

    # Extract time, coordinates and depth of events
    time = ev.origins[0].time

    # Define time stamp
    yr = str(time.year).zfill(4)
    mn = str(time.month).zfill(2)
    dy = str(time.day).zfill(2)
    hr = str(time.hour).zfill(2)
    mi = str(time.minute).zfill(2)
    sc = str(time.second).zfill(2)
    tstamp = yr+mn+dy+'_'+hr+mi+sc

    filename1 = datapath_iris1 + '/' + tstamp + '.pkl'
    filename2 = datapath_iris2 + '/' + tstamp + '.pkl'

    # Check if files exist for station 1
    if os.path.isfile(filename1): 
        # Load data 1
        try:
            # Load file
            infile = open(filename1,'rb')
            (ev1,zz1) = pickle.load(infile)
            infile.close()
        except:
            print("Problem loading file: " + filename1)
            return []
    
    else:
        print("Data unavailable: " + network1 + '_' + \
                station1 + ' on ' + tstamp )
        return []

    # Check if files exist for station 2
    if os.path.isfile(filename2): 
        # Load data 1
        try:
            # Load file
            infile = open(filename2,'rb')
            (ev2,zz2) = pickle.load(infile)
            infile.close()
        except:
            print("Problem loading file: " + filename2)
            return []
    
    else:
        print("Data unavailable: " + network2 + '_' + \
                station2+ ' on ' + tstamp )
        return []

    # Defining time windows
    tw11 = time + epi_dist1/vr_max
    tw12 = time + epi_dist1/vr_min

    tw21 = time + epi_dist2/vr_max
    tw22 = time + epi_dist2/vr_min

    tww1 = min(tw11,tw21)
    tww2 = max(tw12,tw22)

    try:
        # Phase spectrum of cross-correlation (with and without additional signal isolation)
        if (epi_dist1 < epi_dist2):
            phases_iso, amplitudes_iso = get_phases_isolate(zz1,zz2,tww1,tww2)
            tcorr, phases, amplitudes = get_phases(zz1,zz2,tww1,tww2)

        else:
            phases_iso, amplitudes_iso = get_phases_isolate(zz2,zz1,tww1,tww2)
            tcorr, phases, amplitudes = get_phases(zz2,zz1,tww1,tww2)

    except:
        return []


    # Package of results to return
    structure = [phases, amplitudes, \
                 phases_iso, amplitudes_iso, \
                 np.abs(epi_dist1 - epi_dist2), ev]

    return structure
# #####################################################################################
# #####################################################################################

# Function to esimate phase spectrum
# #####################################################################################
def get_phases(z1,z2,tww1,tww2):

    z1c = z1.copy()
    z2c = z2.copy()

    z1c.trim(tww1,tww2)
    z2c.trim(tww1,tww2)

    t_corr = signal.correlate(z2c.data,z1c.data,mode='full')

    if (len(t_corr) < (2*N + 1)):
        print("Phase velocity estimation failed --> CC too short")
        return([])

    # N = 1200 # Number of samples 
    t_corr = t_corr[len(z1c.data):len(z1c.data)+N]

    yf = fft(t_corr)

    phases = np.angle(yf)

    amplitudes = np.abs(yf)

    return t_corr, phases, amplitudes
# #####################################################################################

# Function to estimate phase spectrum with additional signal isolation
# #####################################################################################
def get_phases_isolate(z1,z2,tww1,tww2):

    z1c = z1.copy()
    z2c = z2.copy()

    z1c.trim(tww1,tww2)
    z2c.trim(tww1,tww2)

    t_corr = signal.correlate(z2c.data,z1c.data,mode='full')

    if (len(t_corr) < (2*N + 1)):
        print("Phase velocity estimation failed --> CC too short")
        return([])

    # N = 1200 # Number of samples
    t_corr = t_corr[len(z1c.data):len(z1c.data)+N]
    xf = fftfreq(N, 1.0)

    phases = np.zeros(len(keep_periods_SNR))
    amplitudes = np.zeros(len(keep_periods_SNR))

    times = np.arange(1,N+1)

    # Loop over frequencies of interest
    for iT in range(len(keep_periods_SNR)):
        
        # Narrow pass band around frequency of interest
        fl = 0.85*(1./keep_periods_SNR[iT])
        fh = 1.15*(1./keep_periods_SNR[iT])
        z1c2 = z1c.copy()
        z1c2.data = np.copy(t_corr)
        z1c2.filter('bandpass', freqmin=fl, freqmax=fh, corners=2, zerophase=True)

        # Calculate signal envelope (and locate maximum)
        za1 = np.abs(hilbert(z1c2.data))
        imax = np.argmax(za1) + 1

        # Time domain window for signal isolation 
        g_window = norm.pdf(times,loc=imax,scale=keep_periods_SNR[iT])

        # Filtered CC with time domain windowing
        t_corr_iso = z1c2.data*g_window

        yf = fft(t_corr_iso)
        
        # Only keep amplitude and phase at desired period
        ix = np.argmin( np.abs(xf - 1./keep_periods_SNR[iT]) )
        phases[iT] = np.angle(yf[ix])
        amplitudes[iT] = np.abs(yf[ix])

    return phases, amplitudes
# #####################################################################################








########################
# User defined variables
########################


# Theoretical velocity limit for download time window
vr_min = 2.5
vr_max = 4.5

# Threshold for signal-to-noise ratio for entire frequency band
snr_threshold = 2

# Number of samples to keep in CC
# You can experiment with this, but 1200 gives exact period values as listed in the 'keep' array
# Otherwise there can be a slight period "mismatch" in the fft
# The downside is that this will throw away data for shorter recordings (close events)
# We could experiment with zero-padding CCs, or more-sophisticated period matching (something to discuss)  
N = 1200

# Periods for isolated CC processing
keep_periods_SNR = [120.,100.,80.,60.,50.,40.,30.,25.,20.,15.]

# Path to keep track of QC info
good_data_record_path = "Keep_Data_lists"

# Path to save computed CC
cc_record_path = "CC_process"


# #####################################################################################

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
# #####################################################################################

# Create path to CC results
if not os.path.isdir(cc_record_path): 
    print('Path to '+cc_record_path+' doesn`t exist - creating it')
    os.makedirs(cc_record_path)

# Loop through stations
for istation1 in range(0,len(stations)):

    # Some relic code for testing (if you want to test a specific station)
    # if (stations[istation1] != "CRAG") :
    #     continue


    # Extract pre-computed SNR information for this station
    filename = good_data_record_path + '/' + networks[istation1] + '_' + stations[istation1] + '.pkl'

    infile = open(filename,'rb')
    structure = pickle.load(infile)
    infile.close()

    evs1 = structure[0]
    ev_dists1 = structure[1]
    data_source1 = structure[2]
    ev_times1 = structure[3]
    snrs1 = structure[4]
    snrs_isolate1 = structure[5]

    # Data path for station1
    datapath_iris1 = 'Data/' + networks[istation1] + '_' + stations[istation1]

    # Loop through stations
    for istation2 in range(istation1,len(stations)):

        if (istation1 == istation2):
            continue
        
        # Some relic code for testing (if you want to test a specific station)
        #if (stations[istation2] != "R33M") :
        #    continue

        # Extract pre=computed SNR information for this station
        filename = good_data_record_path + '/' + networks[istation2] + '_' + stations[istation2] + '.pkl'

        infile = open(filename,'rb')
        structure = pickle.load(infile)
        infile.close()

        evs2 = structure[0]
        ev_dists2 = structure[1]
        data_source2 = structure[2]
        ev_times2 = structure[3]
        snrs2 = structure[4]
        snrs_isolate2 = structure[5]

        # Get previously-saved catalog for this 2-station combination
        filename = 'CAT/' +  networks[istation1] + '_' + stations[istation1] + '_' + \
                             networks[istation2] + '_' + stations[istation2] + '.pkl'

        if not os.path.isfile(filename):
            print("No catalog for station pair: " + \
                  networks[istation1] + '_' + stations[istation1] + '<-->' + \
                  networks[istation2] + '_' + stations[istation2] ) 
            continue


        infile = open(filename,'rb')
        sub_cat = pickle.load(infile)
        infile.close()

        # IRIS data path for station2
        datapath_iris2 = 'Data/' + networks[istation2] + '_' + stations[istation2]


        # Find sub EQ catalog that matches pre-calculated SNR threshold
        sub_cat_keep = []
        epi_dists1_sub = []
        epi_dists2_sub = []
        for ev in sub_cat:

            # Check if event satisfies snr criteria
            if not ( (ev in evs1) and (ev in evs2) ):
                continue

            # Index of event
            iev1 = evs1.index(ev)
            iev2 = evs2.index(ev)

            epi_dists1_sub.append(ev_dists1[iev1])
            epi_dists2_sub.append(ev_dists2[iev2])

            # Update subcatalog
            sub_cat_keep.append(ev)

        nev = len(sub_cat_keep)
        # nev is the number of events that satisfy broadband SNR criteria. 
        # You can decide to apply criteria for minimum number of events for a station pair.  
        # if (nev < 5) :
        #     continue


        # Parallel code for CC processing (n_jobs is number of CPUs you have available)
        # Parallel code divides # of events for station pair into n_jobs
        parallel_results = Parallel(n_jobs=4)(delayed(parallel_process_cc)\
            (iev, ev, nev, stations[istation1], stations[istation2],\
             networks[istation1], networks[istation2], \
             datapath_iris1, datapath_iris2, epi_dists1_sub[iev], epi_dists2_sub[iev], \
             vr_min, vr_max) for iev, ev in enumerate(sub_cat_keep) )


        # Book keeping of results of parallel code
        phases_list = []
        amplitudes_list = []
        station_dist_diffs = []
        ev_list = []
        phases_list_iso = []
        amplitudes_list_iso = []


        for result in parallel_results:

            # Parallel code failed (for specific event)
            if result == []:
                continue

            phases         = result[0]
            amplitudes     = result[1]
            phases_iso     = result[2]
            amplitudes_iso = result[3]
            dist_diff      = result[4]
            ev             = result[5]

            phases_list.append(phases)
            amplitudes_list.append(amplitudes)
            station_dist_diffs.append(dist_diff)
            phases_list_iso.append(phases_iso)
            amplitudes_list_iso.append(amplitudes_iso)
            ev_list.append(ev)

            print(phases_iso[0:5])

        if (len(phases_list) == 0):
            print("No CC processed for: " + \
                  networks[istation1] + '_' + stations[istation1] + '<-->' + \
                  networks[istation2] + '_' + stations[istation2] ) 
            continue
            
        # Distance between stations 
        epi_dist3, az3, baz3 = epi(stlats[istation1], stlons[istation1], stlats[istation2], stlons[istation2])
        epi_dist3 /= 1000 


        # Save results of CC processing for this station pair
        filename = cc_record_path + '/' + networks[istation1] + '_' + stations[istation1] + '_' + \
                                          networks[istation2] + '_' + stations[istation2] + '_CC.pkl'

        structure = ( epi_dist3, ev_list, station_dist_diffs, \
                    phases_list, phases_list_iso, \
                    amplitudes_list, amplitudes_list_iso )

        # Save CC results 
        outfile = open(filename,'wb')
        pickle.dump(structure,outfile)
        outfile.close()

































































































