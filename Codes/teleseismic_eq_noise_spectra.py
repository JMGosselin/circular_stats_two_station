
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
from obspy.core import read

from obspy.clients.fdsn import Client

from obspy.geodetics import gps2dist_azimuth as epi
from obspy.geodetics import kilometers2degrees as k2d

from scipy.fftpack import fft, fftfreq, fftshift

from matplotlib.ticker import ScalarFormatter, NullFormatter


#####################################################################################
# Calculate Signal to Noise ratio
def zsnr(trZ, time, t1, t2):

    # Copy Z trace to signal and noise traces
    trSig = trZ.copy()
    trNze = trZ.copy()

    # Trim twin seconds around arrival
    trSig.trim(time + t1, time + t2)
    trNze.trim(time,time + t1)

    # Calculate root mean square (RMS)
    srms = np.sqrt(np.mean(np.square(trSig.data)))
    nrms = np.sqrt(np.mean(np.square(trNze.data)))

    # Calculate signal/noise ratio in dB
    snr = 10*np.log10(srms*srms/nrms/nrms)

    return snr, trNze, trSig 
#####################################################################################

#####################################################################################
def get_snr_isolate(trZ, time, t1, t2, keep_periods):

    snrs = []

    for iT, T in enumerate(keep_periods):

        # Upper and lower bounds of narrow-bandpass
        fl = 0.85*(1./T)
        fh = 1.15*(1./T)
        
        # Copy Z trace
        trc = trZ.copy()
        
        # Filter
        trc.filter('bandpass', freqmin=fl, freqmax=fh, corners=2, zerophase=True)
        
        # Compute SNR
        snr, trNze, trSig = zsnr(trc, time, t1, t2)

        snrs.append(snr)

    return snrs
#####################################################################################






########################
# User defined variables
########################


# Theoretical velocity limit for download time window
vr_min = 2.5
vr_max = 4.5

# lat/lon of centre point of study region (simplifies catalog searches)
GC  = [54.4588, -124.2944] 

# Only download data if station are at least this far apart (in km)
min_station_dist = 200.

# Only download data for station pairs with at least this much overlap in recording period
min_years = 2.

# Arbitrary amplitude scaling for plotting purposes (can modify to suit plotting needs)
scale = 200.

# Threshold for signal-to-noise ratio for entire frequency band
snr_threshold = 2

# Periods for plotting
N = 1200
xf = fftfreq(N, 1.0)
keep_periods = (1./xf)[1:101]
freqs = 1./keep_periods

# Periods for SNR processing
keep_periods_SNR = [120.,100.,80.,60.,50.,40.,30.,25.,20.,15.]

# Path to keep track of QC info
good_data_record_path = "Keep_Data_lists"

# Data locations
resp_path         = "FTP_RESP/"

# Station response file locations
save_path_daylong = "Daylongs/"

# Path to save figures
plot_path = "Noise_plots/"





# #####################################################################################
# #####################################################################################
# #####################################################################################


# Create empty lists for station info...
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

# Create path to QC records
if not os.path.isdir(good_data_record_path): 
    print('Path to '+good_data_record_path+' doesn`t exist - creating it')
    os.makedirs(good_data_record_path)

# Create path to plots directory 
if not os.path.isdir(plot_path): 
    print('Path to '+plot_path+' doesn`t exist - creating it')
    os.makedirs(plot_path)

# Loop through stations and search for events
for istation1 in range(0,len(stations)):

    # Relic for testing purposes
    # if (stations[istation1] != "PGC") :
    #     continue

    # Create figure and subplot
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121)

    # Initialize maximum source-station distance and list of noise traces
    max_dist = 0           # For scalling axes in plots 
    evs = []               # Events that satisfy QC
    noise_traces = []      # Pre-arrival recordings
    signal_traces = []     # post-arrival recordings
    ev_dists = []          # source-station distances
    data_source = []       # Relic (was nedded when combining data from diferent sources - e.g. IRIS vs. CNSN FTP)
    ev_times = []          # Event onset times
    snrs = []              # Broadband signal-to-noise ratios
    time_windows = []      # Absolute noise time window
    snrs_isolate = []      # SNR spectrum

    # Data path for station 
    datapath_iris = 'Data/' + networks[istation1] + '_' + stations[istation1] 

    # Does data path exist for this station?
    if os.path.isdir(datapath_iris): 

        # Loop through files in station path
        for filename in os.listdir(datapath_iris):
            full_filename = os.path.join(datapath_iris, filename)

            try:
                # Load file
                infile = open(full_filename,'rb')
                (ev,zz) = pickle.load(infile)
                infile.close()
            except:
                print("Problem loading file: " + full_filename)
                continue

            # Extract time, coordinates and depth of events
            time = ev.origins[0].time
            lat = ev.origins[0].latitude
            lon = ev.origins[0].longitude
            dep = ev.origins[0].depth

            # Calculate epicentral distance from EQ to station1
            epi_dist1, az1, baz1 = epi(lat, lon, stlats[istation1], stlons[istation1])
            epi_dist1 /= 1000 

            # Update max distance
            if epi_dist1 > max_dist:
                max_dist = epi_dist1 

            # Check SNR
            snr, noise, signal = zsnr(zz, time, epi_dist1/vr_max, epi_dist1/vr_min)
            if ( (snr < snr_threshold) or np.isnan(snr) ):
                continue

            # SNR are unique periods
            snr_isolate = get_snr_isolate(zz, time, epi_dist1/vr_max, epi_dist1/vr_min, keep_periods_SNR)

            # Keep noise trace for later processing 
            noise_traces.append(noise)
            signal_traces.append(signal)
            ev_dists.append(epi_dist1)
            data_source.append('IRIS')
            ev_times.append(time)
            snrs.append(snr)
            evs.append(ev) 
            time_windows.append([time, time + epi_dist1/vr_max])
            snrs_isolate.append(snr_isolate)

            # Plot waveform
            plt.plot(zz.times(),scale*zz.data/np.max(zz.data) + epi_dist1,'b-',alpha=0.5,linewidth=0.1) 

    # Plot surface wave expected time window
    plt.plot([0,max_dist/vr_max],[0,max_dist],'k:')
    plt.plot([0,max_dist/vr_min],[0,max_dist],'k:')

    # Other plot stuff
    ax.set_ylim([0,max_dist+scale])
    ax.set_xlabel('Time since event (s)')
    ax.set_ylabel('Source-station distance (km)')
    plt.title(stations[istation1])

    # Initialize list of noise spectra
    noise_spectra = []
    signal_spectra = []

    # Create subplot
    ax = fig.add_subplot(222)

    ###########################################################
    ###########################################################
    # Loop through noise traces, compute and plot power spectra
    for ii, noise in enumerate(noise_traces):
        
        # Minimum noise window length (window is too short if station is close to event)        
        if (len(noise.data) < N):
            noise_spectra.append([0])
            continue

        # Only want the first N seconds of noise
        noise_trim = noise.data[0:N]
        ps = np.abs(np.fft.fft(noise_trim))**2
        
        # Plot noise spectrum at select periods
        ax.loglog(keep_periods,ps[1:101],'b-',alpha=0.5,linewidth=0.1)

        noise_spectra.append(ps[1:101])

    ax.set_ylabel('Noise Power Spectrum')
    ax.set_xlim([20,200])
    ax.set_ylim([10**(-19),1])

    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())

    ax.set_xticks([20,30,40,50,60,70,80,90,100,200])
    ax.set_xticklabels(['20','','','50','','','','','100','200'])

    # Create subplot
    ax = fig.add_subplot(224)

    ###########################################################
    ###########################################################
    # Loop through noise traces, compute and plot power spectra
    for ii, signal in enumerate(signal_traces):
        
        # Minimum noise window length (window is too short if station is close to event)
        noise = noise_traces[ii]
        if (len(noise.data) < N):
            signal_spectra.append([0])
            continue

        # Signal power spectrum
        N_sig = len(signal.data)
        xf_sig = fftfreq(N_sig, 1.0)
        keep_periods_sig = 1./xf_sig
        ps = np.abs(np.fft.fft(signal.data))**2

        # Find indeces of nearest matching periods
        ps_keep = np.zeros(len(freqs))
        for ip, period in enumerate(keep_periods):
            imin = np.argmin(np.abs(keep_periods_sig - period))
            ps_keep[ip] = ps[imin]

        # Plot power spectrum at select periods (same as noise spectrum)
        ax.loglog(keep_periods,ps_keep,'b-',alpha=0.5,linewidth=0.1)

        signal_spectra.append(ps_keep)


    ###########################################################
    
    ax.set_ylabel('Signal Power Spectrum')
    ax.set_xlabel('Period (s)')
    ax.set_xlim([20,200])
    ax.set_ylim([10**(-19),1])

    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())

    ax.set_xticks([20,30,40,50,60,70,80,90,100,200])
    ax.set_xticklabels(['20','','','50','','','','','100','200'])

    plt.savefig('Noise_plots/' + stations[istation1]+'_noise_spectra.png',bbox_inches='tight',dpi=300)
    plt.close('all')

    # Save QC info for later
    filename = good_data_record_path + '/' + networks[istation1] + '_' + stations[istation1] + '.pkl'

    structure = (evs, \
                ev_dists, \
                data_source, \
                ev_times, \
                snrs,\
                snrs_isolate )

    # Save QC data 
    outfile = open(filename,'wb')
    pickle.dump(structure,outfile)
    outfile.close()

    # Relic for testing...
    # fig = plt.figure()

    # for snr in snrs_isolate:
    #     plt.plot(keep_periods_SNR,snr,'r-')

    # plt.show()






