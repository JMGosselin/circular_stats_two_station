
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
import matplotlib.pyplot as plt
import obspy.core
import os.path
import pickle
from math import pi, sin
from obspy import UTCDateTime
from obspy.core import read, Stream, Trace, AttribDict

# from obspy.fdsn import Client
from obspy.clients.fdsn import Client

from obspy.geodetics import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d


########################
# User defined variables
########################

# Theoretical velocity limit for download time window
vr_min = 2.5

# lat/lon of centre point of study region (simplifies catalog searches)
GC  = [54.4588, -124.2944] 

# Only download data if station are at least this far apart (in km)
min_station_dist = 200.

# Only download data for station pairs with at least this much overlap in recording period
min_years = 2.

# Total timde window for data search (this is refined automatically for each station pair)
tmin = UTCDateTime("2000-01-01T00:00:00")
tmax = UTCDateTime("2020-12-01T00:00:00")
tmax_master = tmax

# Angular threshold for accepting events aligned with a station pair (something between 5-10 degrees)
angle_thresh = 8.


########################
# User defined variables
########################









# Uncomment this block if this is not the first time running this script
# The script keeps track of time windows that are not available on IRIS (so re-running the code is faster)
#####################################################################################
# filename = 'unavailable_files.pkl'
# infile = open(filename,'rb')
# (unavailable_files) = pickle.load(infile)
# infile.close()
#####################################################################################

# Comment this if this is not the first time running this script
unavailable_files = []
#####################################################################################


###############################
# End of user defined variables
###############################





# Create empty lists for station data...
networks = []
stations = []
stlats   = []
stlons   = []
start_times = []
end_times = []  


# Read in station information
with open('station_list_test.txt', 'r') as f:
    line = f.readline()

    while line:
        data = line.split(",")
        
        # This tells you order of station info in each row of station file
        networks.append(data[0])
        stations.append(data[1])
        stlats.append(float(data[2]))
        stlons.append(float(data[3]))
        #stelevs.append(float(data[4]))
        start_times.append(data[5])
        end_times.append(data[6])
        #sources.append(data[7]) 

        line = f.readline()


# Get download client
client = Client()

#####################################################################################
# Get total earthquake catalogue using project start and end time
cat = client.get_events(starttime=tmin, endtime=tmax, minmagnitude=5.5, \
   latitude=GC[0], longitude=GC[1], minradius=20., maxradius=150.,  maxdepth=100.)

# Catalog folder
datapath = 'CAT/' 
if not os.path.isdir(datapath): 
    print('Path to '+datapath+' doesn`t exist - creating it')
    os.makedirs(datapath)

# Write catalog to disc
outfile = open(datapath+'master.pkl','wb')
pickle.dump((cat),outfile)
outfile.close()

# Comment out the above, and uncomment this, if you already have a catalog file
#####################################################################################
# filename = datapath+'master.pkl'
# infile = open(filename,'rb')
# (cat) = pickle.load(infile)
# infile.close()
#####################################################################################



# Loop through stations and search for events
for istation1 in range(0,len(stations)):

    # Some relic code for testing (if you want to test a specific station)
    # if (stations[istation1] != "CRAG") :
    #     continue

    for istation2 in range(istation1,len(stations)):

        # Sanity check
        if (istation1 == istation2):
            continue

        # Some relic code for testing (if you want to test a specific station)
        # if (stations[istation2] != "R33M") :
        #     continue

        # Start and end times for station1 recordings
        tmin1 = UTCDateTime(start_times[istation1])
        tmax1 = UTCDateTime(end_times[istation1] )
  
        # Start and end times for station2 recordings
        tmin2 = UTCDateTime(start_times[istation2])
        tmax2 = UTCDateTime(end_times[istation2] )

        # Update tmin and tmax for station pair (overlapping deployment)
        tmin = max(tmin1,tmin2)
        tmax = min(tmax1,tmax2)

        # Currently operating stations sometimes have weird deployment end times
        tmax = min(tmax_master,tmax)

        # No overlapping data
        if ( tmin > tmax ):
            continue

        # Less than 2 years of overlap data
        if ( (tmax - tmin) < min_years*365.*24.*60.*60. ):
            continue

        # Get station separation and orientation
        epi_dist3, az3, baz3 = epi(stlats[istation1], stlons[istation1], stlats[istation2], stlons[istation2])
        epi_dist3 /= 1000    
        
        # Stations are too close to each other
        if (epi_dist3 < min_station_dist):
            continue 

        # Define subset earthquake catalog for station pair
        sub_cat = []

        # Number of events in subset catalog
        count = 0

        # Loop through all events in master catalog
        for iev, ev in enumerate(cat):

            # Extract time, coordinates and depth of events
            time = ev.origins[0].time
            lat = ev.origins[0].latitude
            lon = ev.origins[0].longitude
            dep = ev.origins[0].depth


            # Define time stamp
            yr = str(time.year).zfill(4)
            mn = str(time.month).zfill(2)
            dy = str(time.day).zfill(2)
            hr = str(time.hour).zfill(2)
            mi = str(time.minute).zfill(2)
            sc = str(time.second).zfill(2)
            tstamp = yr+mn+dy+'_'+hr+mi+sc

            # Check event time
            if ( (time < tmin) or (time > tmax) ):
                continue

            # Calculate epicentral distance from EQ to station1
            epi_dist1, az1, baz1 = epi(lat, lon, stlats[istation1], stlons[istation1])
            epi_dist1 /= 1000    

            # Calculate epicentral distance from EQ to station2
            epi_dist2, az2, baz2 = epi(lat, lon, stlats[istation2], stlons[istation2])
            epi_dist2 /= 1000    

            # Update station separation and orientation
            # Station1 is closer 
            if (epi_dist1 < epi_dist2):
                # Calculate epicentral distance from station1 to station2
                epi_dist3, az3, baz3 = epi(stlats[istation1], stlons[istation1], stlats[istation2], stlons[istation2])
                epi_dist3 /= 1000    

            else:
                # Calculate epicentral distance from station2 to station1
                epi_dist3, az3, baz3 = epi(stlats[istation2], stlons[istation2], stlats[istation1], stlons[istation1])
                epi_dist3 /= 1000    

            # Distance to midpoint of study region
            epi_dist4, az4, baz4 = epi(lat, lon, GC[0], GC[1])
            epi_dist4 /= 1000   


            # Check great-circle path criteria
            if ( ( abs(baz1 - baz3) < angle_thresh ) and ( abs(baz2 - baz3) < angle_thresh ) ):

                count = count + 1
                print( str(count) + ' : Event ' + str(iev) + ' out of ' + str(len(cat)) )

                # Add this event to subset catalog
                sub_cat.append(ev)

                # Different ways of difining download time window
                # Currently set as event onset time to arrival time of minimum velocity (+ 2000 km of extra propagation)
                # tw1 = time + min(epi_dist1,epi_dist2)/vr_max
                tw1 = time 
                # tw2 = time + max(epi_dist1,epi_dist2)/vr_min
                tw2 = time + (epi_dist4 + 2000.)/vr_min

                # Attampt data download for each station in station pair 
                for istation in [istation1,istation2]:

                    # Path for station
                    datapath = 'Data/' + networks[istation] + '_' + stations[istation] 

                    # Check if this data file was already downloaded (e.g., for a different station pair)
                    filename = datapath + '/' + tstamp + '.pkl'
                    if os.path.isfile(filename): 
                        print('File already exists')
                        continue

                    # Check if this data file was already attempted (i.e., not available on IRIS)
                    if (filename in unavailable_files):
                        print('File already checked - unavailable')
                        continue

                    # Data folder for station
                    if not os.path.isdir(datapath): 
                        print('Path to '+datapath+' doesn`t exist - creating it')
                        os.makedirs(datapath)

                    # Try BH and HH channels
                    try:
                        st = client.get_waveforms(networks[istation],stations[istation],\
                          '*',"BHZ",tw1,tw2,attach_response=True)
                        
                    except:
                        try:
                            st = client.get_waveforms(networks[istation],stations[istation],\
                              '*',"HHZ",tw1,tw2,attach_response=True)

                        except:
                            print('Data not available')
                            unavailable_files.append(filename)
                            continue


                    # Pre-processing
                    # You may want to pick a different frequency band
                    pre_filt = (0.005, 0.006, 2.0, 5.0) ################################################################
                    try:
                        st.remove_response(output='DISP', pre_filt=pre_filt)
                    except:
                        print('No valid response file for this day')
                        unavailable_files.append(filename)
                        continue           

                    for tr in st:
                        #tr.resample(1) 
                        # Filter from 20 - 200 seconds (you may want to modify this)
                        tr.filter('bandpass', freqmin=0.005, freqmax=0.05, corners=2, zerophase=True)
                        sampling_rate = int(tr.stats.sampling_rate)
                        # To minimize storage, resample at 1 Hz
                        tr.decimate(factor=sampling_rate, strict_length=False, no_filter=True)

                    try:
                        zz = st.select(channel='BHZ')[0]
                    except:
                        try:
                            zz = st.select(channel='HHZ')[0]
                        except:
                            unavailable_files.append(filename)
                            continue


                    # Could save data to disc at this point...
                    outfile = open(filename,'wb')
                    pickle.dump((ev,zz),outfile)
                    outfile.close()


        # Write to disc after every station pair
        filename = 'unavailable_files.pkl'
        # Could save data to disc at this point...
        outfile = open(filename,'wb')
        pickle.dump((unavailable_files),outfile)
        outfile.close()


        # Write sub catalog to disc
        filename = 'CAT/' +  networks[istation1] + '_' + stations[istation1] + '_' + \
                             networks[istation2] + '_' + stations[istation2] + '.pkl'

        # Could save data to disc at this point...
        outfile = open(filename,'wb')
        pickle.dump((sub_cat),outfile)
        outfile.close()
