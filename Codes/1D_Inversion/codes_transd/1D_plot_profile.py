
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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

### Specify station pair name ###
# This code only runs one station pair at a time
station_pair = 'AT_CRAG__TA_R33M_phase_'

# Station separation
epi_dist3 = 456.

# Periods analyzed
period_plot_list = [15.,20.,25.,30.,40.,50.,60.,80.,100.,120.]


# Max depth
max_depth = 200.


fig = plt.figure(figsize=(12,8))

# Number of layers histogram
ax = fig.add_subplot(2,2,1)
filename = station_pair + "misfit.out"
data_misfit = np.loadtxt('results_1d/'+ filename) 
plt.hist(data_misfit[:,1],bins=range(0,20))
ax.set_xlabel('Number of layers')
ax.set_ylabel('Probability')
ax.set_xlim([0,20])
ax.set_yticks([])
ax.set_xticks([0,5,10,15,20])



cmap = mpl.cm.get_cmap("viridis").copy()
cmap.set_under(color='white')

pdiff = np.diff(period_plot_list)
period_plot_lim = [12.5] + list(pdiff/2.+np.array(period_plot_list[0:9])) + [130.]
pvel_axis = np.linspace(2.99,5.01,101)


# Predicted dispersion marginal probability
ax = fig.add_subplot(2,2,3)
filename = station_pair + "dispersion_grid.out"
dispersion_grid = np.loadtxt('results_1d/'+ filename) 
ax.pcolor(period_plot_lim, pvel_axis, np.fliplr(np.transpose(dispersion_grid)), vmin=5, cmap=cmap)

period_plot_lim = np.array(period_plot_lim)
ipi = 3    
ax.plot(period_plot_lim,epi_dist3/period_plot_lim,':',color='k',zorder=10)
ax.plot(period_plot_lim,0.5*epi_dist3/period_plot_lim,'--',color='k',zorder=10)

ax.set_xticks([20.,30.,40.,50.,60.,80.,100.,120.])
ax.set_ylabel('Phase velocity (km/s)')
ax.set_xlabel('Period (s)')
ax.set_xlim([12.5,130.])
ax.set_ylim([3.0,5.0])


# Depth VS marginal probability profile
ax = fig.add_subplot(1,2,2)
filename = station_pair + "profile_grid.out"
profile_grid = np.loadtxt('results_1d/'+ filename) 
ax.imshow(profile_grid, extent = [1.5, 5.5, -max_depth, 0.], aspect='auto',vmin=5, cmap=cmap)
ax.set_yticks([-200,-150,-100, -50, 0])
ax.set_yticklabels(['200','150','100', '50', '0'])
ax.set_xlabel('Vs (km/s)')
ax.set_ylabel('Depth (km)')


# Save figure
plt.savefig('CRAG_R33M_example.png',bbox_inches='tight',dpi=300)



























































