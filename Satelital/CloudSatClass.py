#!/usr/bin/env python

from pyhdf.SD import SD, SDC
from pyhdf.HDF import *
from pyhdf.VS import *

import pprint
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm

#----------------------------------------------------------------------------------------#
path_LIDR = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS-LIDAR.P_R04/'
file_path = ''
file_name = '2008183012329_11573_CS_2B-CLDCLASS-LIDAR_GRANULE_P_R04_E02.hdf'

#----------------------------------------------------------------------------------------#
# Read HDF Files (VD data) Latitude & Longitude

f = HDF(file_path+file_name, SDC.READ)
vs = f.vstart()

Latitude = vs.attach('Latitude')
Longitude = vs.attach('Longitude')

a = Latitude[:]

Latitude.detach()
Longitude.detach()
vs.end()
f.close()

#----------------------------------------------------------------------------------------#
# SDS Data

file = SD(file_path+file_name, SDC.READ)

file_info = file.info()
print(file_info)  # number of sds and metadata

print('---------- CloudLayerBase ----------')

sds_obj = file.select('CloudLayerBase') # select sds

CloudLayerBase = sds_obj.get()

sds_info = sds_obj.info()

print(CloudLayerBase.shape)
print( sds_info )
print( sds_info[0], sds_info[1] )
print( 'sds attributes' )
pprint.pprint( sds_obj.attributes() )

for i in range(10):
    print([CloudLayerBase[i,j] for j in range(10)])

print('---------- CloudLayerTop ----------')

sds_obj = file.select('CloudLayerTop') # select sds

CloudLayerTop = sds_obj.get()

sds_info = sds_obj.info()

print(CloudLayerTop.shape)
print( sds_info )
print( sds_info[0], sds_info[1] )
print( 'sds attributes' )
pprint.pprint( sds_obj.attributes() )

for i in range(10):
    print([CloudLayerTop[i,j] for j in range(10)])

print('---------- CloudPhase ----------')

sds_obj = file.select('CloudPhase') # select sds

CloudPhase = sds_obj.get()

sds_info = sds_obj.info()

print(CloudPhase.shape)
print( sds_info )
print( sds_info[0], sds_info[1] )
print( 'sds attributes' )
pprint.pprint( sds_obj.attributes() )

for i in range(10):
    print([CloudPhase[i,j] for j in range(10)])


cldclass_lidar_start_idx = 000
cldclass_lidar_end_idx = 1000

tag = '2'

#----------------------------------------------------------------------------------------#
# Plot cldclass-lidar cloud layers + nb cloud layers (aggregation 2: merge two layers if
# same cloud phase and separation distance < threshold (3 for instance) )

f = plt.figure()
ax = f.add_subplot(111)

#for i in range(CloudLayerBase.shape[0]):
for i in range(cldclass_lidar_start_idx,cldclass_lidar_end_idx):
    nb_cloud_layer = np.where(CloudLayerBase[i,:] < 0 )[0][0]
    #for j in range(nb_cloud_layer):
    #    print(CloudLayerBase[i,j],CloudLayerTop[i,j])
    #print('-------')
    nb_cloud_layer_after_merging = nb_cloud_layer
    for j in range(nb_cloud_layer-1):
            if ( CloudLayerBase[i,j+1] - CloudLayerTop[i,j] ) < 3.0:
                if CloudPhase[i,j+1] == CloudPhase[i,j]:
                    nb_cloud_layer_after_merging = nb_cloud_layer_after_merging - 1
    for j in range(nb_cloud_layer):
        if CloudLayerBase[i,j] > 0 and CloudLayerTop[i,j] > 0.0:
            bar_xi = i
            bar_width = 1.0
            bar_base = CloudLayerBase[i,j]
            bar_height = CloudLayerTop[i,j] - bar_base
            bar_color = '1.0'
            bar_edgecolor = '1.0'
            if nb_cloud_layer_after_merging == 2:
                bar_color = 'blue'
                bar_edgecolor = 'blue'
            if nb_cloud_layer_after_merging > 2:
                bar_color = 'salmon'
                bar_edgecolor = 'salmon'
            plt.bar(left=bar_xi, height=bar_height, width=bar_width, bottom=bar_base, color=bar_color, edgecolor=bar_edgecolor)

ax.set_facecolor('xkcd:lightblue')

plt.bar(left=cldclass_lidar_start_idx, height=0, width=0, bottom=0, color='1.0', edgecolor='1.0', label='nb cloud layer = 1')
plt.bar(left=cldclass_lidar_start_idx, height=0, width=0, bottom=0, color='blue', edgecolor='blue', label='nb cloud layer = 2')
plt.bar(left=cldclass_lidar_start_idx, height=0, width=0, bottom=0, color='salmon', edgecolor='salmon', label='nb cloud layer > 2')

plt.legend()

plt.xlim(cldclass_lidar_start_idx,cldclass_lidar_end_idx)
plt.ylim(0,20)

plt.title('cldclass-lidar')

plt.savefig('plot_'+tag+'_cldclass_lidar_cloud_layers_and_nb_cloud_layers_aggregated.png', dpi=100, bbox_inches='tight' )
#plt.show()

plt.close()

#----------------------------------------------------------------------------------------#
# Plot cldclass-lidar cloud layers + nb cloud phase

f = plt.figure()
ax = f.add_subplot(111)

#for i in range(CloudLayerBase.shape[0]):
for i in range(cldclass_lidar_start_idx,cldclass_lidar_end_idx):
    nb_cloud_layer = np.where(CloudLayerBase[i,:] < 0 )[0][0]
    for j in range(nb_cloud_layer):
        if CloudLayerBase[i,j] > 0 and CloudLayerTop[i,j] > 0.0:
            bar_xi = i
            bar_width = 1.0
            bar_base = CloudLayerBase[i,j]
            bar_height = CloudLayerTop[i,j] - bar_base
            bar_color = '1.0'
            bar_edgecolor = '1.0'
            if CloudPhase[i,j] == 3:
                bar_color = 'blue'
                bar_edgecolor = 'blue'
            if CloudPhase[i,j] == 2:
                bar_color = 'salmon'
                bar_edgecolor = 'salmon'
            plt.bar(left=bar_xi, height=bar_height, width=bar_width, bottom=bar_base, color=bar_color, edgecolor=bar_edgecolor)

ax.set_facecolor('xkcd:lightblue')

plt.bar(left=cldclass_lidar_start_idx, height=0, width=0, bottom=0, color='1.0', edgecolor='1.0', label='ice clouds')
plt.bar(left=cldclass_lidar_start_idx, height=0, width=0, bottom=0, color='salmon', edgecolor='salmon', label='mixed clouds')
plt.bar(left=cldclass_lidar_start_idx, height=0, width=0, bottom=0, color='blue', edgecolor='blue', label='liquid clouds')

plt.legend()

plt.xlim(cldclass_lidar_start_idx,cldclass_lidar_end_idx)
plt.ylim(0,20)

plt.title('cldclass-lidar')

plt.savefig( 'plot_'+tag+'_cldclass_lidar_cloud_layers_and_thermodynamic_phase_method_01.png', dpi=100, bbox_inches='tight' )
#plt.show()

plt.close()
