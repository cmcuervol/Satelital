#-*- coding: utf-8 -*-

"""
Created by cmcuervol, cmcuervol@gmail.com

Read and process some files of CloudSat
"""
# from pyhdf import SD
from pyhdf.HDF import *
from pyhdf.V   import *
from pyhdf.VS  import *
from pyhdf.SD  import *
import pprint
# from numpy import *         # hay que importar todo para que funcione pyhdf
import numpy as np
import os, locale, sys, zipfile
from netCDF4 import Dataset
import datetime as dt
import matplotlib

matplotlib.use("template")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
from mpl_toolkits.basemap import Basemap

from Gadgets.Gadgets import * # Funciones propias

# locale.setlocale(locale.LC_ALL, ('en_us','utf-8'))

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

Path_fonts = '/home/cmcuervol/Fuentes/' # Atlas
# Path_fuentes = '/home/ccuervo/Fuentes/' # SAL
Path_data = '/home/cmcuervol/A-Train/CALIPSO/opendap.larc.nasa.gov/opendap/CALIPSO/LID_L2_05kmAPro-Standard-V4-20/2018/08/'
Path_fig  = '/home/cmcuervol/A-Train/CALIPSO/Figures/'
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                              Colors and fonts
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*



# Colors for graphics SIATA style
gris70 = (112/255., 111/255., 111/255.)

ColorInfo1  = (82 /255., 183/255.,196/255.)
ColorInfo2  = (55 /255., 123/255.,148/255.)
ColorInfo3  = (43 /255.,  72/255.,105/255.)
ColorInfo4  = (32 /255.,  34/255., 72/255.)
ColorInfo5  = (34 /255.,  71/255., 94/255.)
ColorInfo6  = (31 /255., 115/255.,116/255.)
ColorInfo7  = (39 /255., 165/255.,132/255.)
ColorInfo8  = (139/255., 187/255.,116/255.)
ColorInfo9  = (200/255., 209/255., 93/255.)
ColorInfo10 = (249/255., 230/255., 57/255.)

AzulChimba   = ( 55/255., 132/255., 251/255.)
AzulChimbita = ( 16/255., 108/255., 214/255.)
VerdeChimba  = (  9/255., 210/255.,  97/255.)
Azul         = ( 96/255., 200/255., 247/255.)
Naranja      = (240/255., 108/255.,  34/255.)
RojoChimba   = (240/255.,  84/255., 107/255.)
Verdecillo   = ( 40/255., 225/255., 200/255.)
Azulillo     = ( 55/255., 150/255., 220/255.)
# flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

# Types of fonts Avenir
#
AvenirHeavy = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Heavy.otf')
AvenirBook  = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Book.otf')
AvenirBlack = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Black.otf')
AvenirRoman = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Roman.ttf')



Files = Listador(Path_data, final='.hdf')

V = HDFvars(Path_data+Files[2])
V = ['Aerosol_Layer_Fraction',
 'Aerosol_Multiple_Scattering_Profile_1064',
 'Aerosol_Multiple_Scattering_Profile_532',
 'Atmospheric_Volume_Description',
 'Backscatter_Coefficient_1064',
 'Backscatter_Coefficient_Uncertainty_1064',
 'CAD_Score',
 'Cloud_Layer_Fraction',
 'Column_Feature_Fraction',
 'Column_IAB_Cumulative_Probability',
 'Column_Integrated_Attenuated_Backscatter_532',
 'Column_Optical_Depth_Cloud_532',
 'Column_Optical_Depth_Cloud_Uncertainty_532',
 'Column_Optical_Depth_Stratospheric_Aerosols_1064',
 'Column_Optical_Depth_Stratospheric_Aerosols_532',
 'Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_1064',
 'Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_532',
 'Column_Optical_Depth_Tropospheric_Aerosols_1064',
 'Column_Optical_Depth_Tropospheric_Aerosols_532',
 'Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_1064',
 'Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_532',
 'Day_Night_Flag',
 'Extinction_Coefficient_1064',
 'Extinction_Coefficient_532',
 'Extinction_Coefficient_Uncertainty_532',
 'Extinction_QC_Flag_1064',
 'Extinction_QC_Flag_532',
 'IGBP_Surface_Type',
 'Latitude',
 'Longitude',
 'Minimum_Laser_Energy_532',
 'Molecular_Number_Density',
 'Ozone_Number_Density',
 'Particulate_Depolarization_Ratio_Profile_532',
 'Particulate_Depolarization_Ratio_Uncertainty_532',
 'Perpendicular_Backscatter_Coefficient_532',
 'Perpendicular_Backscatter_Coefficient_Uncertainty_532',
 'Pressure',
 'Profile_ID',
 'Profile_Time',
 'Profile_UTC_Time',
 'Relative_Humidity',
 'Samples_Averaged',
 'Surface_1064_Integrated_Attenuated_Color_Ratio',
 'Surface_1064_Integrated_Depolarization_Ratio',
 'Surface_532_Integrated_Attenuated_Color_Ratio',
 'Surface_532_Integrated_Depolarization_Ratio',
 'Surface_Base_Altitude_1064',
 'Surface_Base_Altitude_532',
 'Surface_Detection_Confidence_1064',
 'Surface_Detection_Confidence_532',
 'Surface_Detection_Flags_1064',
 'Surface_Detection_Flags_532',
 'Surface_Detections_1km_1064',
 'Surface_Detections_1km_532',
 'Surface_Detections_333m_1064',
 'Surface_Detections_333m_532',
 'Surface_Elevation_Statistics',
 'Surface_Integrated_Attenuated_Backscatter_1064',
 'Surface_Integrated_Attenuated_Backscatter_532',
 'Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
 'Surface_Overlying_Integrated_Attenuated_Backscatter_532',
 'Surface_Peak_Signal_1064',
 'Surface_Peak_Signal_532',
 'Surface_Scaled_RMS_Background_1064',
 'Surface_Scaled_RMS_Background_532',
 'Surface_Top_Altitude_1064',
 'Surface_Top_Altitude_532',
 'Surface_Winds',
 'Temperature',
 'Total_Backscatter_Coefficient_532',
 'Total_Backscatter_Coefficient_Uncertainty_532',
 'Tropopause_Height',
 'Tropopause_Temperature']


# DescribeHDFvar(Path_data+Files[2],'Backscatter_Coefficient_1064' )

fecha = dt.datetime.strptime(Files[2].split('.')[1][:-1], '%Y-%m-%dT%H-%M-%SZ')

Lat = DesHDF(Path_data+Files[2], 'Latitude')
Lat = Lat[:,1]
Lon = DesHDF(Path_data+Files[2], 'Longitude')
Lon = Lon[:,1]
Pres = DesHDF(Path_data+Files[2], 'Pressure')
Pres = np.ma.masked_where(Pres==-9999,Pres)
BSC_1064 = DesHDF(Path_data+Files[2], 'Backscatter_Coefficient_1064')
BSC_532 = DesHDF(Path_data+Files[2], 'Total_Backscatter_Coefficient_532')
# BSC_1064 = np.ma.masked_where(BSC_1064==-9999,BSC_1064)
# BSC_532 = np.ma.masked_where(BSC_532==-9999,BSC_532)
BSC_1064 = np.ma.masked_where(BSC_1064==-9999,BSC_1064)
BSC_532 = np.ma.masked_where(BSC_532==-9999,BSC_532)
Temp = DesHDF(Path_data+Files[2], 'Temperature')
Temp = np.ma.masked_where(Temp==-9999,Temp)




def GranulePloter(Data, lats, lons, cmap, norm, cmaplev='default', extend='neither', \
                  scale='linear', labelData='Variable', fecha=dt.datetime(1992,1,15), name='Prueba.png'):
    """
    Plot a CALIPSO data "Granule" is defined as one orbit. A granule starts
    at the first profile that falls on or past the equator on the descending node.

    IMPUTS:
    Data : array 2D of data to plot
    lats : array of latitudes
    lons : array of longitudes
    cmap : color map
    norm : normalization of values for the colormap
    cmaplev : levels of ticks fot the colorbar. If it's 'default', use the default ticks
    extend  : [ 'neither' | 'both' | 'min' | 'max' ] If not 'neither',
             make pointed end(s) for out of range values.
    scale : Type of scale, linear or logarithmic
    lablelData: Name to show in the data plot
    fecha : datetime of the data
    name  : name of the outputfile

    OUTPUTS
    file save in the Path_fig folder
    """

    def OrganizaLon(longitude):
        """
        reacomodate the longitudes values of a satellite track for graph it easily
        """
        n   = np.where(longitude>0)[0] # Indixes with longitudes positives
        a   = Salto(n,1)         # When occurs the change to negative
        idx = n[a]             # Index when is positive again
        copi = longitude.copy()
        copi[idx:] = copi[idx:]-360 # Change to negative te last positive portion of data

        return copi

    #limpiar guevonadas
    plt.cla()
    plt.clf()
    plt.close('all')

    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(2,1,1)

    lons = OrganizaLon(lons)+360
    m = Basemap(ax=ax1,llcrnrlat=-90, llcrnrlon=lons.min(),
                urcrnrlat=90, urcrnrlon=lons.max())
    m.shadedrelief()
    # draw parallels.
    parallels = np.arange(-90.,91.,30.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=7)
    # draw meridians
    meridians = np.arange(-360.,721.,60.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=7)
    ax1.set_title('CloudSat track',fontproperties=AvenirRoman, color=gris70, fontsize=15)
    ax1.set_title(fecha.strftime('%Y-%m-%d %H:%M:%S'), loc='right', fontsize=10, color=gris70)
    x, y = m(lons, lats)
    m.plot(x, y,'--', color=RojoChimba,lw=2)
    # m.scatter(x,y,color=RojoChimba,alpha=0.7)

    ax2 = fig.add_subplot(2,1,2)
    t,h = np.meshgrid(range(Data.shape[0]), range(Data.shape[1]))
    h = h* 30./399  # There are 399 vertical bins, each one approximately 75 m thick
    t = t*5   # The CloudSat data footprint is approximately 5 km along-track.
    #contorno de colores
    if scale == 'linear':
        cs = ax2.contourf(t, h, Data.T[::-1,:], norm=norm, cmap=cmap, levels=np.linspace(norm.vmin, norm.vmax, 256))
    elif scale == 'logarithmic':
        cs = ax2.contourf(t, h, Data.T[::-1,:], norm=norm, cmap=cmap, levels=np.logspace(np.log10(norm.vmin), np.log10(norm.vmax), 256))
    elif scale == 'default': # whitout levels for interpolation
        cs = cs = ax2.contourf(t, h, Data.T[::-1,:], norm=norm, cmap=cmap)
    else:
        print "Parameter scale not valid. This parameter only accept: 'linear', 'logarithmic', or 'default' as values"

    # cs = ax2.pcolormesh(t, h, Data.T[::-1,:], norm=norm, cmap=cmap)
    cbar_ax = fig.add_axes([0.85, 0.1, 0.015, 0.4])
    if cmaplev == 'default':
        cbar= fig.colorbar(cs, cax=cbar_ax, orientation='vertical', extend=extend)
    else:
        cbar= fig.colorbar(cs, cax=cbar_ax, orientation='vertical', extend=extend, ticks=cmaplev)
    cbar.ax.set_ylabel(labelData,fontproperties=AvenirRoman, color=gris70,fontsize=15)

    ax2.set_ylabel('Ray-path[km]',fontproperties=AvenirRoman, color=gris70, fontsize=15)
    ax2.set_xlabel('Distance along-track [km]',fontproperties=AvenirRoman, color=gris70, fontsize=15)

    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.8, top=0.95, wspace=0.2, hspace=0.1)
    plt.savefig(Path_fig + name, transparent=True)

coqueto=newjet()

# # levels_echo = np.array([1E-16, 1E-15,1E-14,1E-13,1E-12,1E-11,1E-10, 1E-9, 1E-8, 1E-7, 1E-6])
# levels_echo = np.logspace(np.log10(1E-16), np.log10(1E-6),11)
# # cmap_echo, norm_echo = cmapeador(stdcolors,levels_echo, 'echo')
# norm_echo = colors.LogNorm(1E-16, 1E-6)
# GranulePloter(REP, Lat, Lon, plt.cm.jet, norm_echo, levels_echo, 'max','logarithmic',\
#               'Received echo powers',fecha, 'ReceivedEchoPowers'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png' )

# levels_BSC_1064 = np.arange(0,1.1,0.25)
levels_BSC_1064 = np.logspace(np.log10(1E-10), np.log10(1E1),11)
norm_BSC_1064   = colors.LogNorm(vmin=1E-10, vmax=10)
GranulePloter(abs(BSC_1064), Lat, Lon, plt.cm.jet, norm_BSC_1064,'default', 'both', 'logarithmic',\
             r'1064 nm Backscatter Coefficient [km$^{-1}$sr$^{-1}$]', fecha, name='Prueba_1064.png')

# levels_BSC_532 = np.arange(0,1.1,0.25)
levels_BSC_532 = np.logspace(np.log10(1E-10), np.log10(1E1),11)
norm_BSC_532   = colors.LogNorm(vmin=1E-10, vmax=10)
GranulePloter(abs(BSC_532), Lat, Lon, plt.cm.jet, norm_BSC_532,'default', 'both', 'logarithmic',\
             r'532 nm Total Backscatter Coefficient [km$^{-1}$sr$^{-1}$]', fecha, name='Prueba_532.png')


levels_Temp = np.arange(-100,101,25)
norm_Temp   = colors.Normalize(vmin=-100, vmax=100)
GranulePloter(Temp, Lat, Lon, coqueto, norm_Temp,levels_Temp, 'both', 'linear',\
             r'Temperature [$^{\circ}$C]', fecha, name='Prueba_Temp.png')

colores = [(201/255.,201/255.,201/255.),
           (249/255.,249/255.,249/255.),
           ( 47/255.,255/255.,255/255.),
           (  0/255., 21/255.,105/255.),
           (  0/255.,106/255., 67/255.),
           (144/255.,255/255.,  0/255.),
           (255/255.,255/255.,  0/255.),
           (255/255.,255/255.,  0/255.),
           (255/255.,204/255.,  0/255.),
           (230/255.,  0/255.,  0/255.),
           (242/255.,114/255.,195/255.),
           (140/255., 13/255.,135/255.),
           (112/255.,111/255.,111/255.),
           (255/255.,255/255.,255/255.)]




print ('Hello world')
