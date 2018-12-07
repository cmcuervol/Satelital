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
import pandas as pd
import os, locale, sys, zipfile
from netCDF4 import Dataset
import datetime as dt
import matplotlib

matplotlib.use("template")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
from matplotlib.ticker import LogFormatterMathtext, LogLocator

from mpl_toolkits.basemap import Basemap

from Gadgets.Gadgets import * # Funciones propias

# locale.setlocale(locale.LC_ALL, ('en_us','utf-8'))

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

Path_fonts = '/home/cmcuervol/Fuentes/' # Atlas
# Path_fuentes = '/home/ccuervo/Fuentes/' # SAL
Path_data = '/home/cmcuervol/A-Train/CALIPSO/opendap.larc.nasa.gov/opendap/CALIPSO/LID_L2_05kmAPro-Standard-V4-20/2018/08/'
Path_fig  = '/home/cmcuervol/A-Train/CALIPSO/Figures/'
Path_dataVal = '/home/cmcuervol/A-Train/CALIPSO/CALIPSO_DATA/'
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

#
# colores = [(201/255.,201/255.,201/255.),
#            (249/255.,249/255.,249/255.),
#            ( 47/255.,255/255.,255/255.),
#            (  0/255., 21/255.,105/255.),
#            (  0/255.,106/255., 67/255.),
#            (144/255.,255/255.,  0/255.),
#            (255/255.,255/255.,  0/255.),
#            (255/255.,255/255.,  0/255.),
#            (255/255.,204/255.,  0/255.),
#            (230/255.,  0/255.,  0/255.),
#            (242/255.,114/255.,195/255.),
#            (140/255., 13/255.,135/255.),
#            (112/255.,111/255.,111/255.),
#            (255/255.,255/255.,255/255.)]

colores = [(  6/255., 48/255.,167/255.),
           ( 22/255.,131/255.,251/255.),
           ( 22/255.,131/255.,251/255.),
           ( 17/255.,127/255.,126/255.),
           (255/255.,253/255., 56/255.),
           (252/255., 13/255., 27/255.),
           (253/255.,129/255.,170/255.),
           ( 70/255., 70/255., 70/255.),
           (180/255.,180/255.,180/255.),
           (255/255.,255/255.,255/255.),
           (255/255.,255/255.,255/255.),
           (255/255.,255/255.,255/255.)]

Melanie = colors.LinearSegmentedColormap.from_list('Melanie',colores)
Melanie.set_under(colores[0])
# flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

# Types of fonts Avenir
#
AvenirHeavy = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Heavy.otf')
AvenirBook  = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Book.otf')
AvenirBlack = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Black.otf')
AvenirRoman = fm.FontProperties(fname=Path_fonts+'AvenirLTStd-Roman.ttf')



Files = Listador(Path_dataVal, final='.hdf')
FilesVal = Listador(Path_dataVal, final='.hdf')

# V = HDFvars(Path_data+Files[-2])
# V = HDFvars(Path_dataVal+Files[1])
V1 = ['Aerosol_Layer_Fraction',
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
V_Val= ['Amplifier_Gain_1064',
 'Attenuated_Backscatter_1064',
 'Calibration_Constant_1064',
 'Calibration_Constant_532',
 'Calibration_Constant_Uncertainty_1064',
 'Calibration_Constant_Uncertainty_532',
 'Day_Night_Flag',
 'Depolarization_Gain_Ratio_532',
 'Depolarization_Gain_Ratio_Uncertainty_532',
 'Earth-Sun_Distance',
 'Frame_Number',
 'IGBP_Surface_Type',
 'Land_Water_Mask',
 'Laser_Energy_1064',
 'Laser_Energy_532',
 'Latitude',
 'Lidar_Mode',
 'Lidar_Submode',
 'Longitude',
 'Molecular_Number_Density',
 'NSIDC_Surface_Type',
 'Noise_Scale_Factor_1064',
 'Noise_Scale_Factor_532_Parallel',
 'Noise_Scale_Factor_532_Perpendicular',
 'Number_Bins_Shift',
 'Off_Nadir_Angle',
 'Ozone_Number_Density',
 'Parallel_Amplifier_Gain_532',
 'Parallel_Background_Monitor_532',
 'Parallel_Column_Reflectance_532',
 'Parallel_Column_Reflectance_Uncertainty_532',
 'Parallel_RMS_Baseline_532',
 'Perpendicular_Amplifier_Gain_532',
 'Perpendicular_Attenuated_Backscatter_532',
 'Perpendicular_Background_Monitor_532',
 'Perpendicular_Column_Reflectance_532',
 'Perpendicular_Column_Reflectance_Uncertainty_532',
 'Perpendicular_RMS_Baseline_532',
 'Pressure',
 'Profile_ID',
 'Profile_Time',
 'Profile_UTC_Time',
 'QC_Flag',
 'QC_Flag_2',
 'RMS_Baseline_1064',
 'Relative_Humidity',
 'Scattering_Angle',
 'Solar_Azimuth_Angle',
 'Solar_Zenith_Angle',
 'Spacecraft_Altitude',
 'Spacecraft_Attitude',
 'Spacecraft_Attitude_Rate',
 'Spacecraft_Position',
 'Spacecraft_Velocity',
 'Subsatellite_Latitude',
 'Subsatellite_Longitude',
 'Subsolar_Latitude',
 'Subsolar_Longitude',
 'Surface_Altitude_Shift',
 'Surface_Elevation',
 'Surface_Wind_Speeds',
 'Temperature',
 'Total_Attenuated_Backscatter_532',
 'Tropopause_Height',
 'Tropopause_Temperature',
 'Viewing_Azimuth_Angle',
 'Viewing_Zenith_Angle']

# DescribeHDFvar(Path_data+Files[-2],'Backscatter_Coefficient_1064' )

# fecha = dt.datetime.strptime(Files[-2].split('.')[1][:-1], '%Y-%m-%dT%H-%M-%SZ')
#
# Lat = DesHDF(Path_data+Files[-2], 'Latitude')
# Lat = Lat[:,1]
# Lon = DesHDF(Path_data+Files[-2], 'Longitude')
# Lon = Lon[:,1]
# Pres = DesHDF(Path_data+Files[-2], 'Pressure')
# Pres = np.ma.masked_where(Pres==-9999,Pres)
# BSC_1064 = DesHDF(Path_data+Files[-2], 'Backscatter_Coefficient_1064')
# BSC_532 = DesHDF(Path_data+Files[-2], 'Total_Backscatter_Coefficient_532')
# # BSC_1064 = np.ma.masked_where(BSC_1064==-9999,BSC_1064)
# # BSC_532 = np.ma.masked_where(BSC_532==-9999,BSC_532)
# BSC_1064 = np.ma.masked_where(BSC_1064<0,BSC_1064)
# BSC_532  = np.ma.masked_where(BSC_532<0,BSC_532)
# Temp = DesHDF(Path_data+Files[-2], 'Temperature')
# Temp = np.ma.masked_where(Temp==-9999,Temp)


plt.rc(    'font',
    size = 40,
    family = fm.FontProperties(
        fname = '{}AvenirLTStd-Book.otf'.format(Path_fonts)
        ).get_name()
)

typColor = '#%02x%02x%02x' % (115,115,115)
plt.rc('axes',labelcolor=typColor, edgecolor=typColor,)#facecolor=typColor)
plt.rc('axes.spines',right=False, top=False, )#left=False, bottom=False)
plt.rc('text',color= typColor)
plt.rc('xtick',color=typColor)
plt.rc('ytick',color=typColor)
# plt.rc('figure.subplot', left=0, right=1, bottom=0, top=1)

H = np.ones((583,),dtype=float)
res= [(5,300),(290,30),(200,60),(55,180),(33,300)]
count = 0
for i,j in res:
    H[count:i+count]*=j
    count+=i
H = H.cumsum()*1E-3 -2

def GranulePloter(Data, H, lats, lons, cmap, norm, cmaplev='default', extend='both', \
                  scale='linear', labelData='Variable', fecha=dt.datetime(1992,1,15), \
                  name='Prueba.png',DEM= None):
    """
    Plot a CALIPSO data "Granule" is defined as one orbit. A granule starts
    at the first profile that falls on or past the equator on the descending node.

    IMPUTS:
    Data : array 2D of data to plot
    H    : Height of data [km]
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
        try:
            n   = np.where(longitude>0)[0] # Indixes with longitudes positives
            a   = Salto(n,1)         # When occurs the change to negative
            idx = n[a]             # Index when is positive again
            copi = longitude.copy()
            copi[idx:] = copi[idx:]-360 # Change to negative te last positive portion of data
            return copi
        except:
            return longitude

    #limpiar guevonadas
    plt.cla()
    plt.clf()
    plt.close('all')

    fig = plt.figure(figsize=(13,10))
    # ax2 = fig.add_subplot(1,1,1)
    ax2 = fig.add_axes([0,0,1,1])
    ax2.set_title(fecha.strftime('%Y-%m-%d %H:%M:%S'), loc='right', fontsize=10, color=gris70)
    # t,h = np.meshgrid(range(Data.shape[0]), range(Data.shape[1]))
    delta = (H[1:]-H[:-1])/2
    delta = np.append(delta, delta[-1])
    H = H-delta
    H = np.append(H, H[-1]+delta[-1])
    t,h = np.meshgrid(range(Data.shape[0]+1), H)
    # h = h* 30./Data.shape[1] -2  # There are Data.shape[1] vertical bins, each one approximately ∆m thick
    # t = t*5   # The CloudSat data footprint is approximately 5 km along-track.
    #contorno de colores
    if scale == 'linear':
        cs = ax2.pcolormesh(t, h, Data.T[::-1,:], norm=norm, cmap=cmap, levels=np.linspace(norm.vmin, norm.vmax, 256))
    elif scale == 'logarithmic':
        # cs = ax2.pcolormesh(t, h, Data.T[::-1,:], norm=norm, cmap=cmap, levels=np.logspace(np.log10(norm.vmin), np.log10(norm.vmax), 256))
        cs = ax2.pcolormesh(t, h, Data.T[::-1,:], norm=norm, cmap=cmap)
    elif scale == 'default': # whitout levels for interpolation
        cs = ax2.pcolormesh(t, h, Data.T[::-1,:], norm=norm, cmap=cmap)
    else:
        print ("Parameter scale not valid. This parameter only accept: 'linear', 'logarithmic', or 'default' as values")

    # cs = ax2.pcolormesh(t, h, Data.T[::-1,:], norm=norm, cmap=cmap)
    cbar_ax = fig.add_axes([1.02, 0.2, 0.015, 0.6])
    if cmaplev == 'default':
        cbar= plt.colorbar(cs, cax=cbar_ax, norm=norm,orientation='vertical', extend=extend)
    elif scale == 'logarithmic':
        minorTicks = np.hstack([np.arange(1,10,1)*log for log in np.logspace(-10,1,12)])
        minorTicks = minorTicks[(minorTicks >= norm.vmin) & (minorTicks <=norm.vmax)]
        cbar = plt.colorbar(cs, cax=cbar_ax, norm =norm, orientation='vertical', extend=extend, format = LogFormatterMathtext(10) ,ticks=LogLocator(10))
        cbar.ax.yaxis.set_ticks(cs.norm(minorTicks), minor=True)
        cbar.ax.tick_params(which='minor',width=1,length=4)
        cbar.ax.tick_params(which='major',width=1,length=6)
    else:
        cbar= plt.colorbar(cs, cax=cbar_ax, norm=norm, orientation='vertical', extend=extend, ticks=cmaplev)
    cbar.ax.set_ylabel(labelData)

    ax2.set_ylabel('Altitude [km]')
    labels = [str(x) for x in ax2.get_yticks()]
    labels[0] = ''
    ax2.set_yticklabels(labels)
    # ax2.set_xlabel('Distance along-track [km]',fontproperties=AvenirRoman, color=gris70, fontsize=15)
    ax2.text(1.05,-0.075,'Lat\nLon', transform=ax2.transAxes)
    ticks = ax2.get_xticks()
    ax2.set_xticklabels(map(lambda x: '%.2f\n%.2f' %(lats[int(x)], lons[int(x)]), ticks[:-1]))
    if DEM is not None:
        ax2.plot(DEM,lw=2, color= 'k')
    ax2.set_ylim(bottom=17,top=23)
    ax2.set_xlim(left=0,right=Data.shape[0]+1)
    # plt.subplots_adjust(left=0.125, bottom=0.1, right=0.8, top=0.95, wspace=0.2, hspace=0.1)
    # ax1 = fig.add_subplot(2,1,1)
    ax1 = fig.add_axes([0.55,0.62,0.4,0.35])

    lons = OrganizaLon(lons)+360
    m = Basemap(ax=ax1,llcrnrlat=lat.min(), llcrnrlon=lons.min()-5,
                urcrnrlat=lat.max(), urcrnrlon=lons.max()+5, resolution='h')

    m.shadedrelief()
    # draw parallels.
    # parallels = np.arange(-90.,91.,30.)
    parallels = np.arange(-90.,91.,5.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=7)
    # draw meridians
    # meridians = np.arange(-360.,721.,60.)
    meridians = np.arange(-360.,721.,5.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=7)
    ax1.set_title('CALIPSO track',fontproperties=AvenirRoman, color=gris70, fontsize=15)
    x, y = m(lons, lats)
    m.plot(x, y,'--', color=RojoChimba,lw=2)
    # m.scatter(x,y,color=RojoChimba,alpha=0.7)

    plt.savefig(Path_fig + name, transparent=True, bbox_inches='tight')

coqueto=newjet()


def running_mean(x, N):
    """
    Movil mean with a window
    IMPUTS
    x : array
    N : window
    """
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def ProfilePlot(Data, H, surface,Data2=None, H2=None,scale='log', labelData='Variable', fecha=dt.datetime(1992,1,15),\
                name='Prueba.png'):
    """
    Plot a CALIPSO single profile of data
    IMPUTS:
    Data : array 1D of data to plot
    H    : Height of data [km]
    Data2: array 1D of data to plot of LIDAR ground based
    H2   : Height of data of LIDAR ground based [km]
    scale : Type of scale ["linear", "log", "symlog", "logit"]
    lablelData: Name to show in the data plot
    fecha : datetime of the data
    name  : name of the outputfile

    OUTPUTS
    file save in the Path_fig folder
    """
    plt.cla()
    plt.clf()
    plt.close('all')

    fig = plt.figure(figsize=(10,16))
    ax1  = fig.add_axes([0,0,1,1])
    idx = np.where(H>=surface)[0]
    ax1.plot(Data[idx],H[idx]-H[idx][0], color=AzulChimba, alpha=0.7)
    # ax.plot(Data,H, color=AzulChimba, alpha=0.7)
    # media movil
    a = pd.DataFrame(Data[idx], index=H[idx]-H[idx][0])
    # a = pd.DataFrame(Data, index=H)
    c = pd.rolling_median(a,window=10,min_periods=1,center=True)
    ax1.plot(c.values.ravel(),c.index.values,linewidth=5,color=Azulillo)
    ax1.set_xscale(scale)
    ax1.set_ylabel('Altitude [km]')
    ax1.set_xlabel(labelData)
    if Data2 is not None:
        ax2 = ax1.twiny()
        Data2[Data2<1E-4] = np.nan
        ax2.plot(Data2,H2, color=Naranja)
        ax2.set_xlabel(r'RCS [mVkm$^{2}$]', color=Naranja)
        ax2.tick_params('y', colors=Naranja)
        ax2.set_xscale(scale)
        ax2.set_xlim(1E-4,16)
    plt.savefig(Path_fig + name, transparent=True, bbox_inches='tight')

def ProfilePlotVerificador(Data, H, surface,Data2=None, H2=None,scale='log', \
                           labelData='Variable', fecha=dt.datetime(1992,1,15),\
                           name='Prueba.png'):
    """
    Plot a CALIPSO single profile of data
    IMPUTS:
    Data : array 1D of data to plot
    H    : Height of data [km]
    Data2: array 1D of data to plot of LIDAR ground based
    H2   : Height of data of LIDAR ground based [km]
    scale : Type of scale ["linear", "log", "symlog", "logit"]
    lablelData: Name to show in the data plot
    fecha : datetime of the data
    name  : name of the outputfile

    OUTPUTS
    file save in the Path_fig folder
    """
    plt.cla()
    plt.clf()
    plt.close('all')

    fig = plt.figure(figsize=(10,16))
    ax1  = fig.add_axes([0,0,1,1])
    idx = np.where(H>=surface)[0]
    ax1.plot(Data[idx],H[idx]-H[idx][0], color=AzulChimba, alpha=0.7)
    # ax.plot(Data,H, color=AzulChimba, alpha=0.7)
    # media movil
    a = pd.DataFrame(Data[idx], index=H[idx]-H[idx][0])
    # a = pd.DataFrame(Data, index=H)
    c = pd.rolling_median(a,window=10,min_periods=1,center=True)
    ax1.plot(c.values.ravel(),c.index.values,linewidth=2,color=Azulillo)
    ax1.set_xscale(scale)
    ax1.set_ylabel('Altitude [km]')
    ax1.set_xlabel(labelData)
    if Data2 is not None:
        ax1.plot(Data2[idx],H2[idx]-H2[idx][0], color=Naranja,alpha=0.7)
        a = pd.DataFrame(Data2[idx], index=H2[idx]-H2[idx][0])
        c = pd.rolling_median(a,window=10,min_periods=1,center=True)
        ax1.plot(c.values.ravel(),c.index.values,linewidth=2,color=Naranja)
    plt.savefig(Path_fig + name, transparent=True, bbox_inches='tight')


def ProfileMean(Data, H, surface,lats, lons,hdown=18, hup=20,scale='log', \
                           labelData='Variable', name='Prueba.png'):
    """
    Plot a CALIPSO single profile of data
    IMPUTS:
    Data : array 2D of data to plot
    H    : Height of data [km]
    lats : latitude array 1D
    lons : longitude array 1D
    scale: Type of scale ["linear", "log", "symlog", "logit"]
    lablelData: Name to show in the data plot
    name  : name of the outputfile

    OUTPUTS
    file save in the Path_fig folder
    """
    plt.cla()
    plt.clf()
    plt.close('all')

    fig = plt.figure(figsize=(10,16))
    ax1  = fig.add_axes([0,0,1,1])
    idx = np.where(H>=surface)[0]
    H   = H[idx]-H[idx][0]
    dat = Data[:,idx]
    hid = np.where((H>hdown)&(H<hup))[0]
    bsc = np.mean(dat[:,hid],axis=1)

    ax1.plot(bsc,lats, color=AzulChimba, alpha=0.7)
    # media movil
    a = pd.DataFrame(bsc, index=lats)
    c = pd.rolling_median(a,window=10,min_periods=1,center=True)
    ax1.plot(c.index.values,c.values.ravel(),linewidth=2,color=Azulillo)
    ax1.set_yscale(scale)
    ax1.set_ylabel(labelData)
    ax1.text(1.05,-0.075,'Lat\nLon', transform=ax1.transAxes)
    ticks = ax1.get_xticks()
    ax1.set_xticklabels(map(lambda x: '%.2f\n%.2f' %(lats[int(x)], lons[int(x)]), ticks[:-1]))
    ax1.set_ylim(bottom=1E-3, top=1E-1)
    ax1.set_xlim(left=lats.min(), right=lats.max())
    plt.savefig(Path_fig + name, transparent=True, bbox_inches='tight')


def VarSplit(var, lat, lon, bounds, axis=0):
    """
    split variable inside in a determined boundary and filled with a determined value
    IMPUTS
    var : variable to split
    lat : array of latitudes
    lon : array of longitudes
    bounds : boundaries like {'latmin':-35,'latmax':45,'lonmin':-125,'lonmax':-4}
    axis   : axis to split
    RETURNS
    Split : variable inside the bounds
    """
    idx = np.where(((lat>=bounds['latmin']) & (lat<=bounds['latmax'])) & ((lon>=bounds['lonmin']) & (lon<=bounds['lonmax'])))[0]
    # cut = Salto(idx,1)
    if axis == 0:
        Split = var[idx]
    else:
        Split = var[:,idx]

    return Split

Tropical = {'latmin':-30,'latmax':30,'lonmin':-107,'lonmax':-10}
AMVA     = {'latmin':5.930,'latmax':6.590,'lonmin':-75.850,'lonmax':-75.070}
Colombia = {'latmin':-5,'latmax':13,'lonmin':-82,'lonmax':-65}


# 2  : 2018-08-17
# 3  : 2018-08-30
# 6  : 2018-11-17
# 7 : 2018-11-24

idxfile = 6

fecha = dt.datetime.strptime(FilesVal[idxfile].split('.')[1][:-1], '%Y-%m-%dT%H-%M-%SZ')

print(fecha.strftime('%Y-%m-%d_%H-%M-%S'))

Lat = DesHDF(Path_dataVal+FilesVal[idxfile], 'Latitude' ).ravel()
Lon = DesHDF(Path_dataVal+FilesVal[idxfile], 'Longitude').ravel()
Pres = DesHDF(Path_dataVal+FilesVal[idxfile], 'Pressure')
Pres = np.ma.masked_where(Pres==-9999,Pres)

BSC_1064 = DesHDF(Path_dataVal+FilesVal[idxfile], 'Attenuated_Backscatter_1064')
BSC_532  = DesHDF(Path_dataVal+FilesVal[idxfile], 'Total_Attenuated_Backscatter_532')
# BSC_1064 = np.ma.masked_where(BSC_1064==-9999,BSC_1064)
# BSC_532  = np.ma.masked_where(BSC_532 ==-9999,BSC_532)
BSC_1064[BSC_1064<1E-4] = 1.0001E-4
BSC_532 [BSC_532 <1E-4] = 1.0001E-4
# BSC_1064[BSC_1064>1E-1] = 1E-1
# BSC_532 [BSC_532 >1E-1] = 1E-1
BSC_1064 = np.ma.masked_where(BSC_1064>1E-1,BSC_1064)
BSC_532  = np.ma.masked_where(BSC_532 >1E-1,BSC_532)
# BSC_1064 = np.ma.masked_where(BSC_1064<0,BSC_1064)
# BSC_532  = np.ma.masked_where(BSC_532 <0,BSC_532)

Temp = DesHDF(Path_dataVal+FilesVal[idxfile], 'Temperature')
Srf = DesHDF(Path_dataVal+FilesVal[idxfile], 'Surface_Elevation').ravel()
# Surf = DesHDF(Path_dataVal+FilesVal[idxfile], 'Surface_Altitude_Shift')
Temp = np.ma.masked_where(Temp==-9999,Temp)



# levels_BSC_1064 = np.arange(0,1.1,0.25)
# levels_BSC_1064 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_BSC_1064   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(BSC_1064, Lat, Lon, plt.cm.jet, norm_BSC_1064,norm_BSC_1064, 'both', 'logarithmic',\
#              r'1064 nm Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, name='Attenuated_Backscatter_1064_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# # levels_BSC_532 = np.arange(0,1.1,0.25)
# levels_BSC_532 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_BSC_532   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(BSC_532, Lat, Lon, plt.cm.jet, norm_BSC_532,levels_BSC_532, 'both', 'logarithmic',\
#              r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, name='Total_Attenuated_Backscatter_532_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
# levels_Temp = np.arange(-100,101,25)
# norm_Temp   = colors.Normalize(vmin=-100, vmax=100)
# GranulePloter(Temp, Lat, Lon, coqueto, norm_Temp,levels_Temp, 'both', 'linear',\
#              r'Temperature [$^{\circ}$C]', fecha, name='New_Temp.png')

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                               Splited variables
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

bsc_1064 = VarSplit(BSC_1064, Lat, Lon, AMVA)
bsc_532  = VarSplit(BSC_532,  Lat, Lon, AMVA)
lat = VarSplit(Lat, Lat, Lon, AMVA)
lon = VarSplit(Lon, Lat, Lon, AMVA)
srf = VarSplit(Srf, Lat, Lon, AMVA)

# levels_bsc_1064 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_1064   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_1064, H, lat, lon, plt.cm.jet, norm_bsc_1064,norm_bsc_1064, 'both', 'logarithmic',\
#              r'1064 nm Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, name='Attenuated_Backscatter_1064_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png', DEM =srf)
#
# # levels_bsc_532 = np.arange(0,1.1,0.25)
# levels_bsc_532 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_532   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_532,H, lat, lon, plt.cm.jet, norm_bsc_532,levels_bsc_532, 'both', 'logarithmic',\
#              r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, name='Total_Attenuated_Backscatter_532_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png',DEM = srf)

# levels_bsc_1064 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_1064   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_1064, H, lat, lon, Melanie, norm_bsc_1064,norm_bsc_1064, 'max', 'logarithmic',\
#              r'1064 nm Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, \
#              name='Attenuated_Backscatter_1064_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'_Col.png', DEM =srf)
#
# # levels_bsc_532 = np.arange(0,1.1,0.25)
# levels_bsc_532 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_532   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_532, H, lat, lon, Melanie, norm_bsc_532,levels_bsc_532, 'max', 'logarithmic',\
#              r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, \
#              name='Total_Attenuated_Backscatter_532_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'_Col.png',DEM = srf)

# levels_bsc_1064 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_1064   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_1064, H, lat, lon, Melanie, norm_bsc_1064,'default', 'max', 'logarithmic',\
#              r'1064 nm Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, \
#              name='Attenuated_Backscatter_1064_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'_AMVA.png', DEM =srf)
#
# # levels_bsc_532 = np.arange(0,1.1,0.25)
# levels_bsc_532 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_532   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_532, H, lat, lon, Melanie, norm_bsc_532,'default', 'max', 'logarithmic',\
#              r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, \
#              name='Total_Attenuated_Backscatter_532_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'_AMVA.png',DEM = srf)


Perfiles = pd.read_csv('/home/cmcuervol/A-Train/CALIPSO/Perfiles.csv', index_col=0)
# hueco = np.where(srf<1.5)[0] # For 2018-11-24
# hueco = range(91,121) # For 2018-08-17  2  1
# hueco = range(75,116) # For 2018-08-30  3 0
hueco = range(33,61) # For 2018-11-17  6  2
Lidar = Perfiles.iloc[:,2].values
# #
ProfilePlot(np.mean(bsc_532[hueco,:],axis=0), H, srf[hueco].min(), Lidar, Perfiles.index.values,\
            labelData=r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]',\
            name='ProfileCALIPSOvsLIDAR'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
ProfilePlot(np.mean(bsc_532[hueco,:],axis=0), H, srf[hueco].min(), \
            labelData=r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]',\
            name='ProfileCALIPSO_CD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')

# ProfilePlotVerificador(np.mean(bsc_532[hueco,:],axis=0), H, srf[hueco].min(), np.max(bsc_532[hueco,:],axis=0), H,\
#             labelData=r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]',\
#             name='ProfileCALIPSO_RevMax'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
# ProfilePlotVerificador(np.mean(bsc_532[hueco,:],axis=0), H, srf[hueco].min(),np.min(bsc_532[hueco,:],axis=0), H, \
#             labelData=r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]',\
#             name='ProfileCALIPSO_RevMin'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
# ProfileMean(bsc_532[hueco,:], H, srf[hueco].min(),lat[hueco],lon[hueco], 18,20, \
#             labelData=r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]',\
#             name='ProfileCALIPSO_Vmean'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
# ProfileMean(bsc_532, H, srf.min(),lat,lon, 18,20, \
#             labelData=r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]',\
#             name='ProfileCALIPSO_Vmean'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')

# levels_bsc_532 = np.arange(0,1.1,0.25)
# levels_bsc_532 = np.logspace(np.log10(1E-4), np.log10(1E-1),11)
# norm_bsc_532   = colors.LogNorm(vmin=1E-4, vmax=1E-1)
# GranulePloter(bsc_532[hueco,:], H, lat[hueco], lon[hueco], Melanie, norm_bsc_532,'default', 'max', 'logarithmic',\
#              r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, \
#              name='Total_Attenuated_Backscatter_532_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'_AMVA_Vprof.png',DEM = srf[hueco])
# GranulePloter(bsc_532*1E1, H, lat, lon, Melanie, norm_bsc_532,'default', 'max', 'logarithmic',\
#              r'532 nm Total Attenuated Backscatter [km$^{-1}$sr$^{-1}$]', fecha, \
#              name='Total_Attenuated_Backscatter_532_'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'_AMVA_Vprof.png',DEM = srf[hueco])

print ('Hello world')
