#-*- coding: utf-8 -*-

"""
Created on 2018-10-12

@author: cmcuervol, cmcuervol@gmail.com

Make nc form hdf files of CloudSat
"""
from Gadgets.Gadgets import * # Funciones propias
from netCDF4 import Dataset
import numpy as np

path_Data = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1B-CPR.P_R05/'
path_Geop = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2B-GEOPROF.P1_R05/'
path_Rain = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-RAIN-PROFILE.P1_R05/'

path_Ccls = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS.P_R04/'
path_LIDR = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS-LIDAR.P_R04/'
path_TRMM = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2D-CLOUDSAT-TRMM.P_R04/'
path_CWCr = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RO.P_R04/'
path_CWCo = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RVOD.P_R04/'
path_FLUX = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-FLXHR.P2_R04/'



Tropical = [-107,-30,-10,30]
AMVA     = [-75.85,5.9300,-75.07,6.5900]

# years = Listador(path_Data, inicio='2')
#
# years = [years[0]]
# for year in years:
#     days = Listador(path_Data+year)
#     days = [days[0]]
#     for day in days:
#         files = Listador(path_Data+year+'/'+day, final='.hdf')
#         # files = files[0]
#         for i in range(len(files)):
#             Lat = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Latitude')
#             Lon = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Longitude')
#             idx = np.where(((Lat>=Tropical[1]) & (Lat<=Tropical[3])) & ((Lon>=Tropical[0]) & (Lon<=Tropical[2])))
#             print(idx)
#
# #



# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#
year = '2008'
days  = [['001']]
for day in days:
    try:
        files_Data = Listador(path_Data+year+'/'+day, final='hdf' )
        Data = True
    except:
        Data = False
    try:
        files_Geop = Listador(path_Geop+year+'/'+day, final='hdf' )
        Geop = True
    except:
        Geop = False
    try:
        files_Rain = Listador(path_Rain+year+'/'+day, final='hdf' )
        Rain = True
    except:
        Rain = False
    try:
        Unzipeador(path_Ccls+year+'/'+day+'/', erase=False)
        files_Ccls = Listador(path_Ccls+year+'/'+day, final='hdf' )
        CloudClass = True
    except:
        CloudClass = False
    try:
        Unzipeador(path_LIDR+year+'/'+day+'/', erase=False)
        files_LIDR = Listador(path_LIDR+year+'/'+day, final='hdf' )
        LIDAR = True
    except:
        LIDAR = False
    try:
        Unzipeador(path_CWCr+year+'/'+day+'/', erase=False)
        files_CWCr = Listador(path_CWCr+year+'/'+day, final='hdf' )
        CWC_RO = True
    except:
        CWC_RO = False
    try:
        Unzipeador(path_CWCo+year+'/'+day+'/', erase=False)
        files_CWCo = Listador(path_CWCo+year+'/'+day, final='hdf' )
        CWC_RVOD = True
    except:
        CWC_RVOD =False
    # try:
    #     Unzipeador(path_TRMM+year+'/'+day+'/', erase=False)
    #     files_TRMM = Listador(path_TRMM+year+'/'+day, final='hdf' )
    #     TRMM = True
    # except:
    #     TRMM =False
    try:
        Unzipeador(path_FLUX+year+'/'+day+'/', erase=False)
        files_FLUX = Listador(path_FLUX+year+'/'+day, final='hdf' )
        FLUX = True
    except:
        FLUX =False


#
# files_Data = Listador(path_Data+year+'/'+day, final='hdf' )
# files_Geop = Listador(path_Geop+year+'/'+day, final='hdf' )
# files_Rain = Listador(path_Rain+year+'/'+day, final='hdf' )
#
# files_Ccls = Listador(path_Ccls+year+'/'+day, final='hdf' )
# files_LIDR = Listador(path_LIDR+year+'/'+day, final='hdf' )
# files_CWCr = Listador(path_CWCr+year+'/'+day, final='hdf' )
# files_CWCo = Listador(path_CWCo+year+'/'+day, final='hdf' )
# files_TRMM = Listador(path_TRMM+year+'/'+day, final='hdf' )
# files_FLUX = Listador(path_FLUX+year+'/'+day, final='hdf' )
#
# Variables_Data = HDFvars(path_Data+year+'/'+day+'/'+files_Data[0])
# Variables_Geop = HDFvars(path_Geop+year+'/'+day+'/'+files_Geop[0])
# Variables_Rain = HDFvars(path_Rain+year+'/'+day+'/'+files_Rain[0])
#
# Variables_Ccls = HDFvars(path_Ccls+year+'/'+day+'/'+files_Ccls[0])
# Variables_LIDR = HDFvars(path_LIDR+year+'/'+day+'/'+files_LIDR[0])
# Variables_CWCr = HDFvars(path_CWCr+year+'/'+day+'/'+files_CWCr[0])
# Variables_CWCo = HDFvars(path_CWCo+year+'/'+day+'/'+files_CWCo[0])
# Variables_TRMM = HDFvars(path_TRMM+year+'/'+day+'/'+files_TRMM[0])
# Variables_FLUX = HDFvars(path_FLUX+year+'/'+day+'/'+files_FLUX[0])
#
# for var in Variables_Data:
#     DescribeHDFvar(path_Data+year+'/'+day+'/'+files_Data[0], var)
# for var in Variables_Geop:
#     DescribeHDFvar(path_Geop+year+'/'+day+'/'+files_Geop[0], var)
#
# for var in Variables_Rain:
#     DescribeHDFvar(path_Rain+year+'/'+day+'/'+files_Rain[0], var)
#
# for var in Variables_Ccls:
#     DescribeHDFvar(path_Ccls+year+'/'+day+'/'+files_Ccls[0], var)
# for var in Variables_LIDR:
#     DescribeHDFvar(path_LIDR+year+'/'+day+'/'+files_LIDR[0], var)
#
# for var in Variables_CWCr:
#     DescribeHDFvar(path_CWCr+year+'/'+day+'/'+files_CWCr[0], var)
# for var in Variables_CWCo:
#     DescribeHDFvar(path_CWCo+year+'/'+day+'/'+files_CWCo[0], var)
#
# for var in Variables_FLUX:
#     DescribeHDFvar(path_FLUX+year+'/'+day+'/'+files_FLUX[0], var)
#
# for var in Variables_TRMM:
#     DescribeHDFvar(path_TRMM+year+'/'+day+'/'+files_TRMM[0], var)
#


print "Hello world"
