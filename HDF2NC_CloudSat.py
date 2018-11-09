#-*- coding: utf-8 -*-

"""
Created on 2018-10-12

@author: cmcuervol, cmcuervol@gmail.com

Make nc form hdf files of CloudSat
"""
from Gadgets.Gadgets import * # Funciones propias
import numpy as np

path_fig  = '/home/cmcuervol/A-Train/CloudSat/'
# path = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-PRECIP-COLUMN.P1_R05/2006/185/'
# path_Data = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1B-CPR.P_R05/'
# path_Geom = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2B-GEOPROF.P1_R05/'
# path_Rain = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-RAIN-PROFILE.P1_R05/'
# # path_PptC = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-PRECIP-COLUMN.P1_R05/2006/153/'
# path_FRLK = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1b-cpr-fl.p_r04/2012/037/'

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

years = Listador(path_Data, inicio='2')

years = [years[0]]
for year in years:
    days = Listador(path_Data+year)
    days = [days[0]]
    for day in days:
        files = Listador(path_Data+year+'/'+day, final='.hdf')
        # files = files[0]
        for i in range(len(files)):
            Lat = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Latitude')
            Lon = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Longitude')
            idx = np.where(((Lat>=Tropical[1]) & (Lat<=Tropical[3])) & ((Lon>=Tropical[0]) & (Lon<=Tropical[2])))
            print(idx)

#


for i  in range(len(files_Data)):
# for i  in range(1):
    Lat = HDFread(path_Data+files_Data[i], 'Latitude')
    Lon = HDFread(path_Data+files_Data[i], 'Longitude')





#
year = '2008'
day  = '001'
path_Data = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1B-CPR.P_R05/'
path_Geop = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2B-GEOPROF.P1_R05/'
path_Rain = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-RAIN-PROFILE.P1_R05/'

path_Ccls = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS.P_R04/'
path_LIDR = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS-LIDAR.P_R04/'
path_TRMM = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2D-CLOUDSAT-TRMM.P_R04/'
path_CWCr = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RO.P_R04/'
path_CWCo = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RVOD.P_R04/'
path_FLUX = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-FLXHR.P2_R04/'


# Unzipeador(path_Ccls+year+'/'+day+'/', erase=False)
# Unzipeador(path_LIDR+year+'/'+day+'/', erase=False)
# Unzipeador(path_CWCr+year+'/'+day+'/', erase=False)
# Unzipeador(path_CWCo+year+'/'+day+'/', erase=False)
# Unzipeador(path_TRMM+year+'/'+day+'/', erase=False)
# Unzipeador(path_FLUX+year+'/'+day+'/', erase=False)
#

files_Data = Listador(path_Data+year+'/'+day, final='hdf' )
files_Geop = Listador(path_Geop+year+'/'+day, final='hdf' )
files_Rain = Listador(path_Rain+year+'/'+day, final='hdf' )

files_Ccls = Listador(path_Ccls+year+'/'+day, final='hdf' )
files_LIDR = Listador(path_LIDR+year+'/'+day, final='hdf' )
files_CWCr = Listador(path_CWCr+year+'/'+day, final='hdf' )
files_CWCo = Listador(path_CWCo+year+'/'+day, final='hdf' )
files_TRMM = Listador(path_TRMM+year+'/'+day, final='hdf' )
files_FLUX = Listador(path_FLUX+year+'/'+day, final='hdf' )

Variables_Data = HDFvars(path_Data+year+'/'+day+'/'+files_Data[0])
Variables_Geop = HDFvars(path_Geop+year+'/'+day+'/'+files_Geop[0])
Variables_Rain = HDFvars(path_Rain+year+'/'+day+'/'+files_Rain[0])

Variables_Ccls = HDFvars(path_Ccls+year+'/'+day+'/'+files_Ccls[0])
Variables_LIDR = HDFvars(path_LIDR+year+'/'+day+'/'+files_LIDR[0])
Variables_CWCr = HDFvars(path_CWCr+year+'/'+day+'/'+files_CWCr[0])
Variables_CWCo = HDFvars(path_CWCo+year+'/'+day+'/'+files_CWCo[0])
Variables_TRMM = HDFvars(path_TRMM+year+'/'+day+'/'+files_TRMM[0])
Variables_FLUX = HDFvars(path_FLUX+year+'/'+day+'/'+files_FLUX[0])

for var in Variables_Data:
    DescribeHDFvar(path_Data+year+'/'+day+'/'+files_Data[0], var)
for var in Variables_Geop:
    DescribeHDFvar(path_Geop+year+'/'+day+'/'+files_Geop[0], var)

for var in Variables_Rain:
    DescribeHDFvar(path_Rain+year+'/'+day+'/'+files_Rain[0], var)

for var in Variables_Ccls:
    DescribeHDFvar(path_Ccls+year+'/'+day+'/'+files_Ccls[0], var)
for var in Variables_LIDR:
    DescribeHDFvar(path_LIDR+year+'/'+day+'/'+files_LIDR[0], var)

for var in Variables_CWCr:
    DescribeHDFvar(path_CWCr+year+'/'+day+'/'+files_CWCr[0], var)
for var in Variables_CWCo:
    DescribeHDFvar(path_CWCo+year+'/'+day+'/'+files_CWCo[0], var)

for var in Variables_FLUX:
    DescribeHDFvar(path_FLUX+year+'/'+day+'/'+files_FLUX[0], var)

for var in Variables_TRMM:
    DescribeHDFvar(path_TRMM+year+'/'+day+'/'+files_TRMM[0], var)








files_Data = Listador(path_Data, final='hdf' )
files_Geom = Listador(path_Geom, final='hdf' )
files_Rain = Listador(path_Rain, final='hdf' )
# files_PptC = Listador(path_PptC, final='hdf' )
# files_FRLK = Listador(path_FRLK, final='zip' )
#
# for zip_file in files_FRLK:
#     os.system('unzip '+ path_FRLK+zip_file+ ' -d '+path_FRLK)
#     # pzf = PyZipFile(path_FRLK+zip_file)
#     # pzf.extractall()
# del files_FRLK
# files_FRLK = Listador(path_FRLK, final='hdf')

NameHDF_Data = files_Data[0]
NameHDF_Geom = files_Geom[0]
NameHDF_Rain = files_Rain[0]
# NameHDF_PptC = files_PptC[0]
# NameHDF_FRLK = files_FRLK[-1]
# NameHDF_FRLK = '2012037001654_30725_CS_1B-CPR-FL_GRANULE_P_R04_E05.hdf'


Variables_CPR = HDFvars(path_Data+NameHDF_Data)
Variables_Geo = HDFvars(path_Geom+NameHDF_Geom)
Variables_Ran = HDFvars(path_Rain+NameHDF_Rain)
# Variables_PPT = HDFvars(path_PptC+NameHDF_PptC)
# Variables_Frl = HDFvars(path_FRLK+NameHDF_FRLK)

print Variables_CPR
print Variables_Geo
print Variables_Ran
# print Variables_PPT
# print Variables_Frl




print "Hello world"
