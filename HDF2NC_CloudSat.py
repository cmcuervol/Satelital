#-*- coding: utf-8 -*-

"""
Created on 2018-10-12

@author: cmcuervol, cmcuervol@gmail.com

Make netCDF from hdf files of CloudSat
"""
from Gadgets.Gadgets import * # Funciones propias
from netCDF4 import Dataset
from netcdftime import utime
import datetime as dt
import numpy as np

def VarSplitFilled(var, lat, lon, bounds, shape, dtype, NoValue = -9999,axis=0):
    """
    split variable inside in a determined boundary and filled with a determined value
    IMPUTS
    var : variable to split
    lat : array of latitudes
    lon : array of longitudes
    bounds : boundaries like {'latmin':-35,'latmax':45,'lonmin':-125,'lonmax':-4}
    shape : shape of normalized varaible
    dtype : data type
    NoValue : No value fill
    axis : axis to split
    RETURNS
    Split : variable inside the bounds with determined shape filled with NoValue
    """
    idx = np.where(((lat>=bounds['latmin']) & (lat<=bounds['latmax'])) & ((lon>=bounds['lonmin']) & (lon<=bounds['lonmax'])))[0]
    # cut = Salto(idx,1)
    Split = np.ones(shape, dtype=dtype)*NoValue
    if axis == 0:
        Split[:idx.shape[0]] = var[idx]
    else:
        Split[:,:idx.shape[0]] = var[:,idx]

    return Split


path_Data = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1B-CPR.P_R05/'
path_Geop = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2B-GEOPROF.P1_R05/'
path_Rain = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-RAIN-PROFILE.P1_R05/'

path_Ccls = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS.P_R04/'
path_LIDR = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS-LIDAR.P_R04/'
path_TRMM = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2D-CLOUDSAT-TRMM.P_R04/'
path_CWCr = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RO.P_R04/'
path_CWCo = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RVOD.P_R04/'
path_FLUX = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-FLXHR.P2_R04/'


# Region = {'latmin':-35,'latmax':45,'lonmin':-125,'lonmax':-4}
Tropical = {'latmin':-30,'latmax':30,'lonmin':-107,'lonmax':-10}
# Tropical = [-107,-30,-10,30]
AMVA     = {'latmin':5.930,'latmax':6.590,'lonmin':-75.850,'lonmax':-75.070}

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
tim = []
lat = []
lon = []
rep = []
ref = []
msk = []
hgt = []


year = '2008'

mxn = 8000 # Number of maximun values for standarized variable
# days  = ['001']
days = Listador(path_Data+year)
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
    # try:
    #     files_Rain = Listador(path_Rain+year+'/'+day, final='hdf' )
    #     Rain = True
    # except:
    #     Rain = False
    # try:
    #     Unzipeador(path_Ccls+year+'/'+day+'/', erase=False)
    #     files_Ccls = Listador(path_Ccls+year+'/'+day, final='hdf' )
    #     CloudClass = True
    # except:
    #     CloudClass = False
    # try:
    #     Unzipeador(path_LIDR+year+'/'+day+'/', erase=False)
    #     files_LIDR = Listador(path_LIDR+year+'/'+day, final='hdf' )
    #     LIDAR = True
    # except:
    #     LIDAR = False
    # try:
    #     Unzipeador(path_CWCr+year+'/'+day+'/', erase=False)
    #     files_CWCr = Listador(path_CWCr+year+'/'+day, final='hdf' )
    #     CWC_RO = True
    # except:
    #     CWC_RO = False
    # try:
    #     Unzipeador(path_CWCo+year+'/'+day+'/', erase=False)
    #     files_CWCo = Listador(path_CWCo+year+'/'+day, final='hdf' )
    #     CWC_RVOD = True
    # except:
    #     CWC_RVOD =False
    # # try:
    # #     Unzipeador(path_TRMM+year+'/'+day+'/', erase=False)
    # #     files_TRMM = Listador(path_TRMM+year+'/'+day, final='hdf' )
    # #     TRMM = True
    # # except:
    # #     TRMM =False
    # try:
    #     Unzipeador(path_FLUX+year+'/'+day+'/', erase=False)
    #     files_FLUX = Listador(path_FLUX+year+'/'+day, final='hdf' )
    #     FLUX = True
    # except:
    #     FLUX =False
    files = Listador(path_Data+year+'/'+day, final='.hdf')
    # files = [files[1]]
    for i in range(len(files)): # Always exist
        tim.append(GranuleTime(files[i]))
        print(tim[i].strftime('%Y-%m-%d %H:%M:%S'))
        # Georreference
        Lat = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Latitude')
        Lon = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Longitude')
        lat.append(VarSplitFilled(Lat, Lat, Lon, Tropical, (mxn, ), Lat.dtype, NoValue=-9999))
        lon.append(VarSplitFilled(Lon, Lat, Lon, Tropical, (mxn, ), Lon.dtype, NoValue=-9999))

        # ======================================================================
        # CPR variables
        REP = DesHDF(path_Data+year+'/'+day+'/'+files[i], 'ReceivedEchoPowers')
        rep.append(VarSplitFilled(REP,Lat, Lon, Tropical, (mxn,REP.shape[1]), REP.dtype, NoValue=-9999))
        # ======================================================================
        # GEOPROF variables
        f_Geop = Listador(path_Geop+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

        REF = DesHDF(path_Geop+year+'/'+day+'/'+f_Geop, 'Radar_Reflectivity')
        MSK = DesHDF(path_Geop+year+'/'+day+'/'+f_Geop, 'CPR_Cloud_mask')
        HGT = DesHDF(path_Geop+year+'/'+day+'/'+f_Geop, 'Height')
        ref.append(VarSplitFilled(REF,Lat, Lon, Tropical,(mxn,REF.shape[1]), REF.dtype, NoValue=-8888))
        msk.append(VarSplitFilled(MSK,Lat, Lon, Tropical,(mxn,MSK.shape[1]), MSK.dtype, NoValue=-9))
        hgt.append(VarSplitFilled(HGT,Lat, Lon, Tropical,(mxn,HGT.shape[1]), HGT.dtype, NoValue=-9999))
        # ======================================================================

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                             Convert list to array
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
def list2array(var):
    """
    convert to array a list and change the second dimension for the last whit np.moveaxis
    IMPUTS
    var : list to convert
    """
    var = np.array(var)
    var = np.moveaxis(var, 1, -1)

    return var

# lat = list2array(lat)
# lon = list2array(lon)
# rep = list2array(rep)
# ref = list2array(ref)
# msk = list2array(msk)
# hgt = list2array(hgt)
lat = np.array(lat)
lon = np.array(lon)
rep = list2array(rep)
ref = list2array(ref)
msk = list2array(msk)
hgt = list2array(hgt)

cdftime = utime('hours since 1800-01-01 00:00:0.0')
date    = [cdftime.date2num(i) for i in tim]

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
                                        # Crear netCDF
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

path_nc = '/home/cmcuervol/A-Train/CloudSat/'
Name_nc = 'Prueba2008'
# Variables de los datos originales
ngeo    = mxn
ntime   = lat.shape[0]
nheight = 125


# Crear el nuevo archivo nc
nw = Dataset(path_nc+Name_nc+'.nc','w',format='NETCDF4')

# Definir dimensiones
ncdim_geo    = nw.createDimension('ngeo',    ngeo)
ncdim_time   = nw.createDimension('ntime',   ntime)
ncdim_height = nw.createDimension('nheight',nheight)


# Crear variables
ncvar_lat  = nw.createVariable('lat', 'float64',('ntime','ngeo'))
ncvar_lon  = nw.createVariable('lon', 'float64',('ntime','ngeo'))
ncvar_time = nw.createVariable('time','float64',('ntime',))

ncvar_REP = nw.createVariable('ReceivedEchoPowers', 'float32',('ntime','nheight','ngeo'))
ncvar_REF = nw.createVariable('Radar_Reflectivity', 'float32',('ntime','nheight','ngeo'))
ncvar_MSK = nw.createVariable('CPR_Cloud_mask',     'float32',('ntime','nheight','ngeo'))
ncvar_HGT = nw.createVariable('Height',             'float32',('ntime','nheight','ngeo'))

print 'netCDF variables created'


# Agregar unidades a las variables
ncvar_lat .units = 'Degrees north'
ncvar_lon .units = 'degrees east'
ncvar_time.units = 'Hours since 1800-01-01'


ncvar_REP.units = 'unitless'
ncvar_REF.units = 'dBZ'
ncvar_MSK.units = 'unitless'
ncvar_HGT.units = 'm'

# Agregar nombres largos, a prueba de bobos
ncvar_lat .longname = 'Array of latitude values'
ncvar_lon .longname = 'Array of longitude values'
ncvar_time.longname = 'Hours since 1800-01-01'
ncvar_REP. longname = 'Received echo powers of Cloud Profiling Radar (CPR).'
ncvar_REF. longname = 'Radar reflectivity'
ncvar_MSK. longname = 'Cloud mask of CPR'
ncvar_HGT. longname = 'Height'


nw.title = 'CloudSat'
nw.description = "This ncfile contents the CloudSat data , from "\
                +str(Tropical['latmin']) + ' to '+str(Tropical['latmax']) + \
                'degrees north, and from '+str(Tropical['lonmin'])+' to '+str(Tropical['lonmax']) + \
                " degrees east. The variables in this file are:\
                Geolocalization variables: \
                    StdPressureLev   :  Pressure levels of temperature and trace gas profiles and geopotential height. The array order is from the surface upward, in conformance with WMO standard. Note that the L3 pressure levels are a subset of the 28 L2 pressure levels, restricted to the range of [1.0, 1000.0] hPa.\
                    H2OPressureLev   :  Pressure levels of water vapor levelprofiles\
                    H2OPressureLay   :  Midpoints of pressure layers of water vapor layer profiles. Layer boundaries are at StdPressureLev.\
                    Latitude         :  Array of latitude values at the center of the grid box (Degrees).\
                    Longitude        :  Array of longitude values at the center of the grid box (Degrees).\
                    CoarseCloudLayer :  Midlayer pressures of the 3 coarse cloud layers (hPa)\
                    FineCloudLayer   :  Midlayer pressures of the 12 fine cloud layers (hPa)\
                \
                Physical variables:\
                    Temperature      :  Atmospheric temperature (K)  [StdPressureLev, lat, lon]\
                    SurfAirTemp      :  Temperature of the atmosphere at the Earth's surface. (K)\
                    SurfSkinTemp     :  Surface skin temperature. (K)\
                    TotH2OVap        :  Total integrated column water vapor burden. (kg/m**2) [lat, lon]\
                    H2O_MMR          :  Water vapor mass mixing ratio at standard pressure levels (gm/kg dry air) [H2OPressureLev, lat, lon]\
                    H2O_MMR_Surf     :  Water vapor mass mixing ratio at the surface (gm/kg dry air) [lat, lon]\
                    RelHum           :  Relative humidity over equilibrium phase (Percent) [H2OPressureLev, lat, lon]\
                    RelHumSurf       :  Relative humidity at the surface over equilibrium phase (Percent) [lat, lon]\
                    CloudFrc         :  Combine layer cloud fraction. (unitless) [lat, lon]\
                    FineCloudFrc     :  Cloud fraction at fine cloud resolution. (unitless) [FineCloudLayer,lat, lon]\
                    CoarseCloudFrc   :  Cloud fraction at coarse cloud resolution. 3 layers: low, medium, high (unitless) [CoarseCloudLayer,lat, lon] Midlayer pressures of the 3 coarse cloud layers. Layer boundaries are at {1100., 680., 440., 10.} hPa\
                    OLR              :  Outgoing long-waveradiation flux. (watts/m**2) [lat, lon]\
                "

# nw.spatial_resolution = '1 degree'
# nw.metadatos = 'https://disc.gsfc.nasa.gov/information/documents?title=AIRS%20Documentation'

# Agregar los datos al archivo
print '******************************************'
print '    writing variables in netCDF file '
print '******************************************'
ncvar_time[:] = date
ncvar_lat [:,:] = lat
ncvar_lon [:,:] = lon
ncvar_REP [:,:,:] = rep
ncvar_REF [:,:,:] = ref
ncvar_MSK [:,:,:] = msk
ncvar_HGT [:,:,:] = hgt

# Si no cierra el archivo es como dejar la BD abierta... se lo tira!
nw.close()




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
