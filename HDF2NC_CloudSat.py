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
from dateutil.relativedelta import relativedelta
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
mxn = 8000 # Number of maximun values for standarized variable
year = '2008'
start = today = dt.datetime(int(year), 1, 1, 0, 0)

for m in range(1,13):
    actual = dt.datetime(int(year), m, 1)
    future = actual+ relativedelta(months=1)
    begining = (actual - start).days +1
    ending   = (future - start).days

    calendar = ["{:03d}".format(item) for item in range(begining,ending+1)] # calendar days for every month
    folder   = Listador(path_Data+year) # day of the entire folder
    days     = list(set(calendar)&set(folder)) # common days calendar and folder
    days.sort()
    # making empty list to decrease RAM use
    tim = []; lat = [];lon = []
    rep = []
    ref = []; msk = []; hgt = []
    mrf = []; clw = []; piw = []; plw = []
    csc = []
    cfr = []; clb = []; clt = []; cph = []; ppf = []; wlt = []
    lwc_r = []; lnc_r = []; ler_r = []; iwc_r = []; inc_r = []; ier_r = []
    lwc_o = []; lnc_o = []; ler_o = []; iwc_o = []; inc_o = []; ier_o = []

    for day in days:
        files = Listador(path_Data+year+'/'+day, final='.hdf')
        # files = [files[1]]
        for i in range(len(files)): # Always exist
            print( '******************************************')
            print(GranuleTime(files[i]).strftime('%Y-%m-%d %H:%M:%S'))
            # Georreference
            Lat = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Latitude')
            Lon = HDFread(path_Data+year+'/'+day+'/'+files[i], 'Longitude')
            # No include non dataset values in the region
            if VarSplitFilled(Lat, Lat, Lon, Tropical, (mxn, ), Lat.dtype, NoValue=-9999).max()==-9999:
                continue

            tim.append(GranuleTime(files[i]))
            lat.append(VarSplitFilled(Lat, Lat, Lon, Tropical, (mxn, ), Lat.dtype, NoValue=-9999))
            lon.append(VarSplitFilled(Lon, Lat, Lon, Tropical, (mxn, ), Lon.dtype, NoValue=-9999))
            # ======================================================================
            # CPR variables
            REP = DesHDF(path_Data+year+'/'+day+'/'+files[i], 'ReceivedEchoPowers')
            rep.append(VarSplitFilled(REP,Lat, Lon, Tropical, (mxn,REP.shape[1]), REP.dtype, NoValue=-9999))
            # ======================================================================
            # GEOPROF variables
            try:
                f_Geop = Listador(path_Geop+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

                REF = DesHDF(path_Geop+year+'/'+day+'/'+f_Geop, 'Radar_Reflectivity')
                MSK = DesHDF(path_Geop+year+'/'+day+'/'+f_Geop, 'CPR_Cloud_mask')
                HGT = DesHDF(path_Geop+year+'/'+day+'/'+f_Geop, 'Height')
                # Change dtypes for generic
                ref.append(VarSplitFilled(REF,Lat, Lon, Tropical,(mxn,REF.shape[1]), REF.dtype, NoValue=-8888))
                msk.append(VarSplitFilled(MSK,Lat, Lon, Tropical,(mxn,MSK.shape[1]), MSK.dtype, NoValue=-9))
                hgt.append(VarSplitFilled(HGT,Lat, Lon, Tropical,(mxn,HGT.shape[1]), HGT.dtype, NoValue=-9999))
            except:
                ref.append(np.ones((mxn,125), 'int16')*-8888)
                msk.append(np.ones((mxn,125), 'int8' )*-9   )
                hgt.append(np.ones((mxn,125), 'int16')*-9999)
            # ======================================================================
            # cloud clasification variables
            try:
                f_Ccls = Listador(path_Ccls+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

                CSC = DesHDF(path_Ccls+year+'/'+day+'/'+f_Ccls, 'cloud_scenario')
                # Change dtypes for generic
                csc.append(VarSplitFilled(CSC,Lat, Lon, Tropical,(mxn,CSC.shape[1]), CSC.dtype, NoValue=-9999))
            except:
                csc.append(np.ones((mxn,125), 'int16')*-9999)

            # ======================================================================
            # Rain profile variables
            try:
                f_rain = Listador(path_Rain+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

                MRF = DesHDF(path_Rain+year+'/'+day+'/'+f_rain, 'modeled_reflectivity')
                CLW = DesHDF(path_Rain+year+'/'+day+'/'+f_rain, 'cloud_liquid_water')
                PIW = DesHDF(path_Rain+year+'/'+day+'/'+f_rain, 'precip_ice_water')
                PLW = DesHDF(path_Rain+year+'/'+day+'/'+f_rain, 'precip_liquid_water')
                # Change dtypes for generic
                mrf.append(VarSplitFilled(MRF,Lat, Lon, Tropical,(mxn,MRF.shape[1]), MRF.dtype, NoValue=-9999))
                clw.append(VarSplitFilled(CLW,Lat, Lon, Tropical,(mxn,CLW.shape[1]), CLW.dtype, NoValue=-9999))
                piw.append(VarSplitFilled(PIW,Lat, Lon, Tropical,(mxn,PIW.shape[1]), PIW.dtype, NoValue=-9999))
                plw.append(VarSplitFilled(PLW,Lat, Lon, Tropical,(mxn,PLW.shape[1]), PLW.dtype, NoValue=-9999))
            except:
                mrf.append(np.ones((mxn,125), 'int16')*-9999)
                clw.append(np.ones((mxn,125), 'int16')*-9999)
                piw.append(np.ones((mxn,125), 'int16')*-9999)
                plw.append(np.ones((mxn,125), 'int16')*-9999)
            # ======================================================================
            # LIDAR clasification variables
            try:
                f_LIDR = Listador(path_LIDR+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

                CFR = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'CloudFraction')
                CLB = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'CloudLayerBase')
                CLT = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'CloudLayerTop')
                # CLt = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'CloudLayerType')
                CPH = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'CloudPhase')
                # LBF = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'LayerBaseFlag')
                # LBT = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'LayerTopFlag')
                # PHL = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'Phase_log')
                PPF = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'PrecipitationFlag')
                WLT = DesHDF(path_LIDR+year+'/'+day+'/'+f_LIDR, 'Water_layer_top')
                # Change dtypes for generic
                cfr.append(VarSplitFilled(CFR,Lat, Lon, Tropical,(mxn,CFR.shape[1]), CFR.dtype, NoValue=-99.))
                clb.append(VarSplitFilled(CLB,Lat, Lon, Tropical,(mxn,CLB.shape[1]), CLB.dtype, NoValue=-99.))
                clt.append(VarSplitFilled(CLT,Lat, Lon, Tropical,(mxn,CLT.shape[1]), CLT.dtype, NoValue=-99.))
                cph.append(VarSplitFilled(CPH,Lat, Lon, Tropical,(mxn,CPH.shape[1]), CPH.dtype, NoValue=0))
                ppf.append(VarSplitFilled(CPH,Lat, Lon, Tropical,(mxn,CPH.shape[1]), CPH.dtype, NoValue=0))
                wlt.append(VarSplitFilled(WLT,Lat, Lon, Tropical,(mxn,WLT.shape[1]), WLT.dtype, NoValue=-9.))
            except:
                cfr.append(np.ones((mxn,10), 'float32')*-99.)
                clb.append(np.ones((mxn,10), 'float32')*-99.)
                clt.append(np.ones((mxn,10), 'float32')*-99.)
                cph.append(np.ones((mxn,10), 'int8')*-99.)
                ppf.append(np.ones((mxn,10), 'int8')*-99.)
                wlt.append(np.ones((mxn,10), 'float32')*-9.)
            # ======================================================================
            # cloud water content radar only
            try:
                f_CWCr = Listador(path_CWCr+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

                LWC_r = DesHDF(path_CWCr+year+'/'+day+'/'+f_CWCr, 'RO_liq_water_content')
                LNC_r = DesHDF(path_CWCr+year+'/'+day+'/'+f_CWCr, 'RO_liq_number_conc')
                LER_r = DesHDF(path_CWCr+year+'/'+day+'/'+f_CWCr, 'RO_liq_effective_radius')

                IWC_r = DesHDF(path_CWCr+year+'/'+day+'/'+f_CWCr, 'RO_ice_water_content')
                INC_r = DesHDF(path_CWCr+year+'/'+day+'/'+f_CWCr, 'RO_ice_number_conc')
                IER_r = DesHDF(path_CWCr+year+'/'+day+'/'+f_CWCr, 'RO_ice_effective_radius')

                # Change dtypes for generic
                lwc_r.append(VarSplitFilled(LWC_r,Lat, Lon, Tropical,(mxn,LWC_r.shape[1]), LWC_r.dtype, NoValue=-4444))
                lnc_r.append(VarSplitFilled(LNC_r,Lat, Lon, Tropical,(mxn,LNC_r.shape[1]), LNC_r.dtype, NoValue=-4444))
                ler_r.append(VarSplitFilled(LER_r,Lat, Lon, Tropical,(mxn,LER_r.shape[1]), LER_r.dtype, NoValue=-4444))

                iwc_r.append(VarSplitFilled(IWC_r,Lat, Lon, Tropical,(mxn,IWC_r.shape[1]), IWC_r.dtype, NoValue=-4444))
                inc_r.append(VarSplitFilled(INC_r,Lat, Lon, Tropical,(mxn,INC_r.shape[1]), INC_r.dtype, NoValue=-4444))
                ier_r.append(VarSplitFilled(IER_r,Lat, Lon, Tropical,(mxn,IER_r.shape[1]), IER_r.dtype, NoValue=-4444))
            except:
                lwc_r.append(np.ones((mxn,125), 'int16')*-4444)
                lnc_r.append(np.ones((mxn,125), 'int16')*-4444)
                ler_r.append(np.ones((mxn,125), 'int16')*-4444)

                iwc_r.append(np.ones((mxn,125), 'int16')*-4444)
                inc_r.append(np.ones((mxn,125), 'int16')*-4444)
                ier_r.append(np.ones((mxn,125), 'int16')*-4444)
            # ======================================================================
            # cloud water content radar versus optical depth
            try:
                f_CWCo = Listador(path_CWCo+year+'/'+day, inicio=files[i][:13], final='.hdf')[0]

                LWC_o = DesHDF(path_CWCo+year+'/'+day+'/'+f_CWCo, 'RVOD_liq_water_content')
                LNC_o = DesHDF(path_CWCo+year+'/'+day+'/'+f_CWCo, 'RVOD_liq_number_conc')
                LER_o = DesHDF(path_CWCo+year+'/'+day+'/'+f_CWCo, 'RVOD_liq_effective_radius')

                IWC_o = DesHDF(path_CWCo+year+'/'+day+'/'+f_CWCo, 'RVOD_ice_water_content')
                INC_o = DesHDF(path_CWCo+year+'/'+day+'/'+f_CWCo, 'RVOD_ice_number_conc')
                IER_o = DesHDF(path_CWCo+year+'/'+day+'/'+f_CWCo, 'RVOD_ice_effective_radius')

                # Change dtypes for generic
                lwc_o.append(VarSplitFilled(LWC_o,Lat, Lon, Tropical,(mxn,LWC_o.shape[1]), LWC_o.dtype, NoValue=-4444))
                lnc_o.append(VarSplitFilled(LNC_o,Lat, Lon, Tropical,(mxn,LNC_o.shape[1]), LNC_o.dtype, NoValue=-4444))
                ler_o.append(VarSplitFilled(LER_o,Lat, Lon, Tropical,(mxn,LER_o.shape[1]), LER_o.dtype, NoValue=-4444))

                iwc_o.append(VarSplitFilled(IWC_o,Lat, Lon, Tropical,(mxn,IWC_o.shape[1]), IWC_o.dtype, NoValue=-4444))
                inc_o.append(VarSplitFilled(INC_o,Lat, Lon, Tropical,(mxn,INC_o.shape[1]), INC_o.dtype, NoValue=-4444))
                ier_o.append(VarSplitFilled(IER_o,Lat, Lon, Tropical,(mxn,IER_o.shape[1]), IER_o.dtype, NoValue=-4444))
            except:
                lwc_o.append(np.ones((mxn,125), 'int16')*-4444)
                lnc_o.append(np.ones((mxn,125), 'int16')*-4444)
                ler_o.append(np.ones((mxn,125), 'int16')*-4444)

                iwc_o.append(np.ones((mxn,125), 'int16')*-4444)
                inc_o.append(np.ones((mxn,125), 'int16')*-4444)
                ier_o.append(np.ones((mxn,125), 'int16')*-4444)

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

    lat = np.array(lat)
    lon = np.array(lon)

    rep = list2array(rep)

    ref = list2array(ref)
    msk = list2array(msk)
    hgt = list2array(hgt)

    mrf = list2array(mrf)
    clw = list2array(clw)
    piw = list2array(piw)
    plw = list2array(plw)

    csc = list2array(csc)

    cfr = list2array(cfr)
    clb = list2array(clb)
    clt = list2array(clt)
    cph = list2array(cph)
    ppf = list2array(ppf)
    wlt = list2array(wlt)

    lwc_r = list2array(lwc_r)
    lnc_r = list2array(lnc_r)
    ler_r = list2array(ler_r)
    iwc_r = list2array(iwc_r)
    inc_r = list2array(inc_r)
    ier_r = list2array(ier_r)

    lwc_o = list2array(lwc_o)
    lnc_o = list2array(lnc_o)
    ler_o = list2array(ler_o)
    iwc_o = list2array(iwc_o)
    inc_o = list2array(inc_o)
    ier_o = list2array(ier_o)



    cdftime = utime('hours since 1800-01-01 00:00:0.0')
    date    = [cdftime.date2num(i) for i in tim]

    # =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
    # =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
                                            # Crear netCDF
    # =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
    # =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

    path_nc = '/home/cmcuervol/A-Train/CloudSat/'
    Name_nc = 'CloudSat-'+year+"-{:02}".format(m)
    # Variables de los datos originales
    ngeo    = mxn
    ntime   = lat.shape[0]
    nheight = 125
    nclass  = 10


    # Crear el nuevo archivo nc
    nw = Dataset(path_nc+Name_nc+'.nc','w',format='NETCDF4')

    # Definir dimensiones
    ncdim_geo    = nw.createDimension('ngeo',    ngeo)
    ncdim_time   = nw.createDimension('ntime',   ntime)
    ncdim_height = nw.createDimension('nheight', nheight)
    ncdim_class  = nw.createDimension('nclass',  nclass)


    # Crear variables
    ncvar_lat  = nw.createVariable('lat', 'float64',('ntime','ngeo'), zlib=True, complevel=9)
    ncvar_lon  = nw.createVariable('lon', 'float64',('ntime','ngeo'), zlib=True, complevel=9)
    ncvar_time = nw.createVariable('time','float64',('ntime',),       zlib=True, complevel=9)

    ncvar_REP = nw.createVariable('ReceivedEchoPowers', 'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_REF = nw.createVariable('Radar_Reflectivity', 'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_MSK = nw.createVariable('CPR_Cloud_mask',     'int8', ('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_HGT = nw.createVariable('Height',             'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)

    ncvar_MRF = nw.createVariable('modeled_reflectivity','int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_CLW = nw.createVariable('cloud_liquid_water',  'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_PIW = nw.createVariable('precip_ice_water',    'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_PLW = nw.createVariable('precip_liquid_water', 'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)

    ncvar_CSC = nw.createVariable('cloud_scenario',      'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)

    ncvar_CFR = nw.createVariable('CloudFraction',    'float32',('ntime','nclass','ngeo'),zlib=True, complevel=9)
    ncvar_CLB = nw.createVariable('CloudLayerBase',   'float32',('ntime','nclass','ngeo'),zlib=True, complevel=9)
    ncvar_CLT = nw.createVariable('CloudLayerTop',    'float32',('ntime','nclass','ngeo'),zlib=True, complevel=9)
    ncvar_CPH = nw.createVariable('CloudPhase',       'int8',   ('ntime','nclass','ngeo'),zlib=True, complevel=9)
    ncvar_PPF = nw.createVariable('PrecipitationFlag','int8',   ('ntime','nclass','ngeo'),zlib=True, complevel=9)
    ncvar_WLT = nw.createVariable('Water_layer_top',  'float32',('ntime','nclass','ngeo'),zlib=True, complevel=9)

    ncvar_LWC_r = nw.createVariable('RO_liq_water_content',   'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_LNC_r = nw.createVariable('RO_liq_number_conc',     'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_LER_r = nw.createVariable('RO_liq_effective_radius','int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_IWC_r = nw.createVariable('RO_ice_water_content',   'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_INC_r = nw.createVariable('RO_ice_number_conc',     'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_IER_r = nw.createVariable('RO_ice_effective_radius','int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)

    ncvar_LWC_o = nw.createVariable('RVOD_liq_water_content',   'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_LNC_o = nw.createVariable('RVOD_liq_number_conc',     'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_LER_o = nw.createVariable('RVOD_liq_effective_radius','int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_IWC_o = nw.createVariable('RVOD_ice_water_content',   'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_INC_o = nw.createVariable('RVOD_ice_number_conc',     'int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)
    ncvar_IER_o = nw.createVariable('RVOD_ice_effective_radius','int16',('ntime','nheight','ngeo'),zlib=True, complevel=9)


    print 'netCDF variables created'


    # Agregar unidades a las variables
    ncvar_lat .units = 'Degrees north'
    ncvar_lon .units = 'degrees east'
    ncvar_time.units = 'Hours since 1800-01-01'


    ncvar_REP.units = 'unitless'
    ncvar_REF.units = 'dBZ'
    ncvar_MSK.units = 'unitless'
    ncvar_HGT.units = 'm'
    ncvar_MRF.units = 'dBZ'
    ncvar_CLW.units = 'mm'
    ncvar_PIW.units = 'mm'
    ncvar_PLW.units = 'mm'
    ncvar_CFR.units = 'unitless'
    ncvar_CLB.units = 'km'
    ncvar_CLT.units = 'km'
    ncvar_CPH.units = 'unitless'
    ncvar_PPF.units = 'unitless'
    ncvar_WLT.units = 'km'

    ncvar_LWC_r.units = 'mg m^{-3}'
    ncvar_LNC_r.units = 'cm^{-3}'
    ncvar_LER_r.units = 'um'
    ncvar_IWC_r.units = 'mg m^{-3}'
    ncvar_INC_r.units = 'cm^{-3}'
    ncvar_IER_r.units = 'um'

    ncvar_LWC_o.units = 'mg m^{-3}'
    ncvar_LNC_o.units = 'cm^{-3}'
    ncvar_LER_o.units = 'um'
    ncvar_IWC_o.units = 'mg m^{-3}'
    ncvar_INC_o.units = 'cm^{-3}'
    ncvar_IER_o.units = 'um'


    # Agregar nombres largos, a prueba de bobos
    ncvar_lat .longname = 'Array of latitude values'
    ncvar_lon .longname = 'Array of longitude values'
    ncvar_time.longname = 'Hours since 1800-01-01'
    ncvar_REP. longname = 'Received echo powers of Cloud Profiling Radar (CPR).'
    ncvar_REF. longname = 'Radar reflectivity'
    ncvar_MSK. longname = 'Cloud mask of CPR'
    ncvar_HGT. longname = 'Height'
    ncvar_MRF. longname = 'Modeled reflectivity'
    ncvar_CLW. longname = 'Cloud liquid water'
    ncvar_PIW. longname = 'Precipitable ice water'
    ncvar_PLW. longname = 'Precipitable liquid water'
    ncvar_CSC. longname = 'Cloud scenario'
    ncvar_CFR. longname = 'Cloud fraction'
    ncvar_CLB. longname = 'Cloud layer base'
    ncvar_CLT. longname = 'Cloud layer top'
    ncvar_CPH. longname = 'Cloud phase'
    ncvar_PPF. longname = 'Precipitation flag'
    ncvar_WLT. longname = 'Water layer top'

    ncvar_LWC_r.longname = 'Radar-only Liquid Water Content'
    ncvar_LNC_r.longname = 'Radar-only Liquid Number Concentration'
    ncvar_LER_r.longname = 'Radar-only Liquid Effective Radius'
    ncvar_IWC_r.longname = 'Radar-only Ice Water Content'
    ncvar_INC_r.longname = 'Radar-only Ice Number Concentration'
    ncvar_IER_r.longname = 'Radar-only Ice Effective Radius'

    ncvar_LWC_o.longname = 'Radar vs optical depth Liquid Water Content'
    ncvar_LNC_o.longname = 'Radar vs optical depth Liquid Number Concentration'
    ncvar_LER_o.longname = 'Radar vs optical depth Liquid Effective Radius'
    ncvar_IWC_o.longname = 'Radar vs optical depth Ice Water Content'
    ncvar_INC_o.longname = 'Radar vs optical depth Ice Number Concentration'
    ncvar_IER_o.longname = 'Radar vs optical depth Ice Effective Radius'


    nw.title = 'CloudSat of '+year+"-{:02}".format(m)
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
    print( '******************************************')
    print( '    writing variables in netCDF file '     )
    print( '******************************************')
    ncvar_time[:] = date
    ncvar_lat [:,:] = lat
    ncvar_lon [:,:] = lon

    ncvar_REP [:,:,:] = rep

    ncvar_REF [:,:,:] = ref
    ncvar_MSK [:,:,:] = msk
    ncvar_HGT [:,:,:] = hgt

    ncvar_MRF [:,:,:] = mrf
    ncvar_CLW [:,:,:] = clw
    ncvar_PIW [:,:,:] = piw
    ncvar_PLW [:,:,:] = plw

    ncvar_CSC [:,:,:] = csc

    ncvar_CFR [:,:,:] = cfr
    ncvar_CLB [:,:,:] = clb
    ncvar_CLT [:,:,:] = clt
    ncvar_CPH [:,:,:] = cph
    ncvar_PPF [:,:,:] = ppf
    ncvar_WLT [:,:,:] = wlt

    ncvar_LWC_r [:,:,:] = lwc_r
    ncvar_LNC_r [:,:,:] = lnc_r
    ncvar_LER_r [:,:,:] = ler_r
    ncvar_IWC_r [:,:,:] = iwc_r
    ncvar_INC_r [:,:,:] = inc_r
    ncvar_IER_r [:,:,:] = ier_r

    ncvar_LWC_o [:,:,:] = lwc_o
    ncvar_LNC_o [:,:,:] = lnc_o
    ncvar_LER_o [:,:,:] = ler_o
    ncvar_IWC_o [:,:,:] = iwc_o
    ncvar_INC_o [:,:,:] = inc_o
    ncvar_IER_o [:,:,:] = ier_o

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
