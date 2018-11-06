#-*- coding: utf-8 -*-

"""
Created on 2018-10-12

@author: cmcuervol, cmcuervol@gmail.com

Make nc form hdf files of CloudSat
"""
# from pyhdf import SD
from pyhdf.HDF import *
from pyhdf.V   import *
from pyhdf.VS  import *
from pyhdf.SD  import *

# from numpy import *         # hay que importar todo para que funcione pyhdf
import numpy as np
import os, locale, sys
from netCDF4 import Dataset
import datetime as dt

def describevg(refnum):
    # Describe the vgroup with the given refnum.
    # Open vgroup in read mode.
    vg = v.attach(refnum)
    print( "----------------")
    print( "name:", vg._name, "class:",vg._class, "tag,ref:",)
    print( vg._tag, vg._refnum)

    # Show the number of members of each main object type.
    print( "members: ", vg._nmembers,)
    print( "datasets:", vg.nrefs(HC.DFTAG_NDG),)
    print( "vdatas:  ", vg.nrefs(HC.DFTAG_VH),)
    print( "vgroups: ", vg.nrefs(HC.DFTAG_VG))

    # Read the contents of the vgroup.
    members = vg.tagrefs()

    # Display info about each member.
    index = -1
    for tag, ref in members:
        index += 1
        print( "member index", index)
        # Vdata tag
        if tag == HC.DFTAG_VH:
            vd = vs.attach(ref)
            nrecs, intmode, fields, size, name = vd.inquire()
            print( "  vdata:",name, "tag,ref:",tag, ref)
            print("    fields:",fields)
            print("    nrecs:",nrecs)
            vd.detach()

        # SDS tag
        elif tag == HC.DFTAG_NDG:
            sds = sd.select(sd.reftoindex(ref))
            name, rank, dims, type, nattrs = sds.info()
            print ("  dataset:",name, "tag,ref:", tag, ref)
            print ("    dims:",dims)
            print ("    type:",type)
            sds.endaccess()

        # VS tag
        elif tag == HC.DFTAG_VG:
            vg0 = v.attach(ref)
            print( "  vgroup:", vg0._name, "tag,ref:", tag, ref)
            vg0.detach()

        # Unhandled tag
        else:
            print ("unhandled tag,ref",tag,ref)

    # Close vgroup
    vg.detach()
# #
# # Open HDF file in readonly mode.
# # filename = sys.argv[1]
# filename = path_FRLK+NameHDF_FRLK
# hdf = HDF(filename)
#
# # Initialize the SD, V and VS interfaces on the file.
# sd = SD(filename)
# vs = hdf.vstart()
# v  = hdf.vgstart()
#
# # Scan all vgroups in the file.
# ref = -1
# while 1:
#     try:
#         ref = v.getid(ref)
#         print ref
#     except HDF4Error,msg:    # no more vgroup
#         break
#     describevg(ref)
#


def Listador(directorio, inicio=None, final=None):
    """
    Return the elements (files and directories) of any directory,
    optionaly with filter by start o end of the name of the element
    IMPUTS
    directorio : route of the directory
    inicio     : start of the elements
    final      : end of the elements
    OUTPUTS
    lf  : list of elements
    """
    lf = []
    lista = os.listdir(directorio)
    lista.sort()
    if inicio == final == None:
        return lista
    else:
        for i in lista:
            if inicio == None:
                if i.endswith(final):
                    lf.append(i)
            if final == None:
                if i.startswith(inicio):
                    lf.append(i)
        return lf



def Diff(f,x, order=1, tiempo=None):
    """calculate the derived df/dx
    IMPUTS
    f: array to derivated respect to x
    x: array
    order: order of the derived
    tiempo: optional. If x is a time array put anything, then cacule the derivated
            in time, x must be a datetime object, and the derivate calcules the
            difference in seconds
    OUTPUTS
    g: derived of f respect to x
    """

    # df_dx = np.zeros((len(f)-order, ), dtype= float)
    g = f
    for i in range(order):
        dg_dx = np.zeros((len(f)-(i+1), ), dtype= float)
        if tiempo==None:
            for i in range(len(f)-(i+1)):
                dg_dx[i]= (g[i+1]-g[i])/(x[i+1]-x[i])
        else:
            for i in range(len(f)-(i+1)):
                dg_dx[i]= (g[i+1]-g[i])/(x[i+1]-x[i]).seconds

        g = dg_dx

    return g


def Salto(x, delta):
    "Determina cuando hay un cambio en un array"
    a = np.array(map(lambda i: x[i+1]-x[i], range(len(x)-1)))
    try:
        pos = np.where(a!=delta)[0][0]
    except:
        pos = 0
    return pos



def HDFvars(File):
    """
    Extract the variable names for an hdf file
    """
    # hdfFile = SD.SD(File, mode=1)
    hdfFile = SD(File, mode=1)
    dsets = hdfFile.datasets()
    k = []
    for key in dsets.keys():
        k.append(key)
    k.sort()
    hdfFile.end() # close the file
    return k



def DesHDF(File, name_var, forma='_A'):
    """
    Extrae los datos de HDF para un archivo(file), y para una varaiable (name_var),
    puede ser Ascendente ('_A') o descendente ('_D')
    """
    # hdfFile = SD.SD(File, mode=1)
    hdfFile = SD(File, mode=1)
    # print 'Reading ', File
    d1 = hdfFile.select(name_var+forma)
    var = d1.get()
    print 'Extracting ', name_var
    d1.endaccess()
    hdfFile.end()

    return var


def HDFread(filename, variable, Class=None):
    """
    Extract the data for non scientific data in V mode of hdf file
    """
    hdf = HDF(filename, HC.READ)

    # Initialize the SD, V and VS interfaces on the file.
    sd = SD(filename)
    vs = hdf.vstart()
    v  = hdf.vgstart()

    # Encontrar el puto id de las Geolocation Fields
    if Class == None:
        ref = v.findclass('SWATH Vgroup')
    else:
        ref = v.findclass(Class)

    # Open all data of the class
    vg = v.attach(ref)
    # All fields in the class
    members = vg.tagrefs()

    nrecs = []
    names = []
    for tag, ref in members:
        # Vdata tag
        vd = vs.attach(ref)
        # nrecs, intmode, fields, size, name = vd.inquire()
        nrecs.append(vd.inquire()[0]) # number of records of the Vdata
        names.append(vd.inquire()[-1])# name of the Vdata
        vd.detach()

    idx = names.index(variable)
    var = vs.attach(members[idx][1])
    V   = var.read(nrecs[idx])
    var.detach()
    # Terminate V, VS and SD interfaces.
    v.end()
    vs.end()
    sd.end()
    # Close HDF file.
    hdf.close()

    return np.array(V).ravel()


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
        files = files[0]
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
