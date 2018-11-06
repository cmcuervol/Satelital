# from pyhdf import SD
from pyhdf.HDF import *
from pyhdf.V   import *
from pyhdf.VS  import *
from pyhdf.SD  import *

# from numpy import *         # hay que importar todo para que funcione pyhdf
import numpy as np
import os, locale, sys, zipfile

from netCDF4 import Dataset
import datetime as dt


# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                                Basic and system functions
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
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


def Unzipeador(path, ext='.zip', erase=True, path_extract=None):
    """
    Extract and remove the zip files of a directory
    IMPUTS
    path : directory with the zip files
    ext  : files extention
    erase: True to erase the zip files after extract
    path_extract : directory to extract the zip files,
                   if is None, the files are extract in the same directory (path)
    """
    files = Listador(path, final=ext)
    for File in files:
        zip_ref = zipfile.ZipFile(path+File)
        if path_extract is None:
            zip_ref.extractall(path)
        else:
            zip_ref.extractall(path_extract)
        print File
        zip_ref.close()
        if erase == True:
            os.remove(path+File)

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                                HDF files management
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

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

def DescribeHDFvar(filename, variable):
    """
    Describes the info and atributes of HDF variable
    IMPUTS
    filename : complete path and filename
    variable : name of the variable to describe, type STR
    """

    file = SD(filename, SDC.READ)
    print('---------- '+variable+ ' ----------')
    sds_obj = file.select(variable)
    Var = sds_obj.get()
    sds_info = sds_obj.info()
    print(Var.shape)
    print( sds_info )
    print( sds_info[0], sds_info[1] )
    print( 'sds attributes' )
    pprint.pprint( sds_obj.attributes() )
    file.end()

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

def HDFread1D(filename, variable):
    """
    Read HDF file in vs in simple mode
    """
    # Read the file
    f  = HDF(filename)
    # Initialize v mode
    vs = f.vstart()
    # extrect data
    var = vs.attach(variable)
    Var = np.array(var[:]).ravel()
    # Close the file
    var.detach()
    vs.end()
    f.close()

    return Var


def GranuleTime(Filename):
    """
    Convert to datetime format the date of granule data for CloudSat
    """
    Year = int(Filename[ :  4])
    JyDy = int(Filename[ 4: 7])
    Hour = int(Filename[ 7: 9])
    Min  = int(Filename[ 9:11])
    Seg  = int(Filename[11:13])
    GranuTime = dt.datetime(Year,1,1,Hour,Min,Seg)+dt.timedelta(days=JyDy-1)
    return GranuTime
