# from pyhdf import SD
from pyhdf.HDF import *
from pyhdf.V   import *
from pyhdf.VS  import *
from pyhdf.SD  import *

# from numpy import *         # hay que importar todo para que funcione pyhdf
import numpy as np
import os, locale, sys, zipfile, pprint

from netCDF4 import Dataset
import datetime as dt

import matplotlib
matplotlib.use("template")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.font_manager as fm
from mpl_toolkits.basemap import Basemap
# import pkg_pakage Buscar para las rutas de los paquetes

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
            if (inicio is not None) & (final is not None):
                if i.startswith(inicio) & i.endswith(final):
                    lf.append(i)
        return lf

def running_mean(x, N):
    """
    Movil mean with a window
    IMPUTS
    x : array
    N : window
    """
    cumsum = numpy.cumsum(numpy.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

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

def DesHDF(File, name_var, forma=''):
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

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#                                    cmaps
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


def cmapeador(colrs=None, levels=None, name='coqueto'):
    """
    Make a new color map with the colors and levels given, adapted from Julian Sepulveda & Carlos Hoyos

    IMPUTS
    colrs  : list of tuples of RGB colors combinations
    levels : numpy array of levels correspond to each color in colors
    name   : name to register the new cmap

    OUTPUTS
    cmap   : color map
    norm   : normalization with the levesl given optimized for the cmap
    """
    if colrs == None:
        colrs = [(255, 255, 255),(0, 255, 255), (0, 0, 255),(70, 220, 45),(44, 141, 29),\
                  (255,255,75),(255,142,0),(255,0,0),(128,0,128),(102,0,102),\
                  (255, 153, 255)]

    if levels == None:
        levels = np.array([0.,1.,5.,10.,20.,30.,45.,60., 80., 100., 150.])
    # print levels
    scale_factor   = ((255-0.)/(levels.max() - levels.min()))
    new_Limits     = list(np.array(np.round((levels-levels.min()) * scale_factor/255.,3),dtype=float))
    Custom_Color   = map(lambda x: tuple(ti/255. for ti in x) , colrs)
    nueva_tupla    = [((new_Limits[i]),Custom_Color[i],) for i in range(len(Custom_Color))]
    cmap_new       = colors.LinearSegmentedColormap.from_list(name,nueva_tupla)
    levels_nuevos  = np.linspace(np.min(levels),np.max(levels),255)
    # print levels_nuevos
    # print new_Limits
    # levels_nuevos  = np.linspace(np.min(levels),np.max(levels),1000)
    norm_new       = colors.BoundaryNorm(boundaries=levels_nuevos, ncolors=256)
    # norm           = colors.BoundaryNorm(boundaries=levels_nuevos, ncolors=1000)

    return cmap_new, norm_new


def newjet(cmap="jet"):
    """
    function to make a newd colorbar with white at center
    IMPUTS
    cmap: colormap to change
    RETURNS
    newcmap : colormap with white as center
    """
    jetcmap = plt.cm.get_cmap(cmap, 11) #generate a jet map with 11 values
    jet_vals = jetcmap(np.arange(11)) #extract those values as an array
    jet_vals[5] = [1, 1, 1, 1] #change the middle value
    newcmap = colors.LinearSegmentedColormap.from_list("newjet", jet_vals)
    return newcmap


class MidpointNormalize(colors.Normalize):
    """
    New Normalization with a new parameter: midpoint
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


class SqueezedNorm(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, mid=0, s1=1, s2=1, clip=False):
        self.vmin = vmin # minimum value
        self.mid  = mid  # middle value
        self.vmax = vmax # maximum value
        self.s1=s1; self.s2=s2
        f = lambda x, zero,vmax,s: np.abs((x-zero)/(vmax-zero))**(1./s)*0.5
        self.g = lambda x, zero,vmin,vmax, s1,s2: f(x,zero,vmax,s1)*(x>=zero) - \
                                             f(x,zero,vmin,s2)*(x<zero)+0.5
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        r = self.g(value, self.mid,self.vmin,self.vmax, self.s1,self.s2)
        return np.ma.masked_array(r)


# fig, (ax, ax2, ax3) = plt.subplots(nrows=3,
#                                    gridspec_kw={"height_ratios":[3,2,1], "hspace":0.25})
#
# x = np.linspace(-13,4, 110)
# norm=SqueezedNorm(vmin=-13, vmax=4, mid=0, s1=1.7, s2=4)
#
# line, = ax.plot(x, norm(x))
# ax.margins(0)
# ax.set_ylim(0,1)
#
# im = ax2.imshow(np.atleast_2d(x).T, cmap="Spectral_r", norm=norm, aspect="auto")
# cbar = fig.colorbar(im ,cax=ax3,ax=ax2, orientation="horizontal")


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, min_val=None, max_val=None, name='shiftedcmap'):
    """
    Function to offset the "center" of a colormap. Useful for data with a \
    negative min and positive max and you want the middle of the colormap's dynamic \
    range to be at zero. Adapted from https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    IMPUTS
    -----
        cmap  : The matplotlib colormap to be altered.
        start : Offset from lowest point in the colormap's range.
                Defaults to 0.0 (no lower ofset). Should be between
                0.0 and `midpoint`.
        midpoint : The new center of the colormap. Defaults to
                   0.5 (no shift). Should be between 0.0 and 1.0. In
                   general, this should be  1 - vmax/(vmax + abs(vmin))
                   For example if your data range from -15.0 to +5.0 and
                   you want the center of the colormap at 0.0, `midpoint`
                   should be set to  1 - 5/(5 + 15)) or 0.75
        stop : Offset from highets point in the colormap's range.
               Defaults to 1.0 (no upper ofset). Should be between
               `midpoint` and 1.0.
        min_val : mimimun value of the dataset,
                  only use when 0.0 is pretend to be the midpoint of the colormap
        max_val : maximun value of the dataset,
                  only use when 0.0 is pretend to be the midpoint of the colormap
        name    : Name of the output cmap

    """
    epsilon = 0.001
    # start, stop = 0.0, 1.0
    if min_val is not None and max_val is not None:
        min_val, max_val = min(0.0, min_val), max(0.0, max_val)
        midpoint = 1.0 - max_val/(max_val + abs(min_val))

    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
    # shifted index to match the data
    shift_index = np.hstack([np.linspace(0.0, midpoint, 128, endpoint=False), \
                            np.linspace(midpoint, 1.0, 129, endpoint=True)])
    for ri, si in zip(reg_index, shift_index):
        if abs(si - midpoint) < epsilon:
            r, g, b, a = cmap(0.5) # 0.5 = original midpoint.
        else:
            r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap
