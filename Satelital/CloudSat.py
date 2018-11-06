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


locale.setlocale(locale.LC_ALL, ('en_us','utf-8'))

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

Path_fuentes = '/home/cmcuervol/Fuentes/' # Atlas
# Path_fuentes = '/home/ccuervo/Fuentes/' # SAL


# Colors for graphics SIATA style
gris70 = (112/255., 111/255., 111/255.)

ColorInfo1 = (82 /255., 183/255.,196/255.)
ColorInfo2 = (55 /255., 123/255.,148/255.)
ColorInfo3 = (43 /255.,  72/255.,105/255.)
ColorInfo4 = (32 /255.,  34/255., 72/255.)
ColorInfo5 = (34 /255.,  71/255., 94/255.)
ColorInfo6 = (31 /255., 115/255.,116/255.)
ColorInfo7 = (39 /255., 165/255.,132/255.)
ColorInfo8 = (139/255., 187/255.,116/255.)
ColorInfo9 = (200/255., 209/255., 93/255.)
ColorInfo10 =(249/255., 230/255., 57/255.)


AzulChimba =( 55/255., 132/255., 251/255.)
AzulChimbita =( 16/255., 108/255., 214/255.)
VerdeChimba =(  9/255., 210/255.,  97/255.)
Azul =( 96/255., 200/255., 247/255.)
Naranja =( 240/255., 108/255.,  34/255.)
RojoChimba =( 240/255., 84/255.,  107/255.)
Verdecillo =( 40/255., 225/255.,  200/255.)
Azulillo =( 55/255., 150/255.,  220/255.)
# flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]


# Types of fonts Avenir
#
AvenirHeavy = fm.FontProperties(fname=Path_fuentes+'AvenirLTStd-Heavy.otf')
AvenirBook  = fm.FontProperties(fname=Path_fuentes+'AvenirLTStd-Book.otf')
AvenirBlack = fm.FontProperties(fname=Path_fuentes+'AvenirLTStd-Black.otf')
AvenirRoman = fm.FontProperties(fname=Path_fuentes+'AvenirLTStd-Roman.ttf')

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def describevg(refnum):
    # Describe the vgroup with the given refnum.
    # Open vgroup in read mode.
    vg = v.attach(refnum)
    print "----------------"
    print "name:", vg._name, "class:",vg._class, "tag,ref:",
    print vg._tag, vg._refnum

    # Show the number of members of each main object type.
    print "members: ", vg._nmembers,
    print "datasets:", vg.nrefs(HC.DFTAG_NDG),
    print "vdatas:  ", vg.nrefs(HC.DFTAG_VH),
    print "vgroups: ", vg.nrefs(HC.DFTAG_VG)

    # Read the contents of the vgroup.
    members = vg.tagrefs()

    # Display info about each member.
    index = -1
    for tag, ref in members:
        index += 1
        print "member index", index
        # Vdata tag
        if tag == HC.DFTAG_VH:
            vd = vs.attach(ref)
            nrecs, intmode, fields, size, name = vd.inquire()
            print "  vdata:",name, "tag,ref:",tag, ref
            print "    fields:",fields
            print "    nrecs:",nrecs
            vd.detach()

        # SDS tag
        elif tag == HC.DFTAG_NDG:
            sds = sd.select(sd.reftoindex(ref))
            name, rank, dims, type, nattrs = sds.info()
            print "  dataset:",name, "tag,ref:", tag, ref
            print "    dims:",dims
            print "    type:",type
            sds.endaccess()

        # VS tag
        elif tag == HC.DFTAG_VG:
            vg0 = v.attach(ref)
            print "  vgroup:", vg0._name, "tag,ref:", tag, ref
            vg0.detach()

        # Unhandled tag
        else:
            print "unhandled tag,ref",tag,ref

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
    Extrae los nombres de las variables para un archivo (File)
    """
    # hdfFile = SD.SD(File, mode=1)
    hdfFile = SD(File, mode=1)
    dsets = hdfFile.datasets()
    k = []
    for key in dsets.keys():
        # print key
        k.append(key)
    k.sort()
    hdfFile.end()
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


def newjet():
    """función para crear un nuevo color bar con cero en la mitad, modeificando el color bar jet"""
    jetcmap = plt.cm.get_cmap("jet", 11) #generate a jet map with 11 values
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

coqueto=newjet()
# Coquetisimo = shiftedColorMap(Coqueto,)
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

# def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
#     """
#     Function to offset the "center" of a colormap. Useful for
#     data with a negative min and positive max and you want the
#     middle of the colormap's dynamic range to be at zero.
#
#     Input
#     -----
#       cmap : The matplotlib colormap to be altered
#       start : Offset from lowest point in the colormap's range.
#           Defaults to 0.0 (no lower offset). Should be between
#           0.0 and `midpoint`.
#       midpoint : The new center of the colormap. Defaults to
#           0.5 (no shift). Should be between 0.0 and 1.0. In
#           general, this should be  1 - vmax / (vmax + abs(vmin))
#           For example if your data range from -15.0 to +5.0 and
#           you want the center of the colormap at 0.0, `midpoint`
#           should be set to  1 - 5/(5 + 15)) or 0.75
#       stop : Offset from highest point in the colormap's range.
#           Defaults to 1.0 (no upper offset). Should be between
#           `midpoint` and 1.0.
#     """
#     cdict = {
#         'red': [],
#         'green': [],
#         'blue': [],
#         'alpha': []
#     }
#
#     # regular index to compute the colors
#     reg_index = np.linspace(start, stop, 257)
#
#     # shifted index to match the data
#     shift_index = np.hstack([
#         np.linspace(0.0, midpoint, 128, endpoint=False),
#         np.linspace(midpoint, 1.0, 129, endpoint=True)
#     ])
#
#     for ri, si in zip(reg_index, shift_index):
#         r, g, b, a = cmap(ri)
#
#         cdict['red'].append((si, r, r))
#         cdict['green'].append((si, g, g))
#         cdict['blue'].append((si, b, b))
#         cdict['alpha'].append((si, a, a))
#
#     newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
#     plt.register_cmap(cmap=newcmap)
#
#     return newcmap
#


# def GranulePlot(Data, lats, lons, cmap, levels, labelData, fecha, name, centredcmap=False,DEM=None):
#     """
#     Plot a CloudSat data "Granule" is defined as one orbit. A granule starts
#     at the first profile that falls on or past the equator on the descending node.
#     (note: Approximately 20 seconds, or 125 profiles, will be appended to the
#     beginning and end of each granule).
#
#     IMPUTS:
#     Data: array 2D of data to plot
#     lats: array of latitudes
#     lons: array of longitudes
#     cmap: color map
#     levels : type of scale, LogNorm = logarithmic, SymLogNorm= symmetric
#              logarithmic, Norm = lineal, Squeezed = compress in some values
#              or array of levels define by user. If is None dont use norm
#     lableData: Name to show in the data plot
#     fecha : datetime of the data
#     name: name of the outputfile
#     centredcmap : calls the function shiftedColorMap to centred the cmap given
#     DEM : Digital elevation model,
#
#     OUTPUTS
#     file save in the path_fig folder
#     """
#
#     def OrganizaLon(longitude):
#         """
#         reacomodate the longitudes values of a satellite track for graph it easily
#         """
#         n   = np.where(longitude>0)[0] # Indixes with longitudes positives
#         a   = Salto(n,1)         # When occurs the change to negative
#         idx = n[a+1]             # Index when is positive again
#         copi = longitude.copy()
#         copi[idx:] = copi[idx:]-360 # Change to negative te last positive portion of data
#
#         return copi
#
#     #limpiar guevonadas
#     plt.cla()
#     plt.clf()
#     plt.close('all')
#
#     fig = plt.figure(figsize=(10,10))
#     ax1 = fig.add_subplot(2,1,1)
#
#     lons = OrganizaLon(lons)+360
#     m = Basemap(ax=ax1,llcrnrlat=-90, llcrnrlon=lons.min(),
#                 urcrnrlat=90, urcrnrlon=lons.max())
#     # m.etopo()
#     # m.drawcoastlines()
#     m.shadedrelief()
#     # m.bluemarble()
#     # m.drawcoastlines()
#     # m.drawcountries()
#     # m.drawmapboundary(fill_color='cyan')
#     # m.fillcontinents(color='coral')
#     # draw parallels.
#     parallels = np.arange(-90.,91.,30.)
#     m.drawparallels(parallels,labels=[1,0,0,0],fontsize=7)
#     # draw meridians
#     meridians = np.arange(-360.,721.,60.)
#     m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=7)
#     ax1.set_title('CloudSat track',fontproperties=AvenirRoman, color=gris70, fontsize=15)
#     # ny = data.shape[0]; nx = data.shape[1]
#     # lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
#     # x, y = m(lons, lats)
#     x, y = m(lons, lats)
#     # m.scatter(x, y,color=RojoChimba)
#     m.plot(x, y,'--', color=RojoChimba,lw=2)
#
#     ax2 = fig.add_subplot(2,1,2)
#     t,h = np.meshgrid(range(Data.shape[0]), range(Data.shape[1]))
#     h = h*0.24  # There are 125 vertical bins, each one approximately 240m thick
#     t = t*1.7   # The CloudSat data footprint is approximately 1.7km along-track by 1.3km across-track.
#
#     # cmap = plt.cm.jet
#
#
#
#     ax1.set_title(fecha.strftime('%Y-%m-%d %H:%M:%S'), loc='right', fontsize=10, color=gris70)
#     minimo = np.nanmin(Data)
#     maximo = np.nanmax(Data)
#     levels = np.linspace(minimo, maximo, 100)
#     if   levels == 'LogNorm':
#         norml = colors.LogNorm(minimo, maximo)
#     elif levels == 'SymLogNorm':
#         norml = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=minimo, vmax=maximo)
#     elif levels == 'Norm':
#         norml = colors.Normalize(minimo, maximo)
#     elif levels is not None: # Exist the option None check if works bad
#         norml = colors.BoundaryNorm(boundaries=levels, ncolors=256)
#
#     if centredcmap == True:
#         cmap = shiftedColorMap(cmap, min_val=minimo, max_val=maximo)
#     #contorno de colores
#     if levels == None:
#         cs = ax2.contourf(t, h, Data.T[::-1,:], cmap=cmap)
#     else:
#         cs = ax2.contourf(t, h, Data.T[::-1,:], norm=norml, cmap=cmap)
#
#     cbar_ax = fig.add_axes([0.85, 0.1, 0.015, 0.4])
#     cbar= fig.colorbar(cs, cax=cbar_ax, orientation='vertical')
#     cbar.ax.set_ylabel(labelData,fontproperties=AvenirRoman, color=gris70,fontsize=15)
#
#     ax2.set_ylabel('Height [km]',fontproperties=AvenirRoman, color=gris70, fontsize=15)
#     ax2.set_xlabel('Distance along-track [km]',fontproperties=AvenirRoman, color=gris70, fontsize=15)
#     if DEM is not None:
#         ax2.plot(DEM*1E-3,lw=2, color= gris70)
#
#     plt.subplots_adjust(left=0.125, bottom=0.1, right=0.8, top=0.95, wspace=0.2, hspace=0.1)
#     plt.savefig(path_fig + name, transparent=True)
#

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def cmapeador(colrs=None, levels=None, name='coqueto'):
    """
    Make a new color map with the colors and levels given, adapted from Julián Sepulveda & Carlos Hoyos

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



def GranulePloter(Data, lats, lons, cmap, norm, cmaplev='default', extend='neither', \
                  scale='linear', labelData='Variable', fecha='15-01-1992', name='Prueba' ,DEM=None):
    """
    Plot a CloudSat data "Granule" is defined as one orbit. A granule starts
    at the first profile that falls on or past the equator on the descending node.
    (note: Approximately 20 seconds, or 125 profiles, will be appended to the
    beginning and end of each granule).

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
    DEM   : Digital elevation model, correspond for the "Granule"

    OUTPUTS
    file save in the path_fig folder
    """

    def OrganizaLon(longitude):
        """
        reacomodate the longitudes values of a satellite track for graph it easily
        """
        n   = np.where(longitude>0)[0] # Indixes with longitudes positives
        a   = Salto(n,1)         # When occurs the change to negative
        idx = n[a+1]             # Index when is positive again
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

    ax2 = fig.add_subplot(2,1,2)
    t,h = np.meshgrid(range(Data.shape[0]), range(Data.shape[1]))
    h = h*0.24  # There are 125 vertical bins, each one approximately 240m thick
    t = t*1.7   # The CloudSat data footprint is approximately 1.7km along-track by 1.3km across-track.
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
    if DEM is not None:
        ax2.plot(DEM*1E-3,lw=2, color= gris70)

    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.8, top=0.95, wspace=0.2, hspace=0.1)
    plt.savefig(path_fig + name, transparent=True)

# GranulePlot(REP, Lat, Lon, 'Received echo powers', 'Preuba.png')

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

path_fig  = '/home/cmcuervol/A-Train/CloudSat/Figures/'
# path = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-PRECIP-COLUMN.P1_R05/2006/185/'
path_Data = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1B-CPR.P_R05/2006/153/'
path_Geop = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2B-GEOPROF.P1_R05/2006/153/'
path_Rain = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-RAIN-PROFILE.P1_R05/2006/153/'
# path_PptC = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-PRECIP-COLUMN.P1_R05/2006/153/'
path_FRLK = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1b-cpr-fl.p_r04/2012/037/'

path_Ccls = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS.P_R04/'
path_LIDR = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS-LIDAR.P_R04/'
path_TRMM = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2D-CLOUDSAT-TRMM.P_R04/'
path_CWCr = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RO.P_R04/'
path_CWCo = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RVOD.P_R04/'

# files_Data = Listador(path_Data, final='hdf' )
# files_Geop = Listador(path_Geop, final='hdf' )
# files_Rain = Listador(path_Rain, final='hdf' )
# files_PptC = Listador(path_PptC, final='hdf' )
# files_FRLK = Listador(path_FRLK, final='zip' )

# files_Ccls_zip = Listador(path_Ccls, final='.zip')
# files_LIDR_zip = Listador(path_LIDR, final='.zip')
# files_TRMM_zip = Listador(path_TRMM, final='.zip')
# files_CWCr_zip = Listador(path_CWCr, final='.zip')
# files_CWCo_zip = Listador(path_CWCo, final='.zip')
#

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


paths_zip = [path_Ccls, path_LIDR, path_TRMM, path_CWCr, path_CWCo]

# # Esto solo se hace una vez, esperar a que estén completos todos los archivos
# for path_product in paths_zip:
#     years = Listador(path_product, inicio='2')
#     for year in years:
#         days = Listador(path_product+year)
#         for day in days:
#             # print path_product+year+'/'+day
#             # Unzipeador(path_product+year+'/'+day)
#


# files_Ccls = Listador(path_Ccls, final='.hdf')
# files_LIDR = Listador(path_LIDR, final='.hdf')
# files_TRMM = Listador(path_TRMM, final='.hdf')
# files_CWCr = Listador(path_CWCr, final='.hdf')
# files_CWCo = Listador(path_CWCo, final='.hdf')
#


# NameHDF_Data = files_Data[0]
# NameHDF_Geop = files_Geop[0]
# NameHDF_Rain = files_Rain[0]
# # NameHDF_PptC = files_PptC[0]
# # NameHDF_FRLK = files_FRLK[-1]
# # NameHDF_FRLK = '2012037001654_30725_CS_1B-CPR-FL_GRANULE_P_R04_E05.hdf'


# Variables_CPR = HDFvars(path_Data+NameHDF_Data)
# Variables_Geo = HDFvars(path_Geop+NameHDF_Geop)
# Variables_Ran = HDFvars(path_Rain+NameHDF_Rain)
# # Variables_PPT = HDFvars(path_PptC+NameHDF_PptC)
# # Variables_Frl = HDFvars(path_FRLK+NameHDF_FRLK)

#
# variables_TRMM = HDFvars(path_TRMM+files_TRMM[0])
#
# Lat = DesHDF(path_TRMM+files_TRMM[0], 'CPR Height', '')
# print Variables_CPR
# print Variables_Geo
# print Variables_Ran
# # print Variables_PPT
# # print Variables_Frl



# Lat = HDFread(path_Data+NameHDF_Data, 'Latitude')
# Lon = HDFread(path_Data+NameHDF_Data, 'Longitude')
# tim = HDFread(path_Data+NameHDF_Data, 'Profile_time')
# UTC = HDFread(path_Data+NameHDF_Data, 'UTC_start')
# DEM = HDFread(path_Data+NameHDF_Data, 'DEM_elevation')
#
# # DEM = DEM.ravel().astype(float)
# DEM = DEM.astype(float)
# DEM[DEM==-9999]=np.nan
#
# REP = DesHDF(path_Data+NameHDF_Data, 'ReceivedEchoPowers' , forma='')
# REP[REP==-9999] = np.nan
#
#




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

path_Data = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/1B-CPR.P_R05/'
path_Geop = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2B-GEOPROF.P1_R05/'
path_Rain = '/home/cmcuervol/A-Train/CloudSat/ftp1.cloudsat.cira.colostate.edu/2C-RAIN-PROFILE.P1_R05/'

path_Ccls = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS.P_R04/'
path_LIDR = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CLDCLASS-LIDAR.P_R04/'
path_TRMM = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2D-CLOUDSAT-TRMM.P_R04/'
path_CWCr = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RO.P_R04/'
path_CWCo = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-CWC-RVOD.P_R04/'
path_FLUX = '/mnt/Sanduchera/ftp1.cloudsat.cira.colostate.edu/2B-FLXHR.P2_R04/'


stdcolors = [(255, 255, 255),(0, 255, 255), (0, 0, 255),(70, 220, 45),\
             (44, 141, 29),(255,255,75),(255,142,0),(255,0,0),(128,0,128),\
             (102,0,102),(255, 153, 255)]

levels_LiqWat = np.array([0, 25,50,75,100, 150, 200, 250, 500, 750, 1000])*1E-1
#
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#
# year = '2006'
# days = Listador(path_Data+year)
# days = [days[1]]
#
# for day in days:
#     try:
#         files_Data = Listador(path_Data+year+'/'+day, final='hdf' )
#         Data = True
#     except:
#         Data = False
#     try:
#         files_Geop = Listador(path_Geop+year+'/'+day, final='hdf' )
#         Geop = True
#     except:
#         Geop = False
#     try:
#         files_Rain = Listador(path_Rain+year+'/'+day, final='hdf' )
#         Rain = True
#     except:
#         Rain = False
#     try:
#         Unzipeador(path_Ccls+year+'/'+day+'/', erase=False)
#         files_Ccls = Listador(path_Ccls+year+'/'+day, final='hdf' )
#         CloudClass = True
#     except:
#         CloudClass = False
#     try:
#         Unzipeador(path_LIDR+year+'/'+day+'/', erase=False)
#         files_LIDR = Listador(path_LIDR+year+'/'+day, final='hdf' )
#         LIDAR = True
#     except:
#         LIDAR = False
#     try:
#         Unzipeador(path_CWCr+year+'/'+day+'/', erase=False)
#         files_CWCr = Listador(path_CWCr+year+'/'+day, final='hdf' )
#         CWC_RO = True
#     except:
#         CWC_RO = False
#     try:
#         Unzipeador(path_CWCo+year+'/'+day+'/', erase=False)
#         files_CWCo = Listador(path_CWCo+year+'/'+day, final='hdf' )
#         CWC_RVOD = True
#     except:
#         CWC_RVOD =False
#     # try:
#     #     Unzipeador(path_TRMM+year+'/'+day+'/', erase=False)
#     #     files_TRMM = Listador(path_TRMM+year+'/'+day, final='hdf' )
#     #     TRMM = True
#     # except:
#     #     TRMM =False
#     try:
#         Unzipeador(path_FLUX+year+'/'+day+'/', erase=False)
#         files_FLUX = Listador(path_FLUX+year+'/'+day, final='hdf' )
#         FLUX = True
#     except:
#         FLUX =False
#
#     for i  in range(len(files_Data)):
#     # for i  in range(1):
#     # en i =1 para el día 154 del 2006
#         # i = 1
#         fecha = GranuleTime(files_Data[i])
#         print fecha
#         if Data == True:
#             Lat = HDFread(path_Data+year+'/'+day+'/'+files_Data[i], 'Latitude')
#             Lon = HDFread(path_Data+year+'/'+day+'/'+files_Data[i], 'Longitude')
#             print path_Data+year+'/'+day+'/'+files_Data[i], 'ReceivedEchoPowers'
#             REP = DesHDF(path_Data+year+'/'+day+'/'+files_Data[i], 'ReceivedEchoPowers' , forma='')
#             REP = np.ma.masked_where(REP==-9999,REP)
#
#             # levels_echo = np.array([1E-16, 1E-15,1E-14,1E-13,1E-12,1E-11,1E-10, 1E-9, 1E-8, 1E-7, 1E-6])
#             levels_echo = np.logspace(np.log10(1E-16), np.log10(1E-6),11)
#             # cmap_echo, norm_echo = cmapeador(stdcolors,levels_echo, 'echo')
#             norm_echo = colors.LogNorm(1E-16, 1E-6)
#             GranulePloter(REP, Lat, Lon, plt.cm.jet, norm_echo, levels_echo, 'max','logarithmic',\
#                           'Received echo powers',fecha, 'ReceivedEchoPowers'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png' )
#
#         if Geop == True:
#             REF = DesHDF(path_Geop+year+'/'+day+'/'+files_Geop[i], 'Radar_Reflectivity' , forma='')
#             REF = np.ma.masked_where(REF==-8888,REF)
#
#             levels_REF = np.arange(-60,61,10)
#             # norm_REF   = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-40, vmax=40)
#             norm_REF   = colors.Normalize(vmin=-60, vmax=60)
#             GranulePloter(REF*1E-2, Lat, Lon, coqueto, norm_REF,levels_REF, 'neither', 'linear',\
#                          'Radar reflectivity [dBZ]', fecha,'Radar_Reflectivity'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#         if Rain == True:
#             CLW = DesHDF(path_Rain+year+'/'+day+'/'+files_Rain[i], 'cloud_liquid_water' ,   forma='').astype(float)
#             PIC = DesHDF(path_Rain+year+'/'+day+'/'+files_Rain[i], 'precip_ice_water' ,     forma='').astype(float)
#             PLC = DesHDF(path_Rain+year+'/'+day+'/'+files_Rain[i], 'precip_liquid_water' ,  forma='').astype(float)
#             MRF = DesHDF(path_Rain+year+'/'+day+'/'+files_Rain[i], 'modeled_reflectivity' , forma='').astype(float)
#             # PPT = DesHDF(path_PptC+files_PptC[i], 'unused' , forma='').astype(float)
#             CLW = np.ma.masked_where(CLW==-9999,CLW)
#             PLC = np.ma.masked_where(PLC==-9999,PLC)
#             PIC = np.ma.masked_where(PIC==-9999,PIC)
#             MRF = np.ma.masked_where(MRF==-9999,MRF)
#
#
#             cmap_LiqWat, norm_LiqWat = cmapeador(stdcolors,levels_LiqWat, 'LiqWat')
#             GranulePloter(CLW*1E-1, Lat, Lon, cmap_LiqWat, norm_LiqWat, np.delete(levels_LiqWat, [1,3]), 'max', 'linear',\
#                           'Cloud liquid water [mm]',fecha, 'cloud_liquid_water'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png' )
#
#             levels_IceWat = levels_LiqWat*1
#             cmap_IceWat, norm_IceWat = cmapeador(stdcolors,levels_IceWat, 'IceWat')
#             GranulePloter(PIC*1E-1, Lat, Lon, cmap_IceWat, norm_IceWat, np.delete(levels_IceWat, [1,3]), 'max','linear',\
#                           'Preciptable ice water [mm]',fecha, 'precip_ice_water'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png' )
#
#             levels_PptWat = levels_LiqWat*1
#             cmap_PptWat, norm_PptWat = cmapeador(stdcolors,levels_PptWat, 'PptWat')
#             GranulePloter(PLC*1E-1, Lat, Lon, cmap_PptWat, norm_PptWat, np.delete(levels_PptWat, [1,3]), 'max','linear',\
#                           'Preciptable liquid water [mm]',fecha, 'precip_liquid_water'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png' )
#
#             levels_MRF = np.arange(-60,61,10)
#             # norm_MRF   = colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-40, vmax=40)
#             norm_MRF   = colors.Normalize(vmin=-60, vmax=60)
#             GranulePloter(MRF*1E-2, Lat, Lon, coqueto, norm_MRF,levels_MRF, 'neither', 'linear',\
#                          'Modeled reflectivity [dBZ]', fecha,'Modeled_Reflectivity'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#         if CloudClass == True:
#             # ==============================================================================
#             #                                 Cloud Calss
#             # ==============================================================================
#
#             # Height = DesHDF(path_Ccls+year+'/'+day+'/'+files_Ccls[0], 'Height')
#             # absoluto = max(abs(Height.min()),abs(Height.max()))
#             # norm_H   = colors.Normalize(vmin=-absoluto, vmax=absoluto)
#             # GranulePloter(Height, Lat, Lon, coqueto, norm_H,'default', 'neither', 'linear',\
#             #              'Height', fecha,'Height'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             # c = Cloud.ravel()
#             # a = list(set(c))
#             Cloud =  DesHDF(path_Ccls+year+'/'+day+'/'+files_Ccls[i], 'cloud_scenario')
#             # Hallar el No cloud al parecer es el menor a 2500
#             norm_Cloud   = colors.Normalize(vmin=Cloud.min(), vmax=Cloud.max())
#             GranulePloter(Cloud, Lat, Lon, coqueto, norm_Cloud,'default', 'neither', 'linear',\
#                          'Cloud', fecha,'Cloud'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#         if LIDAR == True:
#             # ==============================================================================
#             #                              Cloud Calss LIDAR
#             # ==============================================================================
#
#             CloudFrq = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'CloudFraction')
#             CloudFrq = np.ma.masked_where(CloudFrq==-99,CloudFrq)
#             norm_CloudFrq   = colors.Normalize(vmin=np.nanmin(CloudFrq), vmax=np.nanmax(CloudFrq))
#             GranulePloter(CloudFrq, Lat, Lon, plt.cm.jet, norm_CloudFrq,'default', 'neither', 'linear',\
#                          'Cloud fraction', fecha,'CloudFraction'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
#             CloudLayB = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'CloudLayerBase')
#             CloudLayB = np.ma.masked_where(CloudLayB==-99,CloudLayB)
#             norm_CloudLayB   = colors.Normalize(vmin=np.nanmin(CloudLayB), vmax=np.nanmax(CloudLayB))
#             GranulePloter(CloudLayB, Lat, Lon, plt.cm.jet, norm_CloudLayB,'default', 'neither', 'linear',\
#                          'Cloud layer base', fecha,'CloudLayerBase'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             CloudLayT = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'CloudLayerTop')
#             CloudLayT = np.ma.masked_where(CloudLayT==-99,CloudLayT)
#             norm_CloudLayT   = colors.Normalize(vmin=np.nanmin(CloudLayT), vmax=np.nanmax(CloudLayT))
#             GranulePloter(CloudLayT, Lat, Lon, plt.cm.jet, norm_CloudLayT,'default', 'neither', 'linear',\
#                          'Cloud layer top', fecha,'CloudLayerTop'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
#             CloudLayty = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'CloudLayerType')
#             CloudLayty = np.ma.masked_where(CloudLayty==0,CloudLayty)
#             norm_CloudLayty   = colors.Normalize(vmin=np.nanmin(CloudLayty), vmax=np.nanmax(CloudLayty))
#             GranulePloter(CloudLayty, Lat, Lon, plt.cm.jet, norm_CloudLayty,'default', 'neither', 'default',\
#                          'Cloud layer type', fecha,'CloudLayerType'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             CloudPhase = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'CloudPhase')
#             norm_CloudPhase   = colors.Normalize(vmin=CloudPhase.min(), vmax=CloudPhase.max())
#             GranulePloter(CloudPhase, Lat, Lon, plt.cm.jet, norm_CloudPhase,'default', 'neither', 'default',\
#                          'Cloud phase', fecha,'CloudPhase'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             # CloudLayBf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'LayerBaseFlag')
#             # norm_CloudLayBf   = colors.Normalize(vmin=CloudLayBf.min(), vmax=CloudLayBf.max())
#             # GranulePloter(CloudLayBf, Lat, Lon, plt.cm.gray, norm_CloudLayBf,'default', 'neither', 'linear',\
#             #              'Layer base flag', fecha,'LayerBaseFlag'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             # CloudLayTf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'LayerTopFlag')
#             # norm_CloudLayTf   = colors.Normalize(vmin=CloudLayTf.min(), vmax=CloudLayTf.max())
#             # GranulePloter(CloudLayTf, Lat, Lon, plt.cm.gray, norm_CloudLayTf,'default', 'neither', 'linear',\
#             #              'Layer top flag', fecha,'LayerTopFlag'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             CloudPPTf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'PrecipitationFlag')
#             norm_CloudPPTf   = colors.Normalize(vmin=CloudPPTf.min(), vmax=CloudPPTf.max())
#             GranulePloter(CloudPPTf, Lat, Lon, plt.cm.gray, norm_CloudPPTf,'default', 'neither', 'default',\
#                          'Precipitation flag', fecha,'PrecipitationFlag'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             CloudWlayT = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[i], 'Water_layer_top')
#             CloudWlayT = np.ma.masked_where(CloudWlayT==-9,CloudWlayT)
#             norm_CloudWlayT   = colors.Normalize(vmin=np.nanmin(CloudWlayT), vmax=np.nanmax(CloudWlayT))
#             GranulePloter(CloudWlayT, Lat, Lon, plt.cm.gray, norm_CloudWlayT,'default', 'neither', 'default',\
#                          'Water layer top', fecha,'WaterLayerTop'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#         if CWC_RO == True:
#             # ==============================================================================
#             #                          Cloud water content Radar Only
#             # ==============================================================================
#
#             ICEdisWidthpar_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_ice_distrib_width_param')
#             ICEdisWidthpar_RO = np.ma.masked_where(ICEdisWidthpar_RO<=0,ICEdisWidthpar_RO)
#             norm_Ice_withdis   = colors.Normalize(vmin=np.nanmin(ICEdisWidthpar_RO), vmax=np.nanmax(ICEdisWidthpar_RO))
#             GranulePloter(ICEdisWidthpar_RO, Lat, Lon, plt.cm.jet, norm_Ice_withdis,'default', 'neither', 'linear',\
#                          'Ice width distribution', fecha,'IceWidth_dist_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqdisWidthpar_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_liq_distrib_width_param')
#             LiqdisWidthpar_RO = np.ma.masked_where(LiqdisWidthpar_RO<=0,LiqdisWidthpar_RO)
#             norm_Liq_withdis   = colors.Normalize(vmin=np.nanmin(LiqdisWidthpar_RO), vmax=np.nanmax(LiqdisWidthpar_RO))
#             GranulePloter(LiqdisWidthpar_RO, Lat, Lon, plt.cm.jet, norm_Liq_withdis,'default', 'neither', 'linear',\
#                          'Liquid width distribution', fecha,'LiquidWidth_dist_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             IceWaterCont_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_ice_water_content')
#             IceWaterCont_RO = np.ma.masked_where(IceWaterCont_RO<=0,IceWaterCont_RO)
#             norm_IceWat   = colors.Normalize(vmin=np.nanmin(IceWaterCont_RO), vmax=np.nanmax(IceWaterCont_RO))
#             GranulePloter(IceWaterCont_RO, Lat, Lon, plt.cm.jet, norm_IceWat,'default', 'neither', 'linear',\
#                          'Ice water content', fecha,'IceWaterCont_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqWaterCont_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_liq_water_content')
#             LiqWaterCont_RO = np.ma.masked_where(LiqWaterCont_RO<=0,LiqWaterCont_RO)
#             norm_LiqWat   = colors.Normalize(vmin=np.nanmin(LiqWaterCont_RO), vmax=np.nanmax(LiqWaterCont_RO))
#             GranulePloter(LiqWaterCont_RO, Lat, Lon, plt.cm.jet, norm_LiqWat,'default', 'neither', 'linear',\
#                          'Liquid water content', fecha,'LiqWaterCont_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             IceRadius_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_ice_effective_radius')
#             IceRadius_RO = np.ma.masked_where(IceRadius_RO<=0,IceRadius_RO)
#             norm_IceRadius   = colors.Normalize(vmin=np.nanmin(IceRadius_RO), vmax=np.nanmax(IceRadius_RO))
#             GranulePloter(IceRadius_RO, Lat, Lon, plt.cm.jet, norm_IceRadius,'default', 'neither', 'linear',\
#                          'Ice effective radius', fecha,'IceRadius_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqRadius_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_liq_effective_radius')
#             LiqRadius_RO = np.ma.masked_where(LiqRadius_RO<=0,LiqRadius_RO)
#             norm_LiqRadius   = colors.Normalize(vmin=np.nanmin(LiqRadius_RO), vmax=np.nanmax(LiqRadius_RO))
#             GranulePloter(LiqRadius_RO, Lat, Lon, plt.cm.jet, norm_LiqRadius,'default', 'neither', 'linear',\
#                          'Liquid effective radius', fecha,'LiqRadius_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             IceNumber_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_ice_number_conc')
#             IceNumber_RO = np.ma.masked_where(IceNumber_RO<=0,IceNumber_RO)
#             norm_IceNumber   = colors.Normalize(vmin=np.nanmin(IceNumber_RO), vmax=np.nanmax(IceNumber_RO))
#             GranulePloter(IceNumber_RO, Lat, Lon, plt.cm.jet, norm_IceNumber,'default', 'neither', 'linear',\
#                          'Ice number concentration', fecha,'IceNumber_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqNumber_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[i], 'RO_liq_number_conc')
#             LiqNumber_RO = np.ma.masked_where(LiqNumber_RO<=0,LiqNumber_RO)
#             norm_LiqNumber   = colors.Normalize(vmin=np.nanmin(LiqNumber_RO), vmax=np.nanmax(LiqNumber_RO))
#             GranulePloter(LiqNumber_RO, Lat, Lon, plt.cm.jet, norm_LiqNumber,'default', 'neither', 'linear',\
#                          'Liquid number concentration', fecha,'LiqNumber_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
#         if CWC_RVOD == True:
#             # ==============================================================================
#             #                     Cloud water content Radar vs Optical Depth
#             # ==============================================================================
#
#
#             ICEdisWidthpar_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_ice_distrib_width_param')
#             ICEdisWidthpar_RVOD = np.ma.masked_where(ICEdisWidthpar_RVOD<=0,ICEdisWidthpar_RVOD)
#             norm_Ice_withdis   = colors.Normalize(vmin=np.nanmin(ICEdisWidthpar_RVOD), vmax=np.nanmax(ICEdisWidthpar_RVOD))
#             GranulePloter(ICEdisWidthpar_RVOD, Lat, Lon, plt.cm.jet, norm_Ice_withdis,'default', 'neither', 'linear',\
#                          'Ice width distribution', fecha,'IceWidth_dist_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqdisWidthpar_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_liq_distrib_width_param')
#             LiqdisWidthpar_RVOD = np.ma.masked_where(LiqdisWidthpar_RVOD<=0,LiqdisWidthpar_RVOD)
#             norm_Ice_withdis   = colors.Normalize(vmin=np.nanmin(LiqdisWidthpar_RVOD), vmax=np.nanmax(LiqdisWidthpar_RVOD))
#             GranulePloter(LiqdisWidthpar_RVOD, Lat, Lon, plt.cm.jet, norm_Ice_withdis,'default', 'neither', 'linear',\
#                          'Liquid width distribution', fecha,'LiquidWidth_dist_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             IceWaterCont_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_ice_water_content')
#             IceWaterCont_RVOD = np.ma.masked_where(IceWaterCont_RVOD<=0,IceWaterCont_RVOD)
#             norm_IceWat   = colors.Normalize(vmin=np.nanmin(IceWaterCont_RVOD), vmax=np.nanmax(IceWaterCont_RVOD))
#             GranulePloter(IceWaterCont_RVOD, Lat, Lon, plt.cm.jet, norm_IceWat,'default', 'neither', 'linear',\
#                          'Ice water content', fecha,'IceWaterCont_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqWaterCont_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_liq_water_content')
#             LiqWaterCont_RVOD = np.ma.masked_where(LiqWaterCont_RVOD<=0,LiqWaterCont_RVOD)
#             norm_LiqWat   = colors.Normalize(vmin=np.nanmin(LiqWaterCont_RVOD), vmax=np.nanmax(LiqWaterCont_RVOD))
#             GranulePloter(LiqWaterCont_RVOD, Lat, Lon, plt.cm.jet, norm_LiqWat,'default', 'neither', 'linear',\
#                          'Liquid water content', fecha,'LiqWaterCont_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             IceRadius_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_ice_effective_radius')
#             IceRadius_RVOD = np.ma.masked_where(IceRadius_RVOD<=0,IceRadius_RVOD)
#             norm_IceRadius   = colors.Normalize(vmin=np.nanmin(IceRadius_RVOD), vmax=np.nanmax(IceRadius_RVOD))
#             GranulePloter(IceRadius_RVOD, Lat, Lon, plt.cm.jet, norm_IceRadius,'default', 'neither', 'linear',\
#                          'Ice effective radius', fecha,'IceRadius_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqRadius_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_liq_effective_radius')
#             LiqRadius_RVOD = np.ma.masked_where(LiqRadius_RVOD<=0,LiqRadius_RVOD)
#             norm_LiqRadius   = colors.Normalize(vmin=np.nanmin(LiqRadius_RVOD), vmax=np.nanmax(LiqRadius_RVOD))
#             GranulePloter(LiqRadius_RVOD, Lat, Lon, plt.cm.jet, norm_LiqRadius,'default', 'neither', 'linear',\
#                          'Liquid effective radius', fecha,'LiqRadius_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             IceNumber_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_ice_number_conc')
#             IceNumber_RVOD = np.ma.masked_where(IceNumber_RVOD<=0,IceNumber_RVOD)
#             norm_IceNumber   = colors.Normalize(vmin=np.nanmin(IceNumber_RVOD), vmax=np.nanmax(IceNumber_RVOD))
#             GranulePloter(IceNumber_RVOD, Lat, Lon, plt.cm.jet, norm_IceNumber,'default', 'neither', 'linear',\
#                          'Ice number concentration', fecha,'IceNumber_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#             LiqNumber_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[i], 'RVOD_liq_number_conc')
#             LiqNumber_RVOD = np.ma.masked_where(LiqNumber_RVOD<=0,LiqNumber_RVOD)
#             norm_LiqNumber   = colors.Normalize(vmin=np.nanmin(LiqNumber_RVOD), vmax=np.nanmax(LiqNumber_RVOD))
#             GranulePloter(LiqNumber_RVOD, Lat, Lon, plt.cm.jet, norm_LiqNumber,'default', 'neither', 'linear',\
#                          'Liquid number concentration', fecha,'LiqNumber_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

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




#
#
# # # HeightLID = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'Height') # Es igual
# # CloudFrq = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudFraction')
# # CloudLayB = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudLayerBase')
# # CloudLayT = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudLayerTop')
# # CloudLayty = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudLayerType')
# # CloudPhase = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudPhase')
# # CloudLayBf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'LayerBaseFlag')
# # CloudLayTf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'LayerTopFlag')
# # CloudPPTf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'PrecipitationFlag')
# # CloudWlayT = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'Water_layer_top')
# #
#
#
#
# fecha = GranuleTime(files_Data[0])
# Lat = HDFread(path_Data+year+'/'+day+'/'+files_Data[0], 'Latitude')
# Lon = HDFread(path_Data+year+'/'+day+'/'+files_Data[0], 'Longitude')
#
# # ==============================================================================
# #                                 Cloud Calss
# # ==============================================================================
#
# # Height = DesHDF(path_Ccls+year+'/'+day+'/'+files_Ccls[0], 'Height')
# # absoluto = max(abs(Height.min()),abs(Height.max()))
# # norm_H   = colors.Normalize(vmin=-absoluto, vmax=absoluto)
# # GranulePloter(Height, Lat, Lon, coqueto, norm_H,'default', 'neither', 'linear',\
# #              'Height', fecha,'Height'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# # c = Cloud.ravel()
# # a = list(set(c))
# Cloud =  DesHDF(path_Ccls+year+'/'+day+'/'+files_Ccls[0], 'cloud_scenario')
# # Hallar el No cloud al parecer es el menor a 2500
# norm_Cloud   = colors.Normalize(vmin=Cloud.min(), vmax=Cloud.max())
# GranulePloter(Cloud, Lat, Lon, coqueto, norm_Cloud,'default', 'neither', 'linear',\
#              'Cloud', fecha,'Cloud'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# # ==============================================================================
# #                              Cloud Calss LIDAR
# # ==============================================================================
#
# CloudFrq = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudFraction')
# CloudFrq = np.ma.masked_where(CloudFrq==-99,CloudFrq)
# norm_CloudFrq   = colors.Normalize(vmin=np.nanmin(CloudFrq), vmax=np.nanmax(CloudFrq))
# GranulePloter(CloudFrq, Lat, Lon, plt.cm.jet, norm_CloudFrq,'default', 'neither', 'linear',\
#              'Cloud fraction', fecha,'CloudFraction'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
# CloudLayB = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudLayerBase')
# CloudLayB = np.ma.masked_where(CloudLayB==-99,CloudLayB)
# norm_CloudLayB   = colors.Normalize(vmin=np.nanmin(CloudLayB), vmax=np.nanmax(CloudLayB))
# GranulePloter(CloudLayB, Lat, Lon, plt.cm.jet, norm_CloudLayB,'default', 'neither', 'linear',\
#              'Cloud layer base', fecha,'CloudLayerBase'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# CloudLayT = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudLayerTop')
# CloudLayT = np.ma.masked_where(CloudLayT==-99,CloudLayT)
# norm_CloudLayT   = colors.Normalize(vmin=np.nanmin(CloudLayT), vmax=np.nanmax(CloudLayT))
# GranulePloter(CloudLayT, Lat, Lon, plt.cm.jet, norm_CloudLayT,'default', 'neither', 'linear',\
#              'Cloud layer top', fecha,'CloudLayerTop'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
# CloudLayty = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudLayerType')
# CloudLayty = np.ma.masked_where(CloudLayty==0,CloudLayty)
# norm_CloudLayty   = colors.Normalize(vmin=np.nanmin(CloudLayty), vmax=np.nanmax(CloudLayty))
# GranulePloter(CloudLayty, Lat, Lon, plt.cm.jet, norm_CloudLayty,'default', 'neither', 'default',\
#              'Cloud layer type', fecha,'CloudLayerType'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# CloudPhase = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'CloudPhase')
# norm_CloudPhase   = colors.Normalize(vmin=CloudPhase.min(), vmax=CloudPhase.max())
# GranulePloter(CloudPhase, Lat, Lon, plt.cm.gray, norm_CloudPhase,'default', 'neither', 'default',\
#              'Cloud phase', fecha,'CloudPhase'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# # CloudLayBf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'LayerBaseFlag')
# # norm_CloudLayBf   = colors.Normalize(vmin=CloudLayBf.min(), vmax=CloudLayBf.max())
# # GranulePloter(CloudLayBf, Lat, Lon, plt.cm.gray, norm_CloudLayBf,'default', 'neither', 'linear',\
# #              'Layer base flag', fecha,'LayerBaseFlag'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# # CloudLayTf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'LayerTopFlag')
# # norm_CloudLayTf   = colors.Normalize(vmin=CloudLayTf.min(), vmax=CloudLayTf.max())
# # GranulePloter(CloudLayTf, Lat, Lon, plt.cm.gray, norm_CloudLayTf,'default', 'neither', 'linear',\
# #              'Layer top flag', fecha,'LayerTopFlag'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# CloudPPTf = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'PrecipitationFlag')
# norm_CloudPPTf   = colors.Normalize(vmin=CloudPPTf.min(), vmax=CloudPPTf.max())
# GranulePloter(CloudPPTf, Lat, Lon, plt.cm.gray, norm_CloudPPTf,'default', 'neither', 'default',\
#              'Precipitation flag', fecha,'PrecipitationFlag'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# CloudWlayT = DesHDF(path_LIDR+year+'/'+day+'/'+files_LIDR[0], 'Water_layer_top')
# CloudWlayT = np.ma.masked_where(CloudWlayT==-9,CloudWlayT)
# norm_CloudWlayT   = colors.Normalize(vmin=np.nanmin(CloudWlayT), vmax=np.nanmax(CloudWlayT))
# GranulePloter(CloudWlayT, Lat, Lon, plt.cm.gray, norm_CloudWlayT,'default', 'neither', 'default',\
#              'Water layer top', fecha,'WaterLayerTop'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# # ==============================================================================
# #                          Cloud water content Radar Only
# # ==============================================================================
#
# ICEdisWidthpar_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_ice_distrib_width_param')
# ICEdisWidthpar_RO = np.ma.masked_where(ICEdisWidthpar_RO<=0,ICEdisWidthpar_RO)
# norm_Ice_withdis   = colors.Normalize(vmin=np.nanmin(ICEdisWidthpar_RO), vmax=np.nanmax(ICEdisWidthpar_RO))
# GranulePloter(ICEdisWidthpar_RO, Lat, Lon, plt.cm.gray_r, norm_Ice_withdis,'default', 'neither', 'linear',\
#              'Ice width distribution', fecha,'IceWidth_dist_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqdisWidthpar_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_liq_distrib_width_param')
# LiqdisWidthpar_RO = np.ma.masked_where(LiqdisWidthpar_RO<=0,LiqdisWidthpar_RO)
# norm_Liq_withdis   = colors.Normalize(vmin=np.nanmin(LiqdisWidthpar_RO), vmax=np.nanmax(LiqdisWidthpar_RO))
# GranulePloter(LiqdisWidthpar_RO, Lat, Lon, plt.cm.gray_r, norm_Liq_withdis,'default', 'neither', 'linear',\
#              'Liquid width distribution', fecha,'LiquidWidth_dist_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# IceWaterCont_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_ice_water_content')
# IceWaterCont_RO = np.ma.masked_where(IceWaterCont_RO<=0,IceWaterCont_RO)
# norm_IceWat   = colors.Normalize(vmin=np.nanmin(IceWaterCont_RO), vmax=np.nanmax(IceWaterCont_RO))
# GranulePloter(IceWaterCont_RO, Lat, Lon, plt.cm.gray_r, norm_IceWat,'default', 'neither', 'linear',\
#              'Ice water content', fecha,'IceWaterCont_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqWaterCont_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_liq_water_content')
# LiqWaterCont_RO = np.ma.masked_where(LiqWaterCont_RO<=0,LiqWaterCont_RO)
# norm_LiqWat   = colors.Normalize(vmin=np.nanmin(LiqWaterCont_RO), vmax=np.nanmax(LiqWaterCont_RO))
# GranulePloter(LiqWaterCont_RO, Lat, Lon, plt.cm.gray_r, norm_LiqWat,'default', 'neither', 'linear',\
#              'Liquid water content', fecha,'LiqWaterCont_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# IceRadius_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_ice_effective_radius')
# IceRadius_RO = np.ma.masked_where(IceRadius_RO<=0,IceRadius_RO)
# norm_IceRadius   = colors.Normalize(vmin=np.nanmin(IceRadius_RO), vmax=np.nanmax(IceRadius_RO))
# GranulePloter(IceRadius_RO, Lat, Lon, plt.cm.gray_r, norm_IceRadius,'default', 'neither', 'linear',\
#              'Ice effective radius', fecha,'IceRadius_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqRadius_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_liq_effective_radius')
# LiqRadius_RO = np.ma.masked_where(LiqRadius_RO<=0,LiqRadius_RO)
# norm_LiqRadius   = colors.Normalize(vmin=np.nanmin(LiqRadius_RO), vmax=np.nanmax(LiqRadius_RO))
# GranulePloter(LiqRadius_RO, Lat, Lon, plt.cm.gray_r, norm_LiqRadius,'default', 'neither', 'linear',\
#              'Liquid effective radius', fecha,'LiqRadius_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# IceNumber_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_ice_number_conc')
# IceNumber_RO = np.ma.masked_where(IceNumber_RO<=0,IceNumber_RO)
# norm_IceNumber   = colors.Normalize(vmin=np.nanmin(IceNumber_RO), vmax=np.nanmax(IceNumber_RO))
# GranulePloter(IceNumber_RO, Lat, Lon, plt.cm.gray_r, norm_IceNumber,'default', 'neither', 'linear',\
#              'Ice number concentration', fecha,'IceNumber_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqNumber_RO = DesHDF(path_CWCr+year+'/'+day+'/'+files_CWCr[0], 'RO_liq_number_conc')
# LiqNumber_RO = np.ma.masked_where(LiqNumber_RO<=0,LiqNumber_RO)
# norm_LiqNumber   = colors.Normalize(vmin=np.nanmin(LiqNumber_RO), vmax=np.nanmax(LiqNumber_RO))
# GranulePloter(LiqNumber_RO, Lat, Lon, plt.cm.gray_r, norm_LiqNumber,'default', 'neither', 'linear',\
#              'Liquid number concentration', fecha,'LiqNumber_RO'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
#
#
# # ==============================================================================
# #                     Cloud water content Radar vs Optical Depth
# # ==============================================================================
#
#
# ICEdisWidthpar_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_ice_distrib_width_param')
# ICEdisWidthpar_RVOD = np.ma.masked_where(ICEdisWidthpar_RVOD<=0,ICEdisWidthpar_RVOD)
# norm_Ice_withdis   = colors.Normalize(vmin=np.nanmin(ICEdisWidthpar_RVOD), vmax=np.nanmax(ICEdisWidthpar_RVOD))
# GranulePloter(ICEdisWidthpar_RVOD, Lat, Lon, plt.cm.gray_r, norm_Ice_withdis,'default', 'neither', 'linear',\
#              'Ice width distribution', fecha,'IceWidth_dist_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqdisWidthpar_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_ice_distrib_width_param')
# LiqdisWidthpar_RVOD = np.ma.masked_where(LiqdisWidthpar_RVOD<=0,LiqdisWidthpar_RVOD)
# norm_Ice_withdis   = colors.Normalize(vmin=np.nanmin(LiqdisWidthpar_RVOD), vmax=np.nanmax(LiqdisWidthpar_RVOD))
# GranulePloter(LiqdisWidthpar_RVOD, Lat, Lon, plt.cm.gray_r, norm_Ice_withdis,'default', 'neither', 'linear',\
#              'Liquid width distribution', fecha,'LiquidWidth_dist_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# IceWaterCont_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_ice_water_content')
# IceWaterCont_RVOD = np.ma.masked_where(IceWaterCont_RVOD<=0,IceWaterCont_RVOD)
# norm_IceWat   = colors.Normalize(vmin=np.nanmin(IceWaterCont_RVOD), vmax=np.nanmax(IceWaterCont_RVOD))
# GranulePloter(IceWaterCont_RVOD, Lat, Lon, plt.cm.gray_r, norm_IceWat,'default', 'neither', 'linear',\
#              'Ice water content', fecha,'IceWaterCont_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqWaterCont_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_liq_water_content')
# LiqWaterCont_RVOD = np.ma.masked_where(LiqWaterCont_RVOD<=0,LiqWaterCont_RVOD)
# norm_LiqWat   = colors.Normalize(vmin=np.nanmin(LiqWaterCont_RVOD), vmax=np.nanmax(LiqWaterCont_RVOD))
# GranulePloter(LiqWaterCont_RVOD, Lat, Lon, plt.cm.gray_r, norm_LiqWat,'default', 'neither', 'linear',\
#              'Liquid water content', fecha,'LiqWaterCont_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# IceRadius_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_ice_effective_radius')
# IceRadius_RVOD = np.ma.masked_where(IceRadius_RVOD<=0,IceRadius_RVOD)
# norm_IceRadius   = colors.Normalize(vmin=np.nanmin(IceRadius_RVOD), vmax=np.nanmax(IceRadius_RVOD))
# GranulePloter(IceRadius_RVOD, Lat, Lon, plt.cm.gray_r, norm_IceRadius,'default', 'neither', 'linear',\
#              'Ice effective radius', fecha,'IceRadius_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqRadius_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_liq_effective_radius')
# LiqRadius_RVOD = np.ma.masked_where(LiqRadius_RVOD<=0,LiqRadius_RVOD)
# norm_LiqRadius   = colors.Normalize(vmin=np.nanmin(LiqRadius_RVOD), vmax=np.nanmax(LiqRadius_RVOD))
# GranulePloter(LiqRadius_RVOD, Lat, Lon, plt.cm.gray_r, norm_LiqRadius,'default', 'neither', 'linear',\
#              'Liquid effective radius', fecha,'LiqRadius_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# IceNumber_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_ice_number_conc')
# IceNumber_RVOD = np.ma.masked_where(IceNumber_RVOD<=0,IceNumber_RVOD)
# norm_IceNumber   = colors.Normalize(vmin=np.nanmin(IceNumber_RVOD), vmax=np.nanmax(IceNumber_RVOD))
# GranulePloter(IceNumber_RVOD, Lat, Lon, plt.cm.gray_r, norm_IceNumber,'default', 'neither', 'linear',\
#              'Ice number concentration', fecha,'IceNumber_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
# LiqNumber_RVOD = DesHDF(path_CWCo+year+'/'+day+'/'+files_CWCo[0], 'RVOD_liq_number_conc')
# LiqNumber_RVOD = np.ma.masked_where(LiqNumber_RVOD<=0,LiqNumber_RVOD)
# norm_LiqNumber   = colors.Normalize(vmin=np.nanmin(LiqNumber_RVOD), vmax=np.nanmax(LiqNumber_RVOD))
# GranulePloter(LiqNumber_RVOD, Lat, Lon, plt.cm.gray_r, norm_LiqNumber,'default', 'neither', 'linear',\
#              'Liquid number concentration', fecha,'LiqNumber_RVOD'+fecha.strftime('%Y-%m-%d_%H-%M-%S')+'.png')
#
#
#
#
print 'Hello world'
