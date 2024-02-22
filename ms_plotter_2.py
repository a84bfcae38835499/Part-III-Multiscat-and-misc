import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import mpl_toolkits.axisartist as AA
import mpl_toolkits.axisartist.grid_finder as gf
import numpy as np
import pandas as pd
import datetime
import skewaxes
# Default theme
#mpl.rc('axes',edgecolor='white')
import unicodedata
import re

def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


def import_multiscat(fname):    
    """Impot standrad multiscat output into a pandas data frame."""
    d = pd.read_csv(fname, skiprows=7, delim_whitespace=True, 
                    header=None, names=['#','n1','n2','I'])
    d.drop(columns=['#'], inplace=True)
    return(d)


def calculate_entropy(intensities):
    H = 0
    for p in np.nditer(intensities):
        if(np.isnan(p)):
            print("NaN found! skipping...")
        elif(p == 0):
            print("Zero found! skpping...")
        else:
        #    print("p = " + str(p))
            H += p * np.log(p)
    H /= -np.log(intensities.size)
    return(H)

latticeFile = open('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script', 'r')
count = 0 

B1 = [0,0]
B2 = [0,0]
a1 = [0,0]
a2 = [0,0]
Nsuper = 1337
while count < 5:
    # Get next line from file
    line = latticeFile.readline()
    if(line.startswith("Nsuper = ")):
        line = line[len("Nsuper = "):]
        split = line.split()
        Nsuper = float(split[0])
        print("Nsuper = " + split[0])
        count += 1

    if line.startswith("a1 = "):
        line = line[len("a1 = "):]
        print("line = " + line)
        split = line.split()

        a1[0] = float(split[0])
        a1[1] = float(split[1])
        count += 1
 
    if line.startswith("a2 = "):
        line = line[len("a2 = "):]
        print("line = " + line)
        split = line.split()
        a2[0] = float(split[0])
        a2[1] = float(split[1])
        count += 1

    if line.startswith("b1 = "):
        line = line[len("b1 = "):]
        print("line = " + line)
        split = line.split()

        B1[0] = float(split[0])
        B1[1] = float(split[1])
        count += 1
 
    if line.startswith("b2 = "):
        line = line[len("b2 = "):]
        print("line = " + line)
        split = line.split()
        B2[0] = float(split[0])
        B2[1] = float(split[1])
        count += 1

Babs = np.sqrt(B1[0]**2+B1[1]**2)
print("B1 = ") 
print(B1)
print("B2 = ") 
print(B2)
print("Babs = ") 
print(Babs)
latticeFile.close()
b1 = B1 / Babs
b2 = B2 / Babs

print("b1 = ") 
print(b1)
print("b2 = ") 
print(b2)
print("Babs = ")

d = import_multiscat('diffrac10001.out')
print(d)

#print(d2)
#print(d2.values)
nOccCh = len(d.index)
plotValues = np.zeros((nOccCh))
plotCoordsX = np.zeros((nOccCh))
plotCoordsY = np.zeros((nOccCh))
n1min = 0
n1max = 0
n2min = 0
n2max = 0
valmin = 1
valmax = 0
smolVal = 1e-10
for k in range(0,nOccCh):
    row = d.iloc[k]
    n1 = getattr(row,'n1')
    if(n1 < n1min):
        n1min = n1
    if(n1 > n1max):
        n1max = n1
    n2 = getattr(row,'n2')
    if(n2 < n2min):
        n2min = n2
    if(n2 > n2max):
        n2max = n2
    I = getattr(row,'I')
    if(I == 0):
        print(f"Zero found, setting to {smolVal}")
        I = smolVal
    if(I < valmin):
        valmin = I
    if(I > valmax):
        valmax = I
    #print(f"k = {k}, n1 = {n1}, n2 = {n2}, I = {I}")
    plotValues[k] = I
    plotCoordsX[k] = b1[0] * n1 + b2[0] * n2
    plotCoordsY[k] = b1[1] * n1 + b2[1] * n2
plotCoordsArray = np.array([])
plotCoordsArray = np.array(np.column_stack((plotCoordsX, plotCoordsY)))
print(plotCoordsArray.shape)
print(f"Valmin = {valmin}, valmax = {valmax}")
print("===")
print(f"n1min = {n1min}, n1max = {n1max}, n2min = {n2min}, n2max = {n2max}")
print("===")
    

print("Number of occupied channels = " + str(nOccCh))
H = calculate_entropy(plotValues)
print("Entropy = " + str(H))

#packages to import
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

#sets the colour scale
useLog = False
if(useLog):
    mapper = cm.ScalarMappable(cmap='magma', norm=mpl.colors.LogNorm(valmin,valmax))
else:
    mapper = cm.ScalarMappable(cmap='magma', norm=mpl.colors.Normalize(valmin,valmax))

#create the figure with set figure size
fig = plt.figure(figsize=(10,8))

#puts a grid on the figure
gs = gridspec.GridSpec(16, 10) 

#creates two subplots
ax = plt.subplot2grid((16,20), (0,17), colspan=1, rowspan=16)
ax2 = plt.subplot2grid((16,20), (0,0), colspan=16, rowspan=16)

#creates a colourbar on the first subplot
if(useLog):
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap='magma', norm=mpl.colors.LogNorm(valmin,valmax), orientation='vertical')
else:
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap='magma', norm=mpl.colors.Normalize(valmin,valmax), orientation='vertical')
cb1.set_label('P($n_1$,$n_2$)')

vor = Voronoi(points=plotCoordsArray,furthest_site=False)

#plots the voronoi diagram on the second subplot
voronoi_plot_2d(vor, show_vertices =False, show_points =False, ax=ax2)

#gives the order values for colour scale (which come from a data frame in a different part of the code)
order = plotValues

#colours the voronoi cells    
for r in range(len(vor.point_region)):
    region = vor.regions[vor.point_region[r]]
    if not -1 in region:
        polygon = [vor.vertices[i] for i in region]
        plt.fill(*zip(*polygon), color=mapper.to_rgba(order[r]))
        #plt.fill(*zip(*polygon))


import matplotlib.patheffects as pe
pathefts1 = [pe.Stroke(linewidth=1, foreground='w'), pe.Normal()]
pathefts2 = [pe.Stroke(linewidth=2, foreground='w'), pe.Normal()]

plt.arrow(0,0,b1[0]*Nsuper*0.7,b1[1]*Nsuper*0.7,width=0.1,color='r',zorder=7,path_effects=pathefts2)
plt.arrow(0,0,b2[0]*Nsuper*0.7,b2[1]*Nsuper*0.7,width=0.1,color='g',zorder=7,path_effects=pathefts2)
plt.annotate("b1", (b1[0]*Nsuper,b1[1]*Nsuper),color='r',fontsize=12,weight='bold',path_effects=pathefts1,zorder=8)
plt.annotate("b2", (b2[0]*Nsuper,b2[1]*Nsuper),color='g',fontsize=12,weight='bold',path_effects=pathefts1,zorder=8)

a1col = [1, 0.5, 0.6]
a2col = [0.2, 0.7, 0.3]
plt.arrow(0,0,a1[0],a1[1],width=0.1,color=a1col,zorder=6,path_effects=pathefts2)
plt.arrow(0,0,a2[0],a2[1],width=0.1,color=a2col,zorder=6,path_effects=pathefts2)
plt.annotate("a1", (a1[0],a1[1]+0.2),color=a1col,fontsize=12,weight='bold',path_effects=pathefts1)
plt.annotate("a2", (a2[0],a2[1]+0.5),color=a2col,fontsize=12,weight='bold',path_effects=pathefts1)

for k in range(0,nOccCh):
    row = d.iloc[k]
    n1 = int(getattr(row,'n1'))
    n2 = int(getattr(row,'n2'))
    n1n2 = str(n1) + ',' + str(n2)
    plt.annotate(n1n2,((b1[0]*n1+b2[0]*n2)*Nsuper-0.3,(b1[1]*n1+b2[1]*n2)*Nsuper-0.12),fontsize=6,zorder=10)
ax2.set_aspect('equal')

scatFile = open('scatCond.in', 'r')
E = -69
theta = -69
phi = -69
line = scatFile.readline()
line = scatFile.readline()
vals = line.split(",")
E = float(vals[0])
theta = float(vals[1])
phi = float(vals[2])
titelstr = "$E$ = " + str(E) + " meV, $\\theta$ = " + str(theta) + "$\\degree$, $\\phi$ =" + str(phi) + "$\\degree$"
print(titelstr)
scatFile.close()

ax2.set_title(titelstr)

captiontxt="Entropy = " + "{:.6f}".format(H)
plt.figtext(0.5, -0.07, captiontxt, wrap=True, horizontalalignment='center', fontsize=12,transform=ax2.transAxes)
filenametxt=""
filenametxt="Italicised comment here"
plt.figtext(0.5, -0.11, filenametxt, wrap=True, horizontalalignment='center', fontsize=12,fontstyle='italic',transform=ax2.transAxes)

if(filenametxt == ""):
    savestr = "Figures/Diffraction/" + datetime.datetime.now().strftime('Diffraction_%Y-%m-%d_%H-%M') + "_Hex.png"
else:
    savestr = "Figures/Diffraction/" +slugify(filenametxt) +datetime.datetime.now().strftime('_%Y-%m-%d_%H-%M') + ".png"
print(savestr)
plt.savefig(fname=savestr,dpi=300)
plt.show()