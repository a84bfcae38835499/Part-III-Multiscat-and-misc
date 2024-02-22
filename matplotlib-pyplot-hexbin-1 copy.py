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
while count < 2:
    # Get next line from file
    line = latticeFile.readline()
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
print(f"Valmin = {valmin}, valmax = {valmax}")
print("===")
print(f"n1min = {n1min}, n1max = {n1max}, n2min = {n2min}, n2max = {n2max}")
print("===")
    

print("Number of occupied channels = " + str(nOccCh))
H = calculate_entropy(plotValues)
print("Entropy = " + str(H))

n1Span = n1max-n1min
xSpan = n1Span*0.5
n2Span = n2max-n2min
print(xSpan)
ySpan = int(1+np.sqrt(3)*(xSpan))
if(ySpan%2==0):
    ySpan = ySpan - 1

fig = plt.figure()
#h = ax.hexbin(x, y, gridsize=(int(np.sqrt(3)*scalefact), int(scalefact)))
#print("x coords = ")
#print(plotCoordsX)
#print("y coords = ")
#print(plotCoordsY)
# from curved coordinate to rectlinear coordinate.
fudge = 1 #I have absolutely no idea why this is needed
def tr(x, y):
    matTrans = np.matrix([[ 1/2, np.sqrt(3)/2],[0,1]])
    x, y = np.asarray(x), np.asarray(y)
    return x * matTrans[0,0] + y * matTrans[0,1],x * matTrans[1,0] + y * matTrans[1,1]

# from rectlinear coordinate to curved coordinate.
def inv_tr(x, y):
    matTrans = np.matrix([[ 1/2, np.sqrt(3)/2],[0,1]])
    invMat = np.linalg.inv(matTrans)
    x, y = np.asarray(x), np.asarray(y)
    return x * invMat[0,0] + y * invMat[0,1],x * invMat[1,0] + y * invMat[1,1]

ax = fig.add_subplot(projection='skewx')
# Add the grid
plt.grid(True)
ax.grid(which='major', axis='both', linestyle='--',color=[0., 0., 0.])

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
# An example of a slanted line at constant X
l = ax.axvline(0, color='C0')

# Disables the log-formatting that comes with semilogy
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    
useLog = True
if(useLog):
    h = ax.hexbin(plotCoordsX,plotCoordsY,C=plotValues,gridsize=(ySpan, int(xSpan)),cmap='magma',norm=mpl.colors.LogNorm(valmin,valmax))
    #plt.gca().set_aspect('equal')
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.LogNorm(vmin=valmin,vmax=valmax), cmap='magma'),
             ax=ax, orientation='vertical', label='P($n_1$,$n_2$)')
else:
    h = ax.hexbin(plotCoordsX,plotCoordsY,C=plotValues,gridsize=(ySpan, int(xSpan)),cmap='magma',norm=mpl.colors.Normalize(valmin,valmax))
    #plt.gca().set_aspect('equal')
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=valmin,vmax=valmax), cmap='magma'),
             ax=ax, orientation='vertical', label='P($n_1$,$n_2$)')

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

ax.set_title(titelstr)

captiontxt="Entropy = " + "{:.6f}".format(H)
plt.figtext(0.5, 0.04, captiontxt, wrap=True, horizontalalignment='center', fontsize=12)

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

filenametxt=""
filenametxt=""
plt.figtext(0.5, 0.005, filenametxt, wrap=True, horizontalalignment='center', fontsize=12,fontstyle='italic')

if(filenametxt == ""):
    savestr = "Figures/Diffraction/" + datetime.datetime.now().strftime('Diffraction_%Y-%m-%d_%H-%M') + "_Hex.png"
else:
    savestr = "Figures/Diffraction/" +slugify(filenametxt) +datetime.datetime.now().strftime('_%Y-%m-%d_%H-%M') + ".png"
print(savestr)
plt.savefig(fname=savestr)
plt.show()