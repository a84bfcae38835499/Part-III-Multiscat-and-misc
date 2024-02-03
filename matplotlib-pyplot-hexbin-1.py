import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import datetime
# Default theme

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
    if line.startswith("B1 = "):
        line = line[len("B1 = "):]
        print("line = " + line)
        split = line.split()

        B1[0] = float(split[0])
        B1[1] = float(split[1])
        count += 1
 
    if line.startswith("B2 = "):
        line = line[len("B2 = "):]
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

d = import_multiscat('diffrac10001.out')
print(d)

#print(d2)
#print(d2.values)
nOccCh = len(d.index)
plotValues = np.zeros((nOccCh))
plotCoordsX = np.zeros((nOccCh))
plotCoordsY = np.zeros((nOccCh))
nmin = 0
nmax = 0
for k in range(0,nOccCh):
    row = d.iloc[k]
    n1 = getattr(row,'n1')
    if(n1 < nmin):
        nmin = n1
    if(n1 > nmax):
        nmax = n1
    n2 = getattr(row,'n2')
    I = getattr(row,'I')
    #print(f"k = {k}, n1 = {n1}, n2 = {n2}, I = {I}")
    plotValues[k] = I
    plotCoordsY[k] = -B1[0] * n1 - B2[0] * n2 #for some reason everything gets inverted?? and x and y are swapped from what I'd expect?????
    plotCoordsX[k] = -B1[1] * n1 - B2[1] * n2

print("Number of occupied channels = " + str(nOccCh))
H = calculate_entropy(plotValues)
print("Entropy = " + str(H))

fig, ax = plt.subplots(figsize=(4, 4))
#h = ax.hexbin(x, y, gridsize=(int(np.sqrt(3)*scalefact), int(scalefact)))
#print("x coords = ")
#print(plotCoordsX)
#print("y coords = ")
#print(plotCoordsY)

xSpan = nmax-nmin
xSpan = xSpan
print(xSpan)

h = ax.hexbin(plotCoordsX/Babs,plotCoordsY/Babs,C=plotValues,gridsize=(int(np.sqrt(3)*(xSpan)), int(xSpan)),cmap='magma',bins='log',norm=mpl.colors.Normalize(0,1))
 
fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.LogNorm(vmin=1e-5,vmax=1), cmap='magma'),
             ax=ax, orientation='vertical', label='P($n_1$,$n_2$)')

savestr = "Figures/Diffraction/" + datetime.datetime.now().strftime('Diffraction_%Y-%m-%d_%H-%M') + "_Hex.png"
plt.savefig(fname=savestr)
plt.show()