import numpy as np
import matplotlib.pyplot as plt

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

for k in range(0,nOccCh):
    row = d.iloc[k]
    n1 = getattr(row,'n1')
    n2 = getattr(row,'n2')
    I = getattr(row,'I')
    print(f"k = {k}, n1 = {n1}, n2 = {n2}, I = {I}")
    plotValues[k] = I
    plotCoordsX[k] = B1[0] * n1 + B2[0] * n2
    plotCoordsY[k] = B1[1] * n1 + B2[1] * n2

print("Number of occupied channels = " + str(nOccCh))

scalefact = 5
fig, ax = plt.subplots(figsize=(4, 4))
#h = ax.hexbin(x, y, gridsize=(int(np.sqrt(3)*scalefact), int(scalefact)))
print("x coords = ")
print(plotCoordsX)
print("y coords = ")
print(plotCoordsY)
h = ax.hexbin(plotCoordsX/np.sqrt(B1[0]**2+B1[1]**2),plotCoordsY/np.sqrt(B1[0]**2+B1[1]**2),C=plotValues,gridsize=(int(np.sqrt(3)*scalefact), int(scalefact)))
hx, hy = h.get_offsets().T
plt.show()