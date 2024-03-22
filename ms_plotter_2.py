import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import datetime
import unicodedata
import re
import matplotlib.patheffects as pe
import matplotlib.patches as patches
import math

#packages to import
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
import matplotlib.cm as cm

fileprefix = '3x3ikbt_04'
fileprefix = 'restest_16_50'
fileprefix = '3x3highdefect_adatom'
fileprefix = '5x5MoS2'
fileprefix = '1x1pristine'
fileprefix = '1x1unpristine'
plotFigure = False
useLog = False

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
    d = pd.read_csv(fname, skiprows=1, delim_whitespace=True, 
                    header=None, names=['#','n1','n2','I'])
    d.drop(columns=['#'], inplace=True)
    nlines = sum(1 for line in open(fname))
    topOfSeperation = 0
    dfs = []
    print("nlines = " + str(nlines))
    for index in range(0,nlines-2):
        val1 = d.loc[index,'n1']
        val2 = d.loc[index+1,'n1']
        difference = val2 - val1
        #print("val2 - val1 = "+ str(val2 - val1))
        if(difference < 0):
            print("\nSeperation found! Splitting dataframe between " + str(topOfSeperation) + ", " + str(index))
            dslice = d.iloc[topOfSeperation:index+1,:]
            topOfSeperation = index+1
            dfs.append(dslice)
    print("===")
    dslice = d.iloc[topOfSeperation:,:]
    dfs.append(dslice)
    print("\nFinshed!")
    return(dfs)

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

def find_mean_stdv(values):
    mean = 0
    meanSq = 0
    #print("input = ")
    #print(values)
    #print("type of input = " + str(type(values)))
    #n = len(values)
    n = np.size(values)
    #print("n = " + str(n))
    for v in np.nditer(values):
        mean += v
        meanSq += v**2
    mean /= n
    meanSq /= n
    #print(f"Mean**2, msq = {mean**2}, {meanSq}")
    stdv = np.sqrt(meanSq - mean**2)
    return mean, stdv


latticeFile = open(fileprefix + '.info_for_vivian_python_nice_plotting_hexagon_script', 'r')
count = 0 

B1 = np.array([0.,0.])
B2 = np.array([0.,0.])
a1 = np.array([0.,0.])
a2 = np.array([0.,0.])
Nsuper = 1337
while count < 10:
    # Get next line from file
    line = latticeFile.readline()
    if(line.startswith("Nsuper = ")):
        line = line[len("Nsuper = "):]
        split = line.split()
        Nsuper = float(split[0])
        #print("Nsuper = " + split[0])
        count += 1
    
    if(line.startswith("Unit cell area = ")):
        line = line[len("Unit cell area = "):]
        split = line.split()
        cellArea = float(split[0])
        #print("Nsuper = " + split[0])
        count += 1

    if line.startswith("a1 = "):
        line = line[len("a1 = "):]
        #print("line = " + line)
        split = line.split()

        a1[0] = float(split[0])
        a1[1] = float(split[1])
        count += 1
 
    if line.startswith("a2 = "):
        line = line[len("a2 = "):]
        #print("line = " + line)
        split = line.split()
        a2[0] = float(split[0])
        a2[1] = float(split[1])
        count += 1

    if line.startswith("b1 = "):
        line = line[len("b1 = "):]
        #print("line = " + line)
        split = line.split()

        B1[0] = float(split[0])
        B1[1] = float(split[1])
        count += 1
 
    if line.startswith("b2 = "):
        line = line[len("b2 = "):]
        #print("line = " + line)
        split = line.split()
        B2[0] = float(split[0])
        B2[1] = float(split[1])
        count += 1
    if line.startswith("Theta = "):
        line = line[len("Theta = "):]
        split = line.split()
        Theta = float(split[0])
        #print("Theta = " + split[0])
        count += 1
    if line.startswith("Ensemble size = "):
        line = line[len("Ensemble size = "):]
        split = line.split()
        Nensemble = int(split[0])
        #print("Ensemble size = " + split[0])
        count += 1
    if line.startswith("Positional entropy = "):
        line = line[len("Positional entropy = "):]
        split = line.split()
        entropyInstr = split[0]
        print("Positional entropy = " + split[0])
        count += 1
    if line.startswith("Defect density in cm^-2 = "):
        line = line[len("Defect density in cm^-2 = "):]
        split = line.split()
        defectDensity = float(split[0])
        #print("Defect density in cm^-2 = " + split[0])
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

scatFile = open('scatCond.in', 'r')

Nscat = sum(1 for _ in scatFile)
Nscat -= 1

Es = []
thetas = []
phis = []

scatFile.close()
scatFile = open('scatCond.in', 'r')

line = scatFile.readline()
for index_s in range(Nscat):
    line = scatFile.readline()
    print("Line = " + line)
    vals = line.split(",")
    Es.append( float(vals[0]) )
    thetas.append( float(vals[1]) )
    phis.append( float(vals[2]) )

scatFile.close()

print("[][][][][][][][]")
print("Number of scattering conditions = " + str(Nscat))
print("[][][][][][][][]\n")

dfss = []

for index_n in range(Nensemble):
    importname =  fileprefix + str(10001+index_n) + '.out'
    print("importing file : " + importname)
    dfs = import_multiscat(importname)
    #dfs has scattering varying scattering conditions for one potential
    #print("===")
    #print(dfs)
    dfss.append(dfs)

nOccChArr = [0]*Nscat

for index_s in range(Nscat):
    print("index_s = " + str(index_s))
    rubbis = dfss[0]
    trialdf = rubbis[index_s]
    nOccChArr[index_s] = len(trialdf.index)
    #Assume that all ensemble members for a particular scatering cond
    #have the same nubmer of channels because otherwise I will actually go insane

entropiesOut = [0.]*Nscat
entropiesOutUnc = [0.]*Nscat
intensityArr = [None]*Nscat
coordXArr = [None]*Nscat
coordYArr = [None]*Nscat
brightSpotXArr = [None]*Nscat
brightSpotYArr = [None]*Nscat
kAvg = np.array([0.,0.])
kAbsAvgArr = [0.]*Nscat
kAbsAvgUncArr = [0.]*Nscat

n1Arr = [None]*Nscat
n2Arr = [None]*Nscat
n1minArr = [0]*Nscat
n1maxArr = [0]*Nscat
n2minArr = [0]*Nscat
n2maxArr = [0]*Nscat

valminArr = [1.]*Nscat
valmaxArr = [0.]*Nscat
smolVal = 1e-100
vanityVal = 0
vanity = False

for index_s in range(Nscat):
    print("Now processing numscatcond = " + str(index_s))
    print("Trail df length = " + str(nOccChArr[index_s]))
    intensityArr[index_s] =(np.zeros((nOccChArr[index_s])))
    coordXArr[index_s] =(np.zeros((nOccChArr[index_s]))) #this is what paranoia looks like
    coordYArr[index_s] =(np.zeros((nOccChArr[index_s])))
    cX = [None]*Nscat
    cY = [None]*Nscat
    iI = [None]*Nscat
    n1Arr[index_s] =(np.zeros((nOccChArr[index_s]),dtype='int'))
    n2Arr[index_s] =(np.zeros((nOccChArr[index_s]),dtype='int'))

    brightSpotXArr[index_s] =(np.zeros((nOccChArr[index_s])))
    brightSpotYArr[index_s] =(np.zeros((nOccChArr[index_s])))

    for index_n in range(Nensemble):
        df = dfss[index_n][index_s]
        intensities = np.zeros((nOccChArr[index_s]))
        plotCoordsX = np.zeros((nOccChArr[index_s]))
        plotCoordsY = np.zeros((nOccChArr[index_s]))
        for ch in range(nOccChArr[index_s]):
            #print("k = " + str(k))
            row = df.iloc[ch]
            n1 = getattr(row,'n1')
            n1Arr[index_s][ch] = int(n1)
            n2 = getattr(row,'n2')
            n2Arr[index_s][ch] = int(n2)
            I = float(getattr(row,'I'))
            #print(f"k = {k}, n1 = {n1}, n2 = {n2}, I = {I}")
            intensities[ch] = I
            plotCoordsX[ch] = b1[0] * n1 + b2[0] * n2
            plotCoordsY[ch] = b1[1] * n1 + b2[1] * n2
        cX[index_n] = (plotCoordsX)
        cY[index_n] = (plotCoordsY)
        iI[index_n] = (intensities)
    #Now do ensemble averageing, as well as calculating of meaningful quantities so they can have error bars too
    
    iAverage = np.zeros((nOccChArr[index_s]))
    entropiesDisposable = np.zeros((Nensemble))
    kAbsDisposable = np.zeros((Nensemble))
    for index_n in range(Nensemble):
        Is = iI[index_n]
        for ch in range(nOccChArr[index_s]):
            row = df.iloc[ch]
            I = Is[ch]
            n1 = float(getattr(row,'n1'))
            n2 = float(getattr(row,'n2'))
            #print(f"n1, n2, I = {n1}, {n2}, {I}")
            kAvg += I*(n1*B1 + n2*B2)/float(Nsuper)
            #print("kAvg = " + str(kAvg))
            iAverage[ch] += I
        entropiesDisposable[index_n] = calculate_entropy(iI[index_n])
        kAbsDisposable[index_n] = np.sqrt(kAvg[0]**2+kAvg[1]**2)
    print("kAbsDisposable = ")
    print(kAbsDisposable)
    kAvg /= float(Nensemble)
    kAvg /= Babs
    print("kAvg (n1n2 units)= ")
    print(kAvg)
    iAverage /= float(Nensemble)
    intensityArr[index_s] = iAverage
    coordXArr[index_s] = cX[0] #this is what paranoia looks like
    coordYArr[index_s] = cY[0]
    entropiesOut[index_s], entropiesOutUnc[index_s] = find_mean_stdv(entropiesDisposable)
    kAbsAvgArr[index_s], kAbsAvgUncArr[index_s] = find_mean_stdv(kAbsDisposable)
    print("kAbsAvg = ")
    print(kAbsAvgArr[index_s])
    print("Unc = ")
    print(kAbsAvgUncArr[index_s])

    #Finally, for plotting, we remove channels of 0 intensity cos they cause issues when plotting
    for index_n in range(Nensemble):
        Is = iI[index_n]
        for ch in range(nOccChArr[index_s]):
            I = Is[ch]
            if(I == 0):
                print(f"Zero found, setting to {smolVal}")
                I = smolVal

#We now have our varrious arrays whose index ranges over the number of scattering conditions
#(these arrays are labelled xxxArr)

#Now we find the min and max value spots, as well as spots brigth enough to be considered plotable (I > vanityVal)
for index_s in range(Nscat):
    pCXS = np.array([]) #These variables stand for something but icr what it is lmao
    pCYS = np.array([])
    Is = intensityArr[index_s]
    for ch in range(nOccChArr[index_s]):
            I = Is[ch]
            n1 = n1Arr[index_s][ch]
            n2 = n2Arr[index_s][ch]
            if(I == 0):
                print(f"Zero found, setting to {smolVal}")
                I = smolVal

            if(I < valminArr[index_s]):
                valminArr[index_s] = I
            if(I > valmaxArr[index_s]):
                valmaxArr[index_s] = I
            if(I > vanityVal):
                pCXS = np.append(pCXS,b1[0] * float(n1) + b2[0] * float(n2))
                pCYS = np.append(pCYS,b1[1] * float(n1) + b2[1] * float(n2))
    brightSpotXArr[index_s] = pCXS
    brightSpotYArr[index_s] = pCYS

#This time is used for file naming, set it now so if you spend ages looking
#at a plot and the timer ticks over, the files after still have the same
#name as the ones before
execTime = datetime.datetime.now().strftime('_%Y-%m-%d_%H-%M')

for index_s in range(Nscat):
    plotValuesAvg = intensityArr[index_s]
    plotCoordsX = coordXArr[index_s]
    plotCoordsY = coordYArr[index_s]
    #sets the colour scale
    if(useLog):
        mapper = cm.ScalarMappable(cmap='magma', norm=mpl.colors.LogNorm(valminArr[index_s],valmaxArr[index_s]))
    else:
        mapper = cm.ScalarMappable(cmap='magma', norm=mpl.colors.Normalize(valminArr[index_s],valmaxArr[index_s]))

    #create the figure with set figure size
    fig = plt.figure(figsize=(10,8))

    #creates two subplots
    ax = plt.subplot2grid((16,20), (0,17), colspan=1, rowspan=16)
    ax2 = plt.subplot2grid((16,20), (0,0), colspan=16, rowspan=16)
    #print(thetas)
    E = Es[index_s]
    theta = thetas[index_s]
    #print(theta)
    phi = phis[index_s]

    
    titelstr = "$E$ = " + str(E) + " meV, $\\theta$ = " + str(theta) + "$\\degree$, $\\phi$ =" + str(phi) + "$\\degree$"
    scatcondstr = str(E) + "_" + str(theta) + "_" + str(phi)
    #print(scatcondstr)

    pathefts1 = [pe.Stroke(linewidth=1, foreground='w'), pe.Normal()]
    pathefts2 = [pe.Stroke(linewidth=2, foreground='w'), pe.Normal()]
    b1col = [0.9, 0.1, 0]
    b2col = [0, 0.7, 0.2]
    a1col = [1, 0.5, 0.6]
    a2col = [0.5, 0.8, 0.6]
    hecol = [0, 0.3, 0.8]

    normDiffI = 0.
    normSpecI = 0.
    nDiffCh = 0
    nSpecCh = 0

    if(not vanity):
        plt.arrow(0,0,b1[0]*Nsuper,b1[1]*Nsuper,width=0.05,color=b1col,zorder=7,path_effects=pathefts2,length_includes_head=True)
        plt.arrow(0,0,b2[0]*Nsuper,b2[1]*Nsuper,width=0.05,color=b2col,zorder=7,path_effects=pathefts2,length_includes_head=True)
        plt.annotate("b1", (b1[0]*Nsuper,b1[1]*Nsuper+0.1),color=b1col,fontsize=8,weight='bold',path_effects=pathefts1,zorder=11)
        plt.annotate("b2", (b2[0]*Nsuper,b2[1]*Nsuper),color=b2col,fontsize=8,weight='bold',path_effects=pathefts1,zorder=11)

        plt.arrow(0,0,a1[0]/np.sqrt(a1[0]**2+a1[1]**2),a1[1]/np.sqrt(a1[0]**2+a1[1]**2),width=0.05,color=a1col,zorder=6,length_includes_head=True,alpha=0.5)
        plt.arrow(0,0,a2[0]/np.sqrt(a1[0]**2+a1[1]**2),a2[1]/np.sqrt(a1[0]**2+a1[1]**2),width=0.05,color=a2col,zorder=6,length_includes_head=True,alpha=0.5)
        plt.annotate("a1", (a1[0]/np.sqrt(a1[0]**2+a1[1]**2),a1[1]/np.sqrt(a1[0]**2+a1[1]**2)),color=a1col,fontsize=8,weight='bold',zorder = 5)
        plt.annotate("a2", (a2[0]/np.sqrt(a1[0]**2+a1[1]**2),a2[1]/np.sqrt(a1[0]**2+a1[1]**2)),color=a2col,fontsize=8,weight='bold',zorder = 5)
    mean1 = 0.
    mean2 = 0.
    weighted1 = 0.
    weighted2 = 0.
    kWeighted = [0.,0.]
    for ch in range(nOccChArr[index_s]):
        n1 = n1Arr[index_s][ch]
        n2 = n2Arr[index_s][ch]
        I = plotValuesAvg[ch]
        mean1 += n1
        mean2 += n2
        weighted1 += n1 * I
        weighted2 += n2 * I
        if(n1%int(Nsuper) == 0 and n2%int(Nsuper)==0):
            n1n2 = str(int(n1/Nsuper)) + ',' + str(int(n2/Nsuper))
            #print("Annotating point " + n1n2)
            if(plotValuesAvg[ch] < (valmaxArr[index_s]-valminArr[index_s])*0.75):
                col = 'w'
            else:
                col = 'k'
            if(not vanity):
                plt.annotate(n1n2,((b1[0]*float(n1)+b2[0]*float(n2)),
                                    (b1[1]*float(n1)+b2[1]*float(n2))),
                            fontsize=4,zorder=12,ha='center',va='center',c=col)
            normSpecI += plotValuesAvg[ch]
            nSpecCh += 1
        else:
            normDiffI += plotValuesAvg[ch]
            nDiffCh += 1

    weighted1 /= nOccChArr[index_s]
    weighted2 /= nOccChArr[index_s]
    weightedX = weighted1*B1[0]+weighted2*B2[0]
    weightedY = weighted1*B1[1]+weighted2*B2[1]
    print("weightedX, weightedY = ")
    print(weightedX, weightedY)
    if(Nensemble == 1):
        kstr = "$|K|$ = " + "{:.3f}".format(kAbsAvgArr[index_s]) + "Å$^{-1}$"
    else:
        kstr = "$|K|$ = " + "{:.3f}".format(kAbsAvgArr[index_s]) + "$\pm$" +  "{:.3f}".format(kAbsAvgUncArr[index_s])+ "Å$^{-1}$"
    print(kstr)
    if(not vanity):
        if(not(math.isclose(theta,0.) & math.isclose(phi,0.))):
            #plt.arrow(0,0,
            #          kAvg[0],kAvg[1],
            #          width=0.01,color='w',zorder=7,head_width=0.5,head_length=0.5,length_includes_head=True,
            #         fill=False,overhang=-1.)
                     #fill=True)
            ax2.arrow(0., 0.,kAvg[0],kAvg[1], zorder=8 ,linestyle=(0,(1, 10)), color='w',linewidth=.75,head_length=0.0,head_width=0.0)
            ax2.scatter(kAvg[0],kAvg[1], marker = '+',zorder=8,color='w',linewidth=1,s=5e2)
    mean1 /= nOccChArr[index_s]
    mean2 /= nOccChArr[index_s]
    meanX = mean1*b1[0]+mean2*b2[0]
    meanY = mean1*b1[1]+mean2*b2[1]
    print("mean x, y = ")
    print(meanX, meanY)
    #normSpecI *= float(nSpecCh+nDiffCh)/float(nSpecCh)
    #normDiffI *= float(nSpecCh+nDiffCh)/float(nDiffCh)
    #print("Number of specular channels : " + str(nSpecCh))
    #print("Number of diffuse (non-diffractive) channels : " + str(nDiffCh))
    #print("Specular intensity proportion : " + str(normSpecI))
    #print("Diffuse intensity proportion  : " + str(normDiffI))

    heliumRot = np.matrix([[np.cos(np.deg2rad(phi)),np.sin(np.deg2rad(phi))],
                        [-np.sin(np.deg2rad(phi)),np.cos(np.deg2rad(phi))]])
    heliumDir = -heliumRot * np.reshape(np.array(a1),(2,1))/np.sqrt(np.dot(a1,a1))
    e = 1.602176634E-19
    hbar = 6.62607015e-34/(2*np.pi)
    m = 4 * 1.66053906660e-27
    #heliumk = np.sin(np.deg2rad(theta))*heliumDir * np.sqrt(2*m*e*E/1000)/hbar
    heliumk = heliumDir * np.sqrt(2*m*e*E/1000)/hbar
    heliumk_n = heliumk/(Babs*1e10)
    if(not vanity):
        #print("heliumk_n =")
        #print(heliumk_n)
        if(not(math.isclose(theta,0.) & math.isclose(phi,0.))):
            arrowl = 1
            plt.arrow(meanX+Nsuper*heliumk_n[0,0]-arrowl*np.cos(np.deg2rad(phi)),meanY+Nsuper*heliumk_n[1,0]+arrowl*np.sin(np.deg2rad(phi)),
            #plt.arrow(0,0,
                      arrowl*np.cos(np.deg2rad(phi)),-arrowl*np.sin(np.deg2rad(phi)),
                      width=0.00,color='b',zorder=7,head_width=1,head_length=arrowl,length_includes_head=True,
                     fill=False,overhang=1.)
        ax2.add_patch(plt.Circle((meanX,meanY), np.sqrt(heliumk_n[0]**2 + heliumk_n[1]**2)*Nsuper, color='b', fill=False,zorder=7,linestyle=(0, (25, 50))))
    ax2.set_title(titelstr)

    #creates a colourbar on the first subplot
    #print("valminArr[index_s],valmaxArr[index_s = " + str(valminArr[index_s]) + "," + str(valmaxArr[index_s]))
    if(useLog):
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap='magma', norm=mpl.colors.LogNorm(valminArr[index_s],valmaxArr[index_s]), orientation='vertical')
    else:
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap='magma', norm=mpl.colors.Normalize(valminArr[index_s],valmaxArr[index_s]), orientation='vertical')
    cb1.set_label('P($n_1$,$n_2$)')

    
    eomean = entropiesOut[index_s]
    eosdtv = entropiesOutUnc[index_s]
    defectstr="$n_{defect}$ = " + "{:.4e}".format(defectDensity) + " cm$^{-2}$, $H_{defect}$ = " + entropyInstr
    if(Nensemble == 1):
        entropytxt = "$H_{diffraction}$ = " + "{:.6f}".format(eomean)
    else:
        entropytxt = "$H_{diffraction}$ = " + "{:.6f}".format(eomean) + "$\pm$" +  "{:.6f}".format(eosdtv)
    
    crossSectionWhole = cellArea * np.log(normSpecI)/np.log(1-Theta)
    simgastr = "$\Sigma_T = $" + "{:.2f}".format(crossSectionWhole) + "Å$^2$"
    intenstr = "Proportion of specular intensity = " + str(int(normSpecI*100))+ "%"
    print(simgastr)

    additionalX = []
    additionalY = []
    additionalVals = []

    padCells = True
    if(padCells):
        paddingCells = 20
        for n in range(-paddingCells,paddingCells):
            for m in range(-paddingCells,paddingCells):
                canPlaceSiteHere = True
                for ch in range(nOccChArr[index_s]):
                    n1 = n1Arr[index_s][ch]
                    n2 = n2Arr[index_s][ch]
                    if(m == n1 and n == n2):
                        canPlaceSiteHere = False
                        
                if(canPlaceSiteHere):
                    #print(f"site added at n1, n2 = {m}, {n}")
                    additionalX.append([b1[0]*float(m)+b2[0]*float(n)])
                    additionalY.append([b1[1]*float(m)+b2[1]*float(n)])
                    additionalVals.append(smolVal)
                    #plt.annotate('+',((b1[0]*float(m)+b2[0]*float(n)),(b1[1]*float(m)+b2[1]*float(n))),fontsize=8,zorder=10,ha='center',va='center',color=[1., 0., 0.])
                #else:
                    #print(f"no site at n1, n2 = {m}, {n}")

    #print(additionalX)
    plotCoordsArray = np.array(np.column_stack((np.append(plotCoordsX,additionalX), np.append(plotCoordsY,additionalY))))
    order = np.append(plotValuesAvg,additionalVals)
    points = plotCoordsArray
    vor = Voronoi(points=points,furthest_site=False)

    #plots the voronoi diagram on the second subplot
    if(vanity):
        lw = 0.
    else:
        lw = 0.5
    voronoi_plot_2d(vor, show_vertices =False, show_points =False, ax=ax2,line_width=lw,line_colors=[0.5, 0.5, 0.5])
        

    #colours the voronoi cells    
    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon), color=mapper.to_rgba(order[r]))
            #plt.fill(*zip(*polygon))
    ax2.set_aspect('equal')
    plt.xticks([])  
    plt.yticks([])
    ymin = min(brightSpotYArr[index_s])-1/2
    ymax = max(brightSpotYArr[index_s])+1/2
    xmin = min(brightSpotXArr[index_s])-1/2
    xmax = max(brightSpotXArr[index_s])+1/2
    ax2.set_ylim(ymin,ymax)
    ax2.set_xlim(xmin,xmax)
    #print("....")
    #print(ymin)
    #print(ymax)
    #print(xmin)
    #print(xmax)
    #print("....")


    filenametxt=""
    if(not vanity):
        plt.figtext(0.5, -0.035, defectstr+", "+entropytxt, wrap=True, horizontalalignment='center', fontsize=12,transform=ax2.transAxes)
        plt.figtext(0.5, -0.07, intenstr + ", " + simgastr+", "+kstr, wrap=True, horizontalalignment='center', fontsize=12,transform=ax2.transAxes)
    else:
        filenametxt = "vanity"
        plt.tight_layout(pad=0)

    if(filenametxt == ""):
        plt.figtext(0.5, -0.11, fileprefix, wrap=True, horizontalalignment='center', fontsize=12,fontstyle='italic',transform=ax2.transAxes)
        savestr = "Figures/Diffraction_multi/" + execTime + "_" + slugify(fileprefix) + "_" + scatcondstr +".png"
    else:
        plt.figtext(0.5, -0.11, filenametxt, wrap=True, horizontalalignment='center', fontsize=12,fontstyle='italic',transform=ax2.transAxes)
        savestr = "Figures/Diffraction_multi/" + execTime + "_" + slugify(filenametxt) + "_" + scatcondstr + ".png"
    print(savestr)
    if(vanity):
        plt.savefig(fname=savestr,dpi=1000)
    else:
        plt.savefig(fname=savestr,dpi=300)
    plt.show()