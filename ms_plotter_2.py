import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import datetime
import unicodedata
import re
import matplotlib.patheffects as pe
import math

#packages to import
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
import matplotlib.cm as cm

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
    dfs = np.array([])
    for index in range(0,nlines-2):
        val1 = d.loc[index,'n1']
        val2 = d.loc[index+1,'n1']
        difference = val2 - val1
        #print("val2 - val1 = "+ str(val2 - val1))
        if(difference < 0):
            print("\nSeperation found! Splitting dataframe between " + str(topOfSeperation) + ", " + str(index))
            dslice = np.array(d.iloc[topOfSeperation:index+1,:])
            topOfSeperation = index+1
            print(np.shape(dslice))

            if(index == 0):
                dfs = np.array(dslice)
                dfs = np.expand_dims(dfs,axis=2)
                print("dfs shape = " + str(np.shape(dfs)))
            else:
                dfs = np.append(dfs,np.expand_dims(dslice,axis=2),axis=2)
    dslice = np.array(d.iloc[topOfSeperation:,:])
    #print(dslice)
    dfs = np.append(dfs,dslice,axis=2)
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
    n = len(values)
    for v in np.nditer(values):
        mean += v
        meanSq += v**2
    mean /= n
    meanSq /= n
    stdv = np.sqrt(meanSq - mean**2)
    return mean, stdv

latticeFile = open('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script', 'r')
count = 0 

B1 = [0,0]
B2 = [0,0]
a1 = [0,0]
a2 = [0,0]
Nsuper = 1337
while count < 9:
    # Get next line from file
    line = latticeFile.readline()
    if(line.startswith("Nsuper = ")):
        line = line[len("Nsuper = "):]
        split = line.split()
        Nsuper = float(split[0])
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
        Nensemble = float(split[0])
        #print("Ensemble size = " + split[0])
        count += 1
    if line.startswith("Positional entropy = "):
        line = line[len("Positional entropy = "):]
        split = line.split()
        entropyIn = float(split[0])
        #print("Positional entropy = " + split[0])
        count += 1
    if line.startswith("Defect density in cm^-2 = "):
        line = line[len("Defect density in cm^-2 = "):]
        split = line.split()
        defectDensity = float(split[0])
        #print("Defect density in cm^-2 = " + split[0])
        count += 1

Babs = np.sqrt(B1[0]**2+B1[1]**2)
#print("B1 = ") 
#print(B1)
#print("B2 = ") 
#print(B2)
#print("Babs = ") 
#print(Babs)
latticeFile.close()
b1 = B1 / Babs
b2 = B2 / Babs

scatFile = open('scatCond.in', 'r')

numScatConds = sum(1 for _ in scatFile)
numScatConds -= 1

Es = []
thetas = []
phis = []

scatFile.close()
scatFile = open('scatCond.in', 'r')

line = scatFile.readline()
for index_n in range(0,numScatConds):
    line = scatFile.readline()
    print("Line = " + line)
    vals = line.split(",")
    Es.append( float(vals[0]) )
    thetas.append( float(vals[1]) )
    phis.append( float(vals[2]) )

print("[][][][][][][][]")
print("Number of scattering conditions = " + str(numScatConds))
print("[][][][][][][][]\n\n")

for index_n in range(0,int(Nensemble)):
    importname =  'diffrac' + str(10001+index_n) + '.out'
    print("importing file : " + importname)
    dfs = import_multiscat(importname)
    #print("===")
    #print(dfs)
    if(index_n == 0):
        dfss = dfs
    else:
        dfss = np.append(dfss,dfs)

for index_s in range(0,numScatConds):
    dfs = dfss[index_s][:]
    print("\n---\nscattering condition = " + str(index_s))
    print("\n#####")
    print(dfs)
    print("\n#####")
    df = dfs[0]
    #print("\n£££££ df = ")
    #print(df)
    nOccCh = len(df.index)
    plotValuesAvg = np.zeros((nOccCh))
    plotCoordsX = np.zeros((nOccCh))
    plotCoordsY = np.zeros((nOccCh))
    entropiesOut = np.zeros((int(Nensemble),1))

    for index_n in range(0,int(Nensemble)):
        print(index_n)
        intensities = np.zeros((nOccCh))
        plotCoords = np.zeros((nOccCh))
        plotCoords = np.zeros((nOccCh))
        df = dfs[index_n]
        nOccCh = len(df.index)
        pCXS = np.array([]) #These variables stand for something but icr what it is lmao
        pCYS = np.array([])
        n1min = 0
        n1max = 0
        n2min = 0
        n2max = 0
        valmin = 1
        valmax = 0
        smolVal = 1e-100
        vanityVal = 0

        print("nOccCh = " + str(nOccCh))
        for k in range(0,nOccCh):
            #print("k = " + str(k))
            row = df.iloc[k]
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
            I = float(getattr(row,'I'))
            if(I == 0):
                print(f"Zero found, setting to {smolVal}")
                I = smolVal
            elif(I > vanityVal):
                pCXS = np.append(pCXS,b1[0] * n1 + b2[0] * n2)
                pCYS = np.append(pCYS,b1[1] * n1 + b2[1] * n2)
            if(I < valmin):
                valmin = I
            if(I > valmax):
                valmax = I
            #print(f"k = {k}, n1 = {n1}, n2 = {n2}, I = {I}")
            intensities[k] = I
            plotCoordsX[k] = b1[0] * n1 + b2[0] * n2
            plotCoordsY[k] = b1[1] * n1 + b2[1] * n2
        #print(f"Valmin = {valmin}, valmax = {valmax}")
        #print("===")
        #print(f"n1min = {n1min}, n1max = {n1max}, n2min = {n2min}, n2max = {n2max}")
        #print("===")

        #print("Number of occupied channels = " + str(nOccCh))
        entropiesOut[index_n] = calculate_entropy(intensities)
        print("Diffraction pattern entropy = " + str(entropiesOut[index_n]))
        plotValuesAvg += intensities / Nensemble

    #sets the colour scale
    useLog = True
    if(useLog):
        mapper = cm.ScalarMappable(cmap='magma', norm=mpl.colors.LogNorm(valmin,valmax))
    else:
        mapper = cm.ScalarMappable(cmap='magma', norm=mpl.colors.Normalize(valmin,valmax))

    #create the figure with set figure size
    fig = plt.figure(figsize=(10,8))

    #creates two subplots
    ax = plt.subplot2grid((16,20), (0,17), colspan=1, rowspan=16)
    ax2 = plt.subplot2grid((16,20), (0,0), colspan=16, rowspan=16)

    pathefts1 = [pe.Stroke(linewidth=1, foreground='w'), pe.Normal()]
    pathefts2 = [pe.Stroke(linewidth=2, foreground='w'), pe.Normal()]
    b1col = [0.9, 0.1, 0]
    b2col = [0, 0.7, 0.2]
    a1col = [1, 0.5, 0.6]
    a2col = [0.5, 0.8, 0.6]
    hecol = [0, 0.3, 0.8]

    plt.arrow(0,0,b1[0]*Nsuper,b1[1]*Nsuper,width=0.05,color=b1col,zorder=7,path_effects=pathefts2,length_includes_head=True)
    plt.arrow(0,0,b2[0]*Nsuper,b2[1]*Nsuper,width=0.05,color=b2col,zorder=7,path_effects=pathefts2,length_includes_head=True)
    plt.annotate("b1", (b1[0]*Nsuper,b1[1]*Nsuper+0.1),color=b1col,fontsize=8,weight='bold',path_effects=pathefts1,zorder=11)
    plt.annotate("b2", (b2[0]*Nsuper,b2[1]*Nsuper),color=b2col,fontsize=8,weight='bold',path_effects=pathefts1,zorder=11)

    plt.arrow(0,0,a1[0]/np.sqrt(a1[0]**2+a1[1]**2),a1[1]/np.sqrt(a1[0]**2+a1[1]**2),width=0.05,color=a1col,zorder=6,length_includes_head=True,alpha=0.5)
    plt.arrow(0,0,a2[0]/np.sqrt(a1[0]**2+a1[1]**2),a2[1]/np.sqrt(a1[0]**2+a1[1]**2),width=0.05,color=a2col,zorder=6,length_includes_head=True,alpha=0.5)
    plt.annotate("a1", (a1[0]/np.sqrt(a1[0]**2+a1[1]**2),a1[1]/np.sqrt(a1[0]**2+a1[1]**2)),color=a1col,fontsize=8,weight='bold',zorder = 5)
    plt.annotate("a2", (a2[0]/np.sqrt(a1[0]**2+a1[1]**2),a2[1]/np.sqrt(a1[0]**2+a1[1]**2)),color=a2col,fontsize=8,weight='bold',zorder = 5)


    for k in range(0,nOccCh):
        row = df.iloc[k]
        n1 = int(getattr(row,'n1'))
        n2 = int(getattr(row,'n2'))
        if(n1%int(Nsuper) == 0 and n2%int(Nsuper)==0):
            n1n2 = str(int(n1/Nsuper)) + ',' + str(int(n2/Nsuper))
            if(plotValuesAvg[k] < (valmax-valmin)*0.9):
                col = 'w'
            else:
                col = 'k'
            plt.annotate(n1n2,((b1[0]*n1+b2[0]*n2),(b1[1]*n1+b2[1]*n2)),fontsize=8,zorder=10,ha='center',va='center',c=col)
    print("index = " + str(index_s))    
    #print(thetas)
    E = Es[index_s]
    theta = thetas[index_s]
    #print(theta)
    phi = phis[index_s]

    
    titelstr = "$E$ = " + str(E) + " meV, $\\theta$ = " + str(theta) + "$\\degree$, $\\phi$ =" + str(phi) + "$\\degree$"
    scatcondstr = str(E) + "_" + str(theta) + "_" + str(phi)
    print(scatcondstr)
    scatFile.close()

    heliumRot = np.matrix([[np.cos(np.deg2rad(phi)),np.sin(np.deg2rad(phi))],
                        [-np.sin(np.deg2rad(phi)),np.cos(np.deg2rad(phi))]])
    heliumDir = -heliumRot * np.reshape(np.array(a1),(2,1))/np.sqrt(np.dot(a1,a1))
    e = 1.602176634E-19
    hbar = 6.62607015e-34/(2*np.pi)
    m = 4 * 1.66053906660e-27
    #heliumk = np.sin(np.deg2rad(theta))*heliumDir * np.sqrt(2*m*e*E/1000)/hbar
    heliumk = heliumDir * np.sqrt(2*m*e*E/1000)/hbar
    heliumk_n = heliumk/(Babs*1e10)
    print("heliumk_n =")
    print(heliumk_n)
    if(not(math.isclose(theta,0.) & math.isclose(phi,0.))):
        plt.arrow(0,0,Nsuper*heliumk_n[0,0],Nsuper*heliumk_n[1,0],width=0.03,color='b',zorder=7,head_width=0.1)
    ax2.add_patch(plt.Circle((0, 0), np.sqrt(heliumk_n[0]**2 + heliumk_n[1]**2)*Nsuper, color='b', fill=False,zorder=7,linestyle=(0, (5, 10))))
    ax2.set_title(titelstr)

    #creates a colourbar on the first subplot
    if(useLog):
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap='magma', norm=mpl.colors.LogNorm(valmin,valmax), orientation='vertical')
    else:
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap='magma', norm=mpl.colors.Normalize(valmin,valmax), orientation='vertical')
    cb1.set_label('P($n_1$,$n_2$)')


    additionalX = []
    additionalY = []
    additionalVals = []

    vanity = True
    if(vanity):
        paddingCells = 30
        for n in range(-paddingCells,paddingCells):
            for m in range(-paddingCells,paddingCells):
                canPlaceSiteHere = True
                for k in range(0,nOccCh):
                    row = df.iloc[k]
                    n1 = int(getattr(row,'n1'))
                    n2 = int(getattr(row,'n2'))
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
    order = np.append(intensities,additionalVals)
    points = plotCoordsArray
    vor = Voronoi(points=points,furthest_site=False)

    #plots the voronoi diagram on the second subplot
    voronoi_plot_2d(vor, show_vertices =False, show_points =False, ax=ax2,line_width=0.5,line_colors=[0.5, 0.5, 0.5])
        

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
    ax2.set_ylim(min(pCYS)-1/2,max(pCYS)+1/2)
    ax2.set_xlim(min(pCXS)-1/2,max(pCXS)+1/2)

    eomean, eosdtv = find_mean_stdv(entropiesOut)
    captiontxt="$n_{defect}$ = " + "{:.4e}".format(defectDensity) + " cm$^{-2}$, $H_{defect}$ = " + "{:.4f}".format(entropyIn)
    if(Nensemble == 1):
        entropytxt = "$H_{diffraction}$ = " + "{:.6f}".format(eomean)
    else:
        entropytxt = "$H_{diffraction}$ = " + "{:.6f}".format(eomean) + "$\pm$" +  "{:.6f}".format(eosdtv)

    plt.figtext(0.5, -0.035, captiontxt, wrap=True, horizontalalignment='center', fontsize=12,transform=ax2.transAxes)
    plt.figtext(0.5, -0.07, entropytxt, wrap=True, horizontalalignment='center', fontsize=12,transform=ax2.transAxes)
    filenametxt=""
    filenametxt="Circular defect - flat"
    plt.figtext(0.5, -0.11, filenametxt, wrap=True, horizontalalignment='center', fontsize=12,fontstyle='italic',transform=ax2.transAxes)

    if(filenametxt == ""):
        savestr = "Figures/Diffraction_multi/" + datetime.datetime.now().strftime('_%Y-%m-%d_%H-%M') + "_" + scatcondstr +".png"
    else:
        savestr = "Figures/Diffraction_multi/" + datetime.datetime.now().strftime('_%Y-%m-%d_%H-%M') +slugify(filenametxt) + "_" + scatcondstr + ".png"
    print(savestr)
    plt.savefig(fname=savestr,dpi=300)
    plt.show()