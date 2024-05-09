from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import datetime
import unicodedata
import re
invMaxTheta = 1
fileprefix = '_5x5_05D'

plt.rcParams["font.family"] = "PT Sans"

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

def calculated_chi_squared(infodict,N,n):
    return (infodict['fvec']**2).sum()/(N-n)

execTime = datetime.datetime.now().strftime('_%Y-%m-%d_%H-%M')
latticeFile = open(fileprefix + '.info_for_vivian_python_nice_plotting_hexagon_script', 'r')
count = 0

B1 = np.array([1e-10,0.])
B2 = np.array([1e-10,0.])
a1 = np.array([1e-10,0.])
a2 = np.array([1e-10,0.])
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
        print("Theta = " + split[0])
        count += 1
    if line.startswith("Ensemble size = "):
        line = line[len("Ensemble size = "):]
        split = line.split()
        Nensemble = int(split[0])
        print("Ensemble size = " + split[0])
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
        print("Defect density in cm^-2 = " + split[0])
        count += 1
ftg = 15/180 * np.pi
invCosThetas = np.array([1/np.cos(x*ftg) for x in range(0,6)])
print(invCosThetas)

n1n2OfInterest = [
                  [0,0],
                  [-1,0],
                  [-2,0],
                  [-3,0],
                  ]
n1n2Colours = [
                  [0., 0., 0.],
                  [0.067, 0.769, 0.451],
                  [0.149, 0.792, 0.7],
                  [0.296, 0.369, 1],
                  ]

print("--------------")
fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("$S$ against $\theta$ for "+str(n1n2))
ax1.set_ylabel("I"+str(n1n2)+"/$I_0$")
ax1.set_xlabel("$\Theta$")
ax1.set_xlim(0,1/3)
ax1.set_ylim(0,1)
ax1.errorbar(x=Thetas,y=I,yerr=IUnc*20,label="I"+str(n1n2)+"\n(error bars $\\times$20)",
                 color=n1n2Colours[i],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)
plt.legend()
plt.savefig(fname="Figures/Blueshift/"+slugify( execTime ) + "_blueshift_vacancy",dpi=300)
print("* Saved figure")
plt.show()