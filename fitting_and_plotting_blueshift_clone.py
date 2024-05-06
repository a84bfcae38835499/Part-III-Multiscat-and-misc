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

Thetas = np.array([x/25 for x in range(0,9)]+[x/25 for x in range(10,30,5)])
N = Thetas.size
print(Thetas)
Thetas_continuum = np.linspace(0.,1.,500)

K = np.array([1.8196899,1.8319695,])
K_Unc = np.array([0,0.0262528,])

#how can I contrive to write this line of code:
#   phi lo mean a Kunc

H_diff = np.array([1.,0.8343576,0.6972710,0.5820808,0.4816842,0.3830502,0.3218771
])
H_diff_Unc = np.array([1e-10,0.0007715,0.0007532,0.0007406,0.0010848,0.0014003,0.0037750
])


plt.legend()
plt.savefig(fname="Figures/Blueshift/"+slugify( execTime +"_yarrrr me hearties"),dpi=300)
plt.show()
print("* Saved figure")
