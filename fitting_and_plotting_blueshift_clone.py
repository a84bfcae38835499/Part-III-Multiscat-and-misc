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

#Thetas = np.array([x/25 for x in range(1,9)]+[x/25 for x in range(10,30,5)])
Thetas = np.array([x/25 for x in range(0,10)])
N = Thetas.size
print(Thetas)
Thetas_continuum = np.linspace(0.,1.,500)


Adatom_H_diff = np.array([0.482082, 0.520447,0.548192,0.572525,0.593844,0.614375,0.639014,0.654021,0.673540,0.688574
])
Adatom_H_diff_Unc = np.array([0, 0.004248,0.003308,0.004460,0.006082,0.006033,0.006615,0.006312,0.005530,0.005966
])

Adatom_K = np.array([1.8196899,1.8319695,1.8632271,1.8917125,1.8843851,1.9361424,1.9477458,1.9530034,1.9903159,2.0129283
])
Adatom_K_Unc = np.array([0,0.0262528,0.0497557,0.0474716,0.0893587,0.0878799,0.0696341,0.0722962,0.0910749,0.0931299
])


Vacancy_H_diff = np.array([0.482082,0.518472,0.546852,0.571712,0.594904,0.621270,0.636488,0.653600,0.674068,0.687810
])
Vacancy_H_diff_Unc = np.array([0,0.004146,0.002800,0.007577,0.004823,0.004684,0.005989,0.005131,0.004122,0.011551

])
                                                                                        #   v extremely anomalous shit goin
Vacancy_K = np.array([1.8196899,1.8337429,1.8547909,1.8771613,1.9104914,1.9690134,1.9689019,2.0050719,1.9985184,2.0403300
])
Vacancy_K_Unc = np.array([0,0.0188508,0.0375486,0.0332651,0.0922589,0.0565959,0.0807228,0.0896489,0.0740899,0.0852711
])
#how can I contrive to write this line of code:
#   phi lo mean a Kunc


fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("Entropy and |K| against $\Theta$")
ax1.set_ylabel("H$_{diff}$")
ax2 = ax1.twinx()
ax2.set_ylabel("|$K$|/$Å^{-1}$")
ax1.set_xlabel("$\Theta$")
ax1.set_xlim(min(Thetas),max(Thetas))
ax1.set_ylim(min(Adatom_H_diff),max(Adatom_H_diff))
ax1.errorbar(x=Thetas,y=Adatom_H_diff,yerr=Adatom_H_diff_Unc,label="Adatom H$_{diff}$",
                 color=[0.,0.,1.],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)
ax2.errorbar(x=Thetas,y=Adatom_K,yerr=Adatom_K_Unc,label="Adatom |$K$|",
                 color=[1,.5,0.],linestyle=(0,(4,10)),marker='+',markersize=10,capsize=2.)

ax1.errorbar(x=Thetas,y=Vacancy_H_diff,yerr=Vacancy_H_diff_Unc,label="Vacancy H$_{diff}$",
                 color=[1,0,0.5],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)
ax2.errorbar(x=Thetas,y=Vacancy_K,yerr=Vacancy_K_Unc,label="Vacancy |$K$|",
                 color=[0.,0.5,0.5],linestyle=(0,(4,10)),marker='+',markersize=10,capsize=2.)
ax1.legend(loc=1)
ax2.legend(loc=0)

plt.legend()
plt.savefig(fname="Figures/Blueshift/"+slugify( execTime ) + "_blueshift",dpi=300)
plt.show()
print("* Saved figure")
