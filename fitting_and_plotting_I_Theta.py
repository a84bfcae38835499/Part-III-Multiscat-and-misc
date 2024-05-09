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

Thetas = np.array([0,1/25,2/25,3/25,4/25,5/25,6/25])

Thetas = [x/25 for x in range(0,10)]

Thetas += [(5*x)/25 for x in range(2,6)]

Thetas = np.array(Thetas)
print(Thetas)

#Thetas = np.array([0,1/25,2/25,3/25,4/25,5/25,6/25,7/25,8/25,9/25])
N = Thetas.size
Thetas_continuum = np.linspace(0.,1.,500)
""" vv This data is for real mos2 """

I10 = np.array([1.,0.8343576,0.6972710,0.5820808,0.4816842,0.3830502,0.3218771
])
I10unc = np.array([1e-10,0.0007715,0.0007532,0.0007406,0.0010848,0.0014003,0.0037750
])

I00 = np.array([1.,0.8343576,0.6732104,0.5442669,0.4312601,0.3237506,0.2276699
])
I00unc = np.array([1e-10,0.0007715,0.0010537,0.0011096,0.0016252,0.0020978,0.0056555
])

I_10 = np.array([1.,0.8343576,0.8768258,0.8629333,0.8210708,0.8091048,0.8938186
])
I_10unc = np.array([1e-10,0.0007715,0.0075451,0.0079450,0.0116372,0.0150213,0.0404958
])

I_20 = np.array([1.,0.8854526,0.7679356,0.6512175,0.5378604,0.4310878,0.3931136
])
I_20unc = np.array([1e-10,0.0016455,0.0015001,0.0015796,0.0023136,0.0029864,0.0080510
])

I_30 = np.array([1.,0.8402811,0.6926791,0.5560454,0.4404330,0.3453218,0.2964200
])
I_30unc = np.array([1e-10,0.0020416,0.0018611,0.0019598,0.0028706,0.0037053,0.0099891
])

""" v This data is for the gaussianz """
gaI00 = np.array([1.,0.9235979,0.8848201,0.8652465,0.7989249,0.7926799,0.7571613,0.7605055,0.6751379,0.6843397
])
gaI00unc = np.array([1e-10,0.0176584,0.0371947,0.0399589,0.0459412,0.0309872,0.0766873,0.0529496,0.0724498,0.0557927
])


""" v This data is the total data for the gaussians at 0 degrees """
gaI_T = np.array([
    1.,0.9625589,0.9264127,0.8953515,0.8593757,0.8377734,0.7931277,0.7697033,0.7398964,0.7224723,0.6815588,0.6187126,0.5258513,0.5605555
])
gaI_Tunc = np.array([
    1e-10,0.0000473,0.0000770,0.0001449,0.0001461,0.0002042,0.0002520,0.0002935,0.0002922,0.0003071,0.0004033,0.0004915,0.0006300,0.0060960
])

gvI_T = np.array([1.,
])
gvI_Tunc = np.array([1e-10,
])

n1n2OfInterest = []
n1n2Colours = [0., 0., 0.]
"""
n1n2OfInterest = [[1,0],
                  [0,0],
                  [-1,0],
                  [-2,0],
                  [-3,0]]
n1n2Colours = [[1, 0.82, 0.149],
                  [0., 0., 0.],
                  [0.067, 0.769, 0.451],
                  [0.149, 0.792, 0.7],
                  [0.296, 0.369, 1],]
n1n2IArr = [I10,I00,I_10,I_20,I_30]
n1n2IUncArr = [I10unc,I00unc,I_10unc,I_20unc,I_30unc]
"""
"""
n1n2OfInterest = [[-1,0]]
n1n2Colours = [[0.067, 0.769, 0.451]]
n1n2IArr = [I_10]
n1n2IUncArr = [I_10unc]

n1n2OfInterest = [[0,0]]
n1n2Colours = [[0., 0., 0.]]
n1n2IArr = [gaI00]
n1n2IUncArr = [gaI00unc]

"""


Ninterest = len(n1n2OfInterest)

def i_tb_L(Theta,S):
    I_I0 = np.power((1 - Theta),2*S/cellArea)
    return I_I0
n_tb_L = 1

def i_tb_M(Theta,S):
    eta = np.sqrt(S/cellArea)
    return (1-Theta*(eta**2)+(Theta**2)*(10-16*eta+7*eta**2)-(Theta**3)*(20-28*eta+10*eta**2)+(Theta**4)*(9-12*eta+4*eta**2))**2
n_tb_M = 1

def i_tb_S(Theta,S):
    eta = np.sqrt(S/cellArea)
    return (1-Theta*(eta**2)+(Theta**2)*(1-3*eta)*(1-eta)-2*(Theta**3)*((1-eta)**2))**2
n_tb_S = 1

def i_tb_T(Theta,S):
    eta = np.sqrt(S/cellArea)
    return 1-2*Theta*eta**2
n_tb_T = 1
#===============================================
Ninterest = 1
for i in range(Ninterest):
    #n1n2 = n1n2OfInterest[i]
    #I = n1n2IArr[i]
    #IUnc = n1n2IUncArr[i]
    I = gaI_T
    IUnc = gaI_Tunc

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #plt.title("Intensity against $\Theta$ for "+str(n1n2))
    plt.title("Intensity against $\Theta$ for Gaussian Adatoms")
    #if(Ninterest > 0):
    #    ax1.set_ylabel("I"+str(n1n2)+"/$I_0$")
    #else:
    ax1.set_ylabel("$I/I_0$")
    ax1.set_xlabel("$\Theta$")
    ax1.set_xlim(0,1)
    ax1.set_ylim(0,1)
    fillclr = [1.,0.3,0.5]
    errscale = 20
    #if(errscale == 1):
    #    ax1.errorbar(x=Thetas,y=I,yerr=IUnc*errscale,label="I"+str(n1n2)+"\n(error bars $\\times$"+str(errscale)+")",
    #                color=n1n2Colours[i],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)
    #else:
    #    ax1.errorbar(x=Thetas,y=I,yerr=IUnc,label="I"+str(n1n2),
    #                color=n1n2Colours[i],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)
    if(errscale == 1):
        ax1.errorbar(x=Thetas,y=I,yerr=IUnc,label="$I/I_0$",
                    color=fillclr,linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)
    else:
        ax1.errorbar(x=Thetas,y=I,yerr=IUnc*errscale,label="$I/I_0$"+"\n(error bars $\\times$"+str(errscale)+")",
                    color=fillclr,linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)


    print("* - - - - - -")

    p_t, c_t, infodict_t, null, null = curve_fit(
        i_tb_T,Thetas,I,
        sigma=IUnc,absolute_sigma=True,
        p0 = 4.,check_finite = True,nan_policy='raise',bounds=[0,cellArea],full_output=True
        )
    S_tiny = p_t[0]
    SE = np.sqrt(np.diag(c_t))
    S_tiny_Unc = SE[0]
    print("S_tiny = "+ "{:.7f}".format(S_tiny) + "±" + "{:.7f}".format(S_tiny_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_t,N,n_tb_S)))

    print("* - - - - - -")

    p_s, c_s, infodict_s, null, null = curve_fit(
        i_tb_S,Thetas,I,
        sigma=IUnc,absolute_sigma=True,
        p0 = cellArea*1.1,check_finite = True,nan_policy='raise',bounds=[cellArea,cellArea*(9/4)],full_output=True
        )
    S_small = p_s[0]
    SE = np.sqrt(np.diag(c_s))
    S_small_Unc = SE[0]
    print("S_small = "+ "{:.7f}".format(S_small) + "±" + "{:.7f}".format(S_small_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_s,N,n_tb_S)))
    print("* - - - - - -")

    p_m, c_m, infodict_m, null, null = curve_fit(
        i_tb_M,Thetas,I,
        sigma=IUnc,absolute_sigma=True,
        p0 = 3*cellArea,check_finite = True,nan_policy='raise',bounds=[(9/4)*cellArea,4*cellArea],full_output=True
        )
    S_medium = p_m[0]
    SE = np.sqrt(np.diag(c_m))
    S_medium_Unc=SE[0]
    print("S_medium = "+ "{:.7f}".format(S_medium) + "±" + "{:.7f}".format(S_medium_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_m,N,n_tb_M)))
    print("* = = = = = =")

    p_l, c_l, infodict_l, null, null = curve_fit(
        i_tb_L,Thetas,I,
        sigma=IUnc,absolute_sigma=True,
        p0 = cellArea*5,check_finite = True,nan_policy='raise',bounds=[4*cellArea,np.inf],full_output=True
        )
    S_large = p_l[0]
    SE = np.sqrt(np.diag(c_l))
    S_large_Unc = SE[0]
    print("S_large = "+ "{:.7f}".format(S_large) + "±" + "{:.7f}".format(S_large_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_l,N,n_tb_L)))
    print("* - - - - - -")


    #===============================================
    clr = np.array([1, 0, .2])
    fD_L = i_tb_L(Thetas_continuum,S_large)
    fD_L_L = i_tb_L(Thetas_continuum,S_large-S_large_Unc)
    fD_L_U = i_tb_L(Thetas_continuum,S_large+S_large_Unc)
    ax1.plot(Thetas_continuum,fD_L,label="Large $S$\n"+
             "S=" + "{:.5f}".format(S_large) + "$\pm$" + "{:.5f}".format(S_large_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_L_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_L_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_L_L, fD_L_U, alpha=0.3, edgecolor=clr, facecolor=clr)
    
    #===============================================
    clr = np.array([.3, .8, .0])
    fD_tb_medium = i_tb_M(Thetas_continuum,S_medium)
    fD_tb_medium_L = i_tb_M(Thetas_continuum,S_medium-S_medium_Unc)
    fD_tb_medium_U = i_tb_M(Thetas_continuum,S_medium+S_medium_Unc)
    ax1.plot(Thetas_continuum,fD_tb_medium,label="Medium $S$\n"+
             "S=" + "{:.5f}".format(S_medium) + "$\pm$" + "{:.5f}".format(S_medium_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_tb_medium_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_tb_medium_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_tb_medium_L, fD_tb_medium_U, alpha=0.3, edgecolor=clr, facecolor=clr)

    #===============================================
    """
    clr = np.array([1.0, .5, .0])
    fD_tb_medium_r = i_tb_M_R(Thetas_continuum,eta_tb_medium_r,xi_tb_medium_r,rho_def_tb_medium_r,delta_def_tb_medium_r)
    fD_tb_medium_r_L = i_tb_M_R(Thetas_continuum,eta_tb_medium_r-eta_tb_medium_r_Unc,xi_tb_medium_r-xi_tb_medium_r_Unc,rho_def_tb_medium_r-rho_def_tb_medium_r_Unc,delta_def_tb_medium_r-delta_def_tb_medium_r_Unc)
    fD_tb_medium_r_U = i_tb_M_R(Thetas_continuum,eta_tb_medium_r+eta_tb_medium_r_Unc,xi_tb_medium_r+xi_tb_medium_r_Unc,rho_def_tb_medium_r+rho_def_tb_medium_r_Unc,delta_def_tb_medium_r+delta_def_tb_medium_r_Unc)
    ax1.plot(Thetas_continuum,fD_tb_medium_r,label="Medium $\Sigma$, reflecting\n"+
             "S=" + "{:.5f}".format(S_tb_medium_r) + "$\pm$" + "{:.5f}".format(S_tb_medium_r_Unc)+"\n"
             "$\\rho$ = " +"{:.5f}".format(rho_def_tb_medium_r) + "$\pm$" + "{:.5f}".format(rho_def_tb_medium_r_Unc)+ "\n" +
             "$\delta$ = " +"{:.5f}".format(delta_def_tb_medium_r) + "$\pm$" + "{:.5f}".format(delta_def_tb_medium_r_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_tb_medium_r_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_tb_medium_r_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_tb_medium_r_L, fD_tb_medium_r_U, alpha=0.3, edgecolor=clr, facecolor=clr)
    """
    #===============================================
    """
    clr = np.array([1.0, .5, .0])
    fD_tb_small_r = i_tb_M_R(Thetas_continuum,eta_tb_small_r,xi_tb_small_r,rho_def_tb_small_r,delta_def_tb_small_r)
    fD_tb_small_r_L = i_tb_M_R(Thetas_continuum,eta_tb_small_r-eta_tb_small_r_Unc,xi_tb_small_r-xi_tb_small_r_Unc,rho_def_tb_small_r-rho_def_tb_small_r_Unc,delta_def_tb_small_r-delta_def_tb_small_r_Unc)
    fD_tb_small_r_U = i_tb_M_R(Thetas_continuum,eta_tb_small_r+eta_tb_small_r_Unc,xi_tb_small_r+xi_tb_small_r_Unc,rho_def_tb_small_r+rho_def_tb_small_r_Unc,delta_def_tb_small_r+delta_def_tb_small_r_Unc)
    ax1.plot(Thetas_continuum,fD_tb_small_r,label="Small $\Sigma$, reflecting\n"+
            "S=" + "{:.5f}".format(S_tb_small_r) + "$\pm$" + "{:.5f}".format(S_tb_small_r_Unc)+"\n"
            "$\\rho$ = " +"{:.5f}".format(rho_def_tb_small_r) + "$\pm$" + "{:.5f}".format(rho_def_tb_small_r_Unc)+ "\n" +
            "$\delta$ = " +"{:.5f}".format(delta_def_tb_small_r) + "$\pm$" + "{:.5f}".format(delta_def_tb_small_r_Unc),
            color=clr)
    ax1.plot(Thetas_continuum,fD_tb_small_r_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_tb_small_r_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_tb_small_r_L, fD_tb_small_r_U, alpha=0.3, edgecolor=clr, facecolor=clr)
    """
    #===============================================
    clr = np.array([0, 0.4, 0.4])
    fD_tb_small = i_tb_M(Thetas_continuum,S_small)
    fD_tb_small_L = i_tb_M(Thetas_continuum,S_small-S_small_Unc)
    fD_tb_small_U = i_tb_M(Thetas_continuum,S_small+S_small_Unc)
    ax1.plot(Thetas_continuum,fD_tb_small,label="Small S\n"+
             "S=" + "{:.5f}".format(S_small) + "$\pm$" + "{:.5f}".format(S_small_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_tb_small_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_tb_small_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_tb_small_L, fD_tb_small_U, alpha=0.3, edgecolor=clr, facecolor=clr)

    #===============================================
    clr = np.array([0, 0.7, 1.0])
    fD_tb_tiny = i_tb_M(Thetas_continuum,S_tiny)
    fD_tb_tiny_L = i_tb_M(Thetas_continuum,S_tiny-S_tiny_Unc)
    fD_tb_tiny_U = i_tb_M(Thetas_continuum,S_tiny+S_tiny_Unc)
    ax1.plot(Thetas_continuum,fD_tb_small,label="Tiny S\n"+
             "S=" + "{:.5f}".format(S_tiny) + "$\pm$" + "{:.5f}".format(S_tiny_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_tb_tiny_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_tb_tiny_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_tb_tiny_L, fD_tb_tiny_U, alpha=0.3, edgecolor=clr, facecolor=clr)

    #===============================================
    plt.legend(loc='upper right')
    #plt.savefig(fname="Figures/Fitted Plots/"+slugify( execTime +"_fitted_" + str(n1n2)),dpi=300)
    plt.savefig(fname="Figures/Fitted Plots/"+slugify( execTime +"_fitted_ga"),dpi=300)
    plt.show()
    print("* Saved figure")
