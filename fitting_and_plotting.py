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
N = Thetas.size
Thetas_continuum = np.linspace(0.,1.,500)

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
n1n2OfInterest = [[-1,0]]
n1n2Colours = [[0.067, 0.769, 0.451]]
n1n2IArr = [I_10]
n1n2IUncArr = [I_10unc]
"""
Ninterest = len(n1n2OfInterest)

def i_L(Theta,S):
    I_I0 = np.power((1 - Theta),2*S/cellArea)
    return I_I0
n_L = 1

def i_S_D_NOO(Theta,S):
    def F(Theta,R,D):
        return (1-Theta)*((1-Theta)**6+(R-D/6)/(R)*Theta*(1-Theta)**3+(R-D/3)/(R)*(Theta**2)*(1-Theta)+(R-D/2)/(R)*(Theta**3))
    I_I0 = F(Theta,cellArea,S-cellArea)**2
    return I_I0
n_S_D_NOO = 1

def i_S_R_2_NOO(Theta,S,phase,rho):
    def F(Theta,R,D):
        return (1-Theta)*((1-Theta)**6+(R-D/6)/(R)*Theta*(1-Theta)**3+(R-D/3)/(R)*(Theta**2)*(1-Theta)+(R-D/2)/(R)*(Theta**3))
    I_I0 = F(Theta,cellArea,S-cellArea)**2 + (rho*Theta)**2 + 2*rho*Theta*F(Theta,cellArea,S-cellArea)*np.cos(phase)
    return I_I0
n_S_R_2_NOO = 3

def i_tb_S_D(Theta,eta):
    return (1-Theta*(eta**2)+(Theta**2)*(1-3*eta)*(1-eta)-2*(Theta**3)*((1-eta)**2))**2
n_tb_S_D = 1

def i_tb_M_D(Theta,eta):
    return (1-Theta*(eta**2)+(Theta**2)*(7*(eta**2)-16*eta+10)+(Theta**3)*(10*(eta**2)-28*eta+20)+(Theta**4)* (4*(eta**2)-12*eta+9))**2
n_tb_M_D = 1

def i_tb_S_S(Theta,xi):
    return Theta*((xi-1)**2)-(Theta**2)*(xi-1)*(5-3*xi)+2*(Theta**3)*((1-xi)**2)
n_tb_S_S = 1

def i_tb_M_S(Theta,xi):
    return Theta*((xi-2)**2)-(Theta**2)*((xi-2)**2)-2*(Theta**3)*(xi-2)*(3*xi-4)+(Theta**4)*((2*xi-3)**2)
n_tb_M_S = 1

def i_tb_S_R(Theta,eta,xi,rho,delta):
    return np.abs(i_tb_S_D(Theta,eta)+rho*np.exp(1.j * delta)*i_tb_S_S(Theta,xi))**2
n_tb_S_R = 4

def i_tb_M_R(Theta,eta,xi,rho,delta):
    return np.abs(i_tb_M_D(Theta,eta)+rho*np.exp(1.j * delta)*i_tb_M_S(Theta,xi))**2
n_tb_M_R = 4
#===============================================

for i in range(Ninterest):
    n1n2 = n1n2OfInterest[i]
    I = n1n2IArr[i]
    IUnc = n1n2IUncArr[i]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title("Intensity against $\Theta$ for "+str(n1n2))
    ax1.set_ylabel("I"+str(n1n2)+"/$I_0$")
    ax1.set_xlabel("$\Theta$")
    ax1.set_xlim(0,1/3)
    ax1.set_ylim(0,1)
    ax1.errorbar(x=Thetas,y=I,yerr=IUnc*20,label="I"+str(n1n2)+"\n(error bars $\\times$20)",
                 color=n1n2Colours[i],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)

    fillclr = [0.3,1.0,0.0]

    print("* - - - - - -")

    p_L, c_L, infodict_L, null, null = curve_fit(
        i_L,Thetas,I,
        #sigma=IUnc,absolute_sigma=True,
        p0 = 1.,check_finite = True,nan_policy='raise',bounds=[0,np.inf],full_output=True
        )
    S_L = p_L[0]
    SE = np.sqrt(np.diag(c_L))
    S_L_Unc = SE[0]
    #print(p_L)
    #print(c_L)
    print("S_L = "+ "{:.7f}".format(S_L) + "±" + "{:.7f}".format(S_L_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_L,N,n_L)))
    print("* - - - - - -")

    p_tb_small, c_tb_small, infodict_tb_small, null, null = curve_fit(
        i_tb_S_D,Thetas,I,
        #sigma=IUnc,absolute_sigma=True,
        p0 = 2.,check_finite = True,nan_policy='raise',bounds=[0,np.inf],full_output=True
        )
    eta_tb_small = p_tb_small[0]
    S_tb_small = cellArea * (eta_tb_small)**2
    SE = np.sqrt(np.diag(c_tb_small))
    eta_tb_small_Unc = SE[0]
    S_tb_small_Unc = cellArea * eta_tb_small_Unc**2

    #print(p_tb_small)
    #print(c_tb_small)
    print("S_tb_small = "+ "{:.7f}".format(S_tb_small) + "±" + "{:.7f}".format(S_tb_small_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_tb_small,N,n_tb_S_D)))

    print("* - - - - - -")

    p_tb_medium, c_tb_medium, infodict_tb_medium, null, null = curve_fit(
        i_tb_M_D,Thetas,I,
        #sigma=IUnc,absolute_sigma=True,
        p0 = 2.,check_finite = True,nan_policy='raise',bounds=[0,np.inf],full_output=True
        )
    eta_tb_medium = p_tb_medium[0]
    S_tb_medium = cellArea * (eta_tb_medium)**2
    SE = np.sqrt(np.diag(c_tb_medium))
    eta_tb_medium_Unc=SE[0]
    S_tb_medium_Unc = cellArea * eta_tb_medium_Unc**2

    #print(p_tb_medium)
    #print(c_tb_medium)
    print("S_tb_medium = "+ "{:.7f}".format(S_tb_medium) + "±" + "{:.7f}".format(S_tb_medium_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_tb_medium,N,n_tb_M_D)))
    print("* = = = = = =")

    p_tb_small_r, c_tb_small_r, infodict_tb_small_r, null, null = curve_fit(
        i_tb_S_R,Thetas,I,
        #sigma=IUnc,absolute_sigma=True,
        p0 = [1.5,1.5,0,np.pi],check_finite = True,nan_policy='raise',full_output=True,
        bounds=[[1,1,0,0],[1.5,1.5,1,2*np.pi]]
        )
    #print(p_tb_small_r)
    #print(c_tb_small_r)

    eta_tb_small_r = p_tb_small_r[0]
    S_tb_small_r = cellArea * eta_tb_small_r**2 #eta
    xi_tb_small_r = p_tb_small_r[1]
    S_def_tb_small_r = cellArea * xi_tb_small_r**2 #xi
    rho_def_tb_small_r = (p_tb_small_r[2]) #rho
    delta_def_tb_small_r = (p_tb_small_r[3]) #delta
    SE = np.sqrt(np.diag(c_tb_small_r))
    eta_tb_small_r_Unc = SE[0]
    S_tb_small_r_Unc = cellArea * eta_tb_small_r_Unc**2
    xi_tb_small_r_Unc = SE[1]
    S_def_tb_small_r_Unc = cellArea * xi_tb_small_r_Unc**2
    rho_def_tb_small_r_Unc = SE[2]
    delta_def_tb_small_r_Unc = SE[3]

    print("S_tb_small_r = "+ "{:.7f}".format(S_tb_small_r) + "±" + "{:.7f}".format(S_tb_small_r_Unc))
    print("S_def_tb_small_r = "+ "{:.7f}".format(S_def_tb_small_r) + "±" + "{:.7f}".format(S_def_tb_small_r_Unc))
    print("rho_def_tb_small_r = "+ "{:.7f}".format(rho_def_tb_small_r) + "±" + "{:.7f}".format(rho_def_tb_small_r_Unc))
    print("delta_def_tb_small_r = "+ "{:.7f}".format(delta_def_tb_small_r) + "±" + "{:.7f}".format(delta_def_tb_small_r_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_tb_small_r,N,n_tb_S_R)))
    print("* = = = = = =")

    p_tb_medium_r, c_tb_medium_r, infodict_tb_medium_r, null, null = curve_fit(
        i_tb_M_R,Thetas,I,
        #sigma=IUnc,absolute_sigma=True,
        p0 = [1.8,1.8,0.,np.pi],check_finite = True,nan_policy='raise',full_output=True,
        bounds=[[1.5,1.5,0,0],[2,2,1,2*np.pi]]
        )
    #print(p_tb_medium_r)
    #print(c_tb_medium_r)

    eta_tb_medium_r = p_tb_medium_r[0]
    S_tb_medium_r = cellArea * eta_tb_medium_r**2 #eta
    xi_tb_medium_r = p_tb_medium_r[1]
    S_def_tb_medium_r = cellArea * xi_tb_medium_r**2 #xi
    rho_def_tb_medium_r = (p_tb_medium_r[2]) #rho
    delta_def_tb_medium_r = (p_tb_medium_r[3]) #delta
    SE = np.sqrt(np.diag(c_tb_medium_r))
    eta_tb_medium_r_Unc = SE[0]
    S_tb_medium_r_Unc = cellArea * eta_tb_medium_r_Unc**2
    xi_tb_medium_r_Unc = SE[1]
    S_def_tb_medium_r_Unc = cellArea * xi_tb_medium_r_Unc**2
    rho_def_tb_medium_r_Unc = SE[2]
    delta_def_tb_medium_r_Unc = SE[3]

    print("S_tb_medium_r = "+ "{:.7f}".format(S_tb_medium_r) + "±" + "{:.7f}".format(S_tb_medium_r_Unc))
    print("S_def_tb_medium_r = "+ "{:.7f}".format(S_def_tb_medium_r) + "±" + "{:.7f}".format(S_def_tb_medium_r_Unc))
    print("rho_def_tb_medium_r = "+ "{:.7f}".format(rho_def_tb_medium_r) + "±" + "{:.7f}".format(rho_def_tb_medium_r_Unc))
    print("delta_def_tb_medium_r = "+ "{:.7f}".format(delta_def_tb_medium_r) + "±" + "{:.7f}".format(delta_def_tb_medium_r_Unc))
    print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_tb_medium_r,N,n_tb_M_R)))
    print("* = = = = = =")


    #===============================================
    clr = np.array([1, 0, .2])
    fD_L = i_L(Thetas_continuum,S_L)
    fD_L_L = i_L(Thetas_continuum,S_L-S_L_Unc)
    fD_L_U = i_L(Thetas_continuum,S_L+S_L_Unc)
    ax1.plot(Thetas_continuum,fD_L,label="Large $\Sigma$\n"+
             "S=" + "{:.5f}".format(S_L) + "$\pm$" + "{:.5f}".format(S_L_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_L_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_L_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_L_L, fD_L_U, alpha=0.3, edgecolor=clr, facecolor=clr)
    
    #===============================================
    clr = np.array([.3, .8, .0])
    fD_tb_medium = i_tb_M_D(Thetas_continuum,eta_tb_medium)
    fD_tb_medium_L = i_tb_M_D(Thetas_continuum,eta_tb_medium-eta_tb_medium_Unc)
    fD_tb_medium_U = i_tb_M_D(Thetas_continuum,eta_tb_medium+eta_tb_medium_Unc)
    ax1.plot(Thetas_continuum,fD_tb_medium,label="Medium $\Sigma$\n"+
             "S=" + "{:.5f}".format(S_tb_medium) + "$\pm$" + "{:.5f}".format(S_tb_medium_Unc),
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
    fD_tb_small = i_tb_M_D(Thetas_continuum,eta_tb_small)
    fD_tb_small_L = i_tb_M_D(Thetas_continuum,eta_tb_small-eta_tb_small_Unc)
    fD_tb_small_U = i_tb_M_D(Thetas_continuum,eta_tb_small+eta_tb_small_Unc)
    ax1.plot(Thetas_continuum,fD_tb_small,label="Small $\Sigma$\n"+
             "S=" + "{:.5f}".format(S_tb_small) + "$\pm$" + "{:.5f}".format(S_tb_small_Unc),
             color=clr)
    ax1.plot(Thetas_continuum,fD_tb_small_L,color=clr/2)
    ax1.plot(Thetas_continuum,fD_tb_small_U,color=clr/2)
    ax1.fill_between(Thetas_continuum, fD_tb_small_L, fD_tb_small_U, alpha=0.3, edgecolor=clr, facecolor=clr)

    #===============================================
    plt.legend()
    plt.savefig(fname=slugify("Figures/Fitted Plots/" + execTime +"_fitting_" + str(n1n2)),dpi=300)
    plt.show()
    print("* Saved figure")
