from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chisquare

invMaxTheta = 1
fileprefix = '_5x5_05D'


def calculated_chi_squared(infodict):
    return (infodict['fvec']**2).sum()/(len(infodict[2]['fvec'])-len(infodict[0]))

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

Thetas = np.array([0.,1/25,2/25,3/25,4/25,5/25,6/25])
N = Thetas.size
Thetas_continuum = np.linspace(0.,1.,500)
I10 = np.array([1.,0.8343576,0.6972710,0.5820808,0.4816842,0.3830502,0.3218771
])
I10unc = np.array([0.,0.0007715,0.0007532,0.0007406,0.0010848,0.0014003,0.0037750
])

I00 = np.array([1.,0.8343576,0.6732104,0.5442669,0.4312601,0.3237506,0.2276699
])
I00unc = np.array([0.,0.0007715,0.0010537,0.0011096,0.0016252,0.0020978,0.0056555
])

I_10 = np.array([1.,0.8343576,0.8768258,0.8629333,0.8210708,0.8091048,0.8938186
])
I_10unc = np.array([0.,0.0007715,0.0075451,0.0079450,0.0116372,0.0150213,0.0404958
])

I_20 = np.array([1.,0.8854526,0.7679356,0.6512175,0.5378604,0.4310878,0.3931136
])
I_20unc = np.array([0.,0.0016455,0.0015001,0.0015796,0.0023136,0.0029864,0.0080510
])

I_30 = np.array([1.,0.8402811,0.6926791,0.5560454,0.4404330,0.3453218,0.2964200
])
I_30unc = np.array([0.,0.0020416,0.0018611,0.0019598,0.0028706,0.0037053,0.0099891
])



def i_L(Theta,S): #This is our y-data, I/I_0. We Fit S to the real data
    I_I0 = np.power((1-invMaxTheta * Theta),(2*S)/(invMaxTheta*cellArea))
    return I_I0
n_L = 1#Number of parameters OBSOLETE u can delete these

def i_S_D_NOO(Theta,S): #This is our y-data, I/I_0. We Fit S to the real data
    def F(Theta,R,D):
        return (1-Theta)*((1-Theta)**6+(R-D/6)/(R)*Theta*(1-Theta)**3+(R-D/3)/(R)*(Theta**2)*(1-Theta)+(R-D/2)/(R)*(Theta**3))
    I_I0 = F(Theta,cellArea,S-cellArea)**2
    return I_I0
n_S_D_NOO = 1

def i_S_R_2_NOO(Theta,S,phase,rho): #This is our y-data, I/I_0. We Fit S to the real data
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

def i_tb_S_R(Theta,eta,xi,rho,phi):
    return np.abs(i_tb_S_D(Theta,eta)+rho*np.exp(1.j * phi)*i_tb_S_S(Theta,xi))**2
n_tb_S_R = 4

def i_tb_M_R(Theta,eta,xi,rho,phi):
    return np.abs(i_tb_M_D(Theta,eta)+rho*np.exp(1.j * phi)*i_tb_M_S(Theta,xi))**2
n_tb_M_R = 4
#===============================================


fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.title("Intensity against defect density")
ax1.set_ylabel("$I[1,0]/I_0$")
ax1.set_xlabel("$\Theta$")
ax1.set_xlim(0,1/3)
ax1.set_ylim(0,1)
ax1.errorbar(x=Thetas,y=I00,yerr=I00unc,label="I[0,0]",color=[0., 0.5, 1.],linestyle=(0,(8,10)),marker='x',markersize=10,capsize=2.)

fillclr = [0.3,1.0,0.0]
print("* - - - - - -")

p_L, c_L, infodict_L, null, null = curve_fit(
    i_L,Thetas,I00,
    p0 = 10,check_finite = True,nan_policy='raise',bounds=[0,np.inf],full_output=True
    )
crossSection_L = p_L[0]
SE = np.sqrt(np.diag(c_L))
crossSection_L_Unc = SE[0]
print(p_L)
print(c_L)
print("crossSection_L = "+ "{:.7f}".format(crossSection_L) + "±" + "{:.7f}".format(crossSection_L_Unc))
print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_L,N,n_L)))
print("* - - - - - -")

p_S_D_NOO, c_S_D_NOO, infodict_S_D_NOO, null, null = curve_fit(
    i_S_D_NOO,Thetas,I00,
    p0 = 10,check_finite = True,nan_policy='raise',bounds=[0,np.inf],full_output=True
    )
crossSection_S_D_NOO = p_S_D_NOO[0]
SE = np.sqrt(np.diag(c_S_D_NOO))
crossSection_S_D_NOO_Unc = SE[0]
print(p_S_D_NOO)
print(c_S_D_NOO)
print("crossSection_S_D_NOO = "+ "{:.7f}".format(crossSection_S_D_NOO) + "±" + "{:.7f}".format(crossSection_S_D_NOO_Unc))
print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_S_D_NOO)))
print("* - - - - - -")

p_S_R_2_NOO, c_S_R_2_NOO, infodict_S_R_2_NOO, null, null = curve_fit(
    i_S_R_2_NOO,Thetas,I00,
    p0 = [10.,0.,1.],check_finite = True,nan_policy='raise',full_output=True,
    bounds=[[0,0,0],[1000,2*np.pi,1.]]
    )
crossSection_S_R_2_NOO = p_S_R_2_NOO[0]
phase_S_R_2_NOO = p_S_R_2_NOO[1]
rho_S_R_2_NOO = p_S_R_2_NOO[2]
SE = np.sqrt(np.diag(c_S_R_2_NOO))
crossSection_S_R_2_NOO_Unc = SE[0]
phase_S_R_2_NOO_Unc = SE[1]
rho_S_R_2_NOO_Unc = SE[2]

print(p_S_R_2_NOO)
print(c_S_R_2_NOO)
print("crossSection_S_R_2_NOO = "+ "{:.7f}".format(crossSection_S_R_2_NOO) + "±" + "{:.7f}".format(crossSection_S_R_2_NOO_Unc))
print("phase_S_R_2_NOO = "+ "{:.7f}".format(phase_S_R_2_NOO) + "±" + "{:.7f}".format(phase_S_R_2_NOO_Unc))
print("rho_S_R_2_NOO = "+ "{:.7f}".format(rho_S_R_2_NOO) + "±" + "{:.7f}".format(rho_S_R_2_NOO_Unc))
print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_S_R_2_NOO)))
print("* = = = = = =")

p_tb_medium, c_tb_medium, infodict_tb_medium, null, null = curve_fit(
    i_tb_M_D,Thetas,I00,
    p0 = 2.,check_finite = True,nan_policy='raise',bounds=[0,np.inf],full_output=True
    )
crossSection_tb_medium = cellArea * (p_S_R_2_NOO[0])**2
SE = np.sqrt(np.diag(c_S_R_2_NOO))
crossSection_tb_medium_Unc = cellArea * SE[0]**2

print(p_tb_medium)
print(c_tb_medium)
print("crossSection_tb_medium = "+ "{:.7f}".format(crossSection_tb_medium) + "±" + "{:.7f}".format(crossSection_tb_medium_Unc))
print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_tb_medium)))
print("* = = = = = =")

p_tb_medium_r, c_tb_medium_r, infodict_tb_medium_r, null, null = curve_fit(
    i_tb_M_R,Thetas,I00,
    p0 = [2.,2.,0.5,0.],check_finite = True,nan_policy='raise',full_output=True,
    bounds=[[1.5,1.5,0.,0.],[2.0,2.0,1.0,np.pi*2]]
    )
print(p_tb_medium_r)
print(c_tb_medium_r)

eta = p_tb_medium_r[0]
crossSection_tb_medium_r = cellArea * eta**2 #eta
xi = p_tb_medium_r[1]
crossSection_def_tb_medium_r = cellArea * xi**2 #xi
rho_def_tb_medium_r = (p_tb_medium_r[2]) #rho
phi_def_tb_medium_r = (p_tb_medium_r[3]) #phi
SE = np.sqrt(np.diag(c_tb_medium_r))
eta_Unc = SE[0]
crossSection_tb_medium_r_Unc = cellArea * eta_Unc**2
xi_Unc = SE[1]
crossSection_def_tb_medium_r_Unc = cellArea * xi_Unc**2
rho_def_tb_medium_r_Unc = SE[2]
phi_def_tb_medium_r_Unc = SE[3]

print("crossSection_tb_medium_r = "+ "{:.7f}".format(crossSection_tb_medium_r) + "±" + "{:.7f}".format(crossSection_tb_medium_r_Unc))
print("crossSection_def_tb_medium_r = "+ "{:.7f}".format(crossSection_def_tb_medium_r) + "±" + "{:.7f}".format(crossSection_def_tb_medium_r_Unc))
print("rho_def_tb_medium_r = "+ "{:.7f}".format(rho_def_tb_medium_r) + "±" + "{:.7f}".format(rho_def_tb_medium_r_Unc))
print("phi_def_tb_medium_r = "+ "{:.7f}".format(phi_def_tb_medium_r) + "±" + "{:.7f}".format(phi_def_tb_medium_r_Unc))
print("chi-squared = " + "{:.7f}".format(calculated_chi_squared(infodict_tb_medium_r)))
print("* = = = = = =")


#===============================================

fD_L = i_L(Thetas_continuum,crossSection_L)
fL_L = i_L(Thetas_continuum,crossSection_L-crossSection_L_Unc)
fL_U = i_L(Thetas_continuum,crossSection_L+crossSection_L_Unc)

fD_S_D_NOO = i_S_D_NOO(Thetas_continuum,crossSection_S_D_NOO)
fD_S_D_NOO_L = i_S_D_NOO(Thetas_continuum,crossSection_S_D_NOO-crossSection_S_D_NOO_Unc)
fD_S_D_NOO_U = i_S_D_NOO(Thetas_continuum,crossSection_S_D_NOO+crossSection_S_D_NOO_Unc)

fD_S_R_2_NOO = i_S_R_2_NOO(Thetas_continuum,crossSection_S_R_2_NOO,phase_S_R_2_NOO,rho_S_R_2_NOO)
fD_S_R_2_NOO_L = i_S_R_2_NOO(Thetas_continuum,crossSection_S_R_2_NOO-crossSection_S_R_2_NOO_Unc,phase_S_R_2_NOO-phase_S_R_2_NOO_Unc,rho_S_R_2_NOO-rho_S_R_2_NOO_Unc)
fD_S_R_2_NOO_U = i_S_R_2_NOO(Thetas_continuum,crossSection_S_R_2_NOO+crossSection_S_R_2_NOO_Unc,phase_S_R_2_NOO+phase_S_R_2_NOO_Unc,rho_S_R_2_NOO+rho_S_R_2_NOO_Unc)

fD_tb_medium_r = i_tb_S_R(Thetas_continuum,np.sqrt(crossSection_def_tb_medium_r/cellArea),np.sqrt(crossSection_def_tb_medium_r/cellArea),rho_def_tb_medium_r,phi_def_tb_medium_r)

ax1.plot(Thetas_continuum,fD_L,label="Large $\Sigma$",color=[1., .0, .0])
ax1.plot(Thetas_continuum,fD_S_D_NOO,label="Small $\Sigma$",color=[.0, .8, .2])
ax1.plot(Thetas_continuum,fD_S_R_2_NOO,label="Small $\Sigma$, reflection",color=[.0, .2, 0.8])
ax1.plot(Thetas_continuum,fD_tb_medium_r,label="Textbook medium $\Sigma$, reflecting",color=[.0, .0, 0.0])
#ax1.fill_between(Thetas_continuum, fLCSD_L, fLCSD_U, alpha=0.3, edgecolor=fillclr, facecolor=fillclr)
plt.legend()
plt.show()
