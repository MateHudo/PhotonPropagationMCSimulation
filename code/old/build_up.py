#%% import
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Users/Admin/Desktop/Python_koda')
from scipy.optimize import curve_fit
import sympy as sym
import scipy as sc
#import pandas as pd
import numpy.random as rndm
#from numpy.random import random as RND
from numpy.random import uniform as UNFM
#from numpy.random import randint 
#from matplotlib.animation import FuncAnimation
import time

import scienceplots
plt.style.use(['science','notebook','grid'])

path_uproot = 'C:\\Users\\Admin\\anaconda3\\Lib\\site-packages'
sys.path.append(path_uproot)
print(sys.path)

sys.version

import uproot

#%%  XCOM, ESTAR data
   # MAC = Mass Absorption Coefficient
E_gamma, MAC_coh, MAC_compton, MAC_fe, MAC_tp_jedro, MAC_tp_el, MAC_tot, MAC_wo_coh = [],[], [],[],[],[], [],[]
with open('C:/Users/Admin/Desktop/MAGISTERIJ/4. semester/ZF - zdravstvnea_fizika/simulacije/1-Bulid-up faktor/attenuation_coefficient_Lead.txt','r') as file:
    for vrstica in file:
        vrstica = vrstica.split()
        #print(vrstica[0])
        E_gamma.append(float(vrstica[0]))
        MAC_coh.append(float(vrstica[1]))
        MAC_compton.append(float(vrstica[2]))
        MAC_fe.append(float(vrstica[3]))
        MAC_tp_jedro.append(float(vrstica[4]))
        MAC_tp_el.append(float(vrstica[5]))
        MAC_tot.append(float(vrstica[6]))
        MAC_wo_coh.append(float(vrstica[7]))
E_gamma,MAC_coh,MAC_compton,MAC_fe = np.array(E_gamma), np.array(MAC_coh),np.array(MAC_compton),np.array(MAC_fe)
MAC_tp_jedro,MAC_tp_el,MAC_tot,MAC_wo_coh = np.array(MAC_tp_jedro),np.array(MAC_tp_el),np.array(MAC_tot),np.array(MAC_wo_coh)

  # interaction names
interaction_types = ['coherent','compton','fe','tp_nucleus','tp_el','total','wo_coherent']

   # AC = Absorption Coefficient
Pb_density = 11.35
AC_coh,AC_compton,AC_fe = Pb_density*MAC_coh, Pb_density*MAC_compton, Pb_density*MAC_fe
AC_tp_jedro,AC_tp_el,AC_tot,AC_wo_coh = Pb_density*MAC_tp_jedro, Pb_density*MAC_tp_el, Pb_density*MAC_tot, Pb_density*MAC_wo_coh


#%% CHATgpt INTERPOLATION (energijska intrpolacija abs koeficientov)
def AC_interpolation(E0, E, AC):
    """
    Linearly interpolates the absorption coefficient at energy E0.
    
    Parameters:
        E0 (float): The energy at which to determine the absorption coefficient.
        E (array-like): Array of energies (must be in ascending order).
        AC (array-like): Array of absorption coefficients corresponding to the energies in E.

    Returns:
        float: Interpolated absorption coefficient at energy E0.
    """
    # Find the index where E0 should be inserted
    index = np.searchsorted(E, E0, side='right') - 1
        
    # If E0 is outside the range of E, return the absorption coefficient at the closest energy
    if index < 0:
        return AC[0]
    elif index >= len(E) - 1:
        return AC[-1]
    # če je genau ena od tableiranih energij (predvsem če začetno energijo tko izbereš...
    elif E[index] == E0:
        return AC[index]
    
    # Perform linear interpolation
    x0, x1 = E[index], E[index + 1]
    y0, y1 = AC[index], AC[index + 1]
    return y0 + (y1 - y0) * (E0 - x0) / (x1 - x0)

# Example usage
E = np.array([1, 2, 3, 4, 5])
AC = np.array([0.5, 0.7, 0.9, 1.2, 1.5])
E0 = 2
# Interpolate the absorption coefficient at energy E0
interpolated_AC = AC_interpolation(E0, E, AC)
print("Interpolated absorption coefficient at energy", E0, "is:", interpolated_AC)

# primer iz xcom podatkov 
E0 = 2.2
interpolated_AC = AC_interpolation(E0, E_gamma, AC_tot)
print("Interpolated absorption coefficient at energy", E0, "is:", interpolated_AC)
E_gamma[np.searchsorted(E_gamma,E0,side='right')-1]
AC_tot[np.searchsorted(E_gamma,E0,side='right')-1]

#AC_tp_jedro
AC_interpolation(10,E_gamma,AC_tot)

#%% Interacition length
def InteractionLength(E0,E,AC):
    ac = AC_interpolation(E0,E,AC)
    s_int = - np.log(UNFM()) / ac
    return s_int
    
E0 = 1  # MeV
print(AC_interpolation(E0,E_gamma,AC_compton))
print(AC_interpolation(E0,E_gamma,AC_fe))
print(AC_interpolation(E0,E_gamma,AC_fe)/AC_interpolation(E0,E_gamma,AC_compton))


N_paths = 10**5
s_fe = [InteractionLength(E0,E_gamma,AC_fe) for _ in range(N_paths)]
s_com = [InteractionLength(E0,E_gamma,AC_compton) for _ in range(N_paths)]

plt.hist(s_com,bins=50,density='True',color='r',alpha=0.4,label='compton')
plt.hist(s_fe,bins=50,density='True',color='b',alpha=0.9,label='FE')
plt.xlim(0,30)
plt.legend()

 #%% (x,y,z) <-> (phi,theta)
def Cartesian_to_spherical(r):
    """
    Converts Cartesian coordinates to spherical coordinates.
    
    Parameters:
        r (array-like): Array containing the Cartesian coordinates (x, y, z).
        
    Returns:
        phi (float): Azimuthal angle in radians.
        theta (float): Polar angle in radians.
    """
    x, y, z = r
    # Calculate azimuthal angle (theta)
    theta = np.arctan2(y, x)
    # Calculate polar angle (phi)
    r_xy = np.sqrt(x**2 + y**2)
    phi = np.arctan2(r_xy, z)

    return theta, phi

omega = (1,1,0)
print(Cartesian_to_spherical(omega))

def Spherical_to_cartesian(direction):
    min_value = 10**(-15)
    theta, phi = direction
    #Sphi = np.sin(phi) if np.sin(phi) 
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    if x == 0:
        x = min_value
    if y == 0:
        y = min_value
    if z==0:
        z = min_value    
    return np.array([x,y,z])


#omega = (1,1,10)
#Cartesian_to_spherical(omega)
dir = (np.pi/2,0)
Spherical_to_cartesian(dir)   

#%% Compton 1
def ComptonEnergyGenerator1(E):
    n = 0
    a = E / 0.511   # MeV
    u1,u2,u3 = UNFM(), UNFM(), UNFM()
    e_min = 1/(1+2*a)
    M1, M2 = 1/2 * (1-e_min**2), -np.log(e_min)
    if u1 > M1/(M1+M2):
        e_proposed = e_min * np.exp(M2*u2)
    else:
        e_proposed = np.sqrt(e_min**2 + 2*M1*u2)
    f = (1-e_proposed)/a/e_proposed
    h= 1 - e_proposed/(1+e_proposed**2) * f * (2-f)
    if u3 < h:
        n += 1
        return e_proposed
    else:
        return False
        #return 'e_poposed failed'

ComptonEnergyGenerator1(E=1)

E = 1
N = 10**5
CEG = [ComptonEnergyGenerator1(E) for i in range(N)]


#%% Compton 2 - boljš (while True)
def ComptonEnergyGenerator2(E):
    #n1,n2 = 0,0
    a = E / 0.511   # MeV
    e_min = 1/(1+2*a)
    M1, M2 = 1/2 * (1-e_min**2), -np.log(e_min)
    while True:
        u1,u2,u3 = UNFM(), UNFM(), UNFM()
        if u1 > M1/(M1+M2):
            e_proposed = e_min * np.exp(M2*u2)
        else:
            e_proposed = np.sqrt(e_min**2 + 2*M1*u2)
        f = (1-e_proposed)/a/e_proposed
        h= 1 - e_proposed/(1+e_proposed**2) * f * (2-f)
        if u3 < h:
            return e_proposed


def ComptonEnergyGeneratorEfficiency(E):
    #n1,n2 = 0,0
    a = E / 0.511   # MeV
    e_min = 1/(1+2*a)
    M1, M2 = 1/2 * (1-e_min**2), -np.log(e_min)
    u1,u2,u3 = UNFM(), UNFM(), UNFM()
    if u1 > M1/(M1+M2):
        e_proposed = e_min * np.exp(M2*u2)
    else:
        e_proposed = np.sqrt(e_min**2 + 2*M1*u2)
    f = (1-e_proposed)/a/e_proposed
    h= 1 - e_proposed/(1+e_proposed**2) * f * (2-f)
    if u3 < h:
        return True
    else:
        return False

E = 1
N = 10**4
CEG = [ComptonEnergyGeneratorEfficiency(E) for _ in range(N)]
print(len(CEG))
accepted_samples = sum(ComptonEnergyGeneratorEfficiency(E) for _ in range(N))
accepted_samples

def EfficiencyAnalysis(E_list,N):
    eff = []
    for E in E_list:
        accepted_samples = sum(ComptonEnergyGeneratorEfficiency(E) for _ in range(N))
        eff.append(accepted_samples/N*100)
    #return eff
    plt.plot(E_list,eff,'r--',marker='o',markersize=5)
    plt.xlabel('energija [MeV]')
    plt.ylabel('raven sprejemnaja generiranih energij [%]')

E_list = [.1,.2,.3,.4,.5,0.8,1,2,3,4,5,6,8,10]
EfficiencyAnalysis(E_list,N=10**3)  # rabš velko cifro za statistiko!!    


    
#%%   Compton plotting distributions
# initial incident photon energy E0  
E0 = .01
alpha = E0/0.511
e_min = 1/(1+2*alpha)  # min energy of photon after CS
   
   # normalising PDF . sympy in scipy nomralizacija
# symbols
e_sym, a_sym,th_sym, cos_th_sym = sym.symbols('e_sym a_sym th_sym cos_th_sym')  # e = eps = E'/e; a_sym = alpha = E/m_el

 # pdf functions that i want to normalise and then also plot vs RN generator histogram  
   # i) dw/d(eps) ... spremelnjivka je eps!!!
pdf_eps = e_sym + 1/e_sym + (1-e_sym)/(a_sym*e_sym) * (-2 + (1-e_sym)/(a_sym*e_sym)) 
pdf_eps_func = sym.lambdify(e_sym,pdf_eps.subs({a_sym:alpha}),'numpy')
eps_plot = np.linspace(e_min,1,100)
pdf_eps_plot = pdf_eps_func(eps_plot)
#plt.plot(eps_plot,pdf_plot,'r--')
 
   # ii) dw/d(cos(th))  .. spremenljikva je cos(th)!!!
e_cos_sym = 1 /(1 + a_sym*(1-cos_th_sym))
sin_th_sym  = sym.sqrt(1 - cos_th_sym**2)
pdf_omega = e_cos_sym**2 * (e_cos_sym + 1/e_cos_sym - sin_th_sym**2)
pdf_cos_th_func = sym.lambdify(cos_th_sym,pdf_omega.subs({a_sym:alpha}),'numpy')
cos_plot = np.linspace(-1,1,100)
kot_plot = np.arccos(cos_plot)
pdf_omega_plot =pdf_cos_th_func(cos_plot)
#plt.plot(cos_plot,pdf_cos_th_func(cos_plot),'r--')
plt.plot(kot_plot,pdf_omega_plot,'r--')

   # iii) dw/d(th)  ... spremenljivka je theta !!!
e_th_sym = 1 /(1 + a_sym*(1-sym.cos(th_sym)))
pdf_th = e_th_sym**2 * (e_th_sym + 1/e_th_sym - sym.sin(th_sym)**2)  * sym.sin(th_sym)
pdf_th_func = sym.lambdify(th_sym,pdf_th.subs({a_sym:alpha}),'numpy')
th_plot = np.linspace(0,np.pi,100)
pdf_theta_plot = pdf_th_func(th_plot)
#plt.plot(th_plot,pdf_theta_plot,'r--')


  # normalisation function - normalizes pdfs so integral is 1
def SymbolicNormalization(pdf,int_sym,int_boundaries,dict_of_consts):
    pdf_int = sym.integrate(pdf.subs(dict_of_consts),(int_sym,int_boundaries[0],int_boundaries[1]))
    C_norm = 1 / pdf_int.evalf()
    return C_norm

def NumericNormalization(pdf,int_sym,boundaries,dict_of_consts):
    pdf_eval = pdf.subs(dict_of_consts)
    pdf_num = sym.lambdify(int_sym,pdf_eval,'numpy')
    num_int, err = sc.integrate.quad(pdf_num,boundaries[0],boundaries[1])
    return 1 / num_int 

# 
  # test
#e_min = 1/(1+2*alpha)  # min energy of photon after CS
C1_sym = SymbolicNormalization(pdf_eps,e_sym,(e_min,1),{a_sym:alpha})
C1_num = NumericNormalization(pdf_eps,e_sym,(e_min,1),{a_sym:alpha})
C2_num = NumericNormalization(pdf_omega,cos_th_sym,(-1,1),{a_sym:alpha}) 
C3_num = NumericNormalization(pdf_th,th_sym,(0,sym.pi),{a_sym:alpha}) 
#print(C1_sym)
#print(C1_num)
print(C2_num)
print(C3_num)


def KleinNishinaPDF(E0=1,N=100):
    # vrne tri porazdelitve: P2 = do/d(omega), P3 = do/d(theta)=P2 * sin(th)*2*pi
    # in P1 = do/d(eps)
    alpha = E0/0.511
    e_min = 1/(1+2*alpha)
    C1_num = NumericNormalization(pdf_eps,e_sym,(e_min,1),{a_sym:alpha})
    C2_num = NumericNormalization(pdf_omega,cos_th_sym,(-1,1),{a_sym:alpha}) 
    C3_num = NumericNormalization(pdf_th,th_sym,(0,sym.pi),{a_sym:alpha}) 
    
    # do/d(eps)
    eps_vals = np.linspace(e_min,1,N)
    P1 = C1_num * pdf_eps_func(eps_vals)
    
    # do/d(omega)
    delta = 0
    th_vals = np.linspace(0+delta,np.pi-delta,N)    
    cos_vals = np.cos(th_vals)
    P2 = C2_num * pdf_cos_th_func(cos_vals)
    
    # do/d(theta)
    P3 = C3_num * pdf_th_func(th_vals)
    
    return [(eps_vals,P1), (th_vals,P2), (th_vals,P3)]

#KleinNishinaPDF(E0)[0][0]

  # plot
#E0 = 1
#xx, yy = KleinNishinaPDF(E0)[1]
#xx = np.arccos(xx)
#plt.plot(xx,yy)
   

   #  generacija energij po CS
N_gama = 10**5
ComptonEnergies = [ComptonEnergyGenerator2(E0) for i in range(N_gama)]
#ComptonEnergies = [ComptonEnergies[i] for i in range(len(ComptonEnergies)) if ComptonEnergies[i] is not False]

  # kot sipanja fotona
def ComptonAngle(e,E0=1):
    a = E0 / 0.511
    cos_th = 1- (1-e)/(a*e)
    angle = np.arccos(cos_th)
    return [angle, angle / np.sin(angle)]
    #return [cos_th*np.sin(angle),cos_th]
    #return [cos_th, cos_th / 2 ] #np.sin(np.arccos(cos_th))]

ComptonAngles = np.array([ComptonAngle(energy,E0) for energy in ComptonEnergies]).T
#če daš E0 =1, je težava !! ??????
#ComptonAngles[0]
 
   # dod: porazdelitev po d(cos(th)) 
#def ComptonCosAngle(e,E0):
     

def FuncNorm(fun,xmin,xmax):
    integral, err = sc.integrate.quad(fun,xmin,xmax)
    return 1 / integral    
#FuncNorm(pdf_cos_th_func,-1,1)

   # teoretična poradelitev - KN
def KN_PDF(kot,E0=1):
    # vrne tri porazdelitve: P1 = do/d(omega), P2 = do/d(theta)=P1 * sin(th)*2*pi
    # in P3 = do/d(eps)
    a = E0 / 0.511
    e = 1/(1+a*(1-np.cos(kot)))
    P1 = (e + 1/e - np.sin(kot)**2)
    P2 = e**2 * P1
    P3 = P2 * np.sin(kot) 
    return P1
    #return [P1, P2, P3]


           
kot_plot = np.linspace(0,np.pi,100)
cos_kot_plot = np.linspace(-1,1,100)
arccos_kot = np.arccos(cos_kot_plot)
energija_plot = np.linspace(e_min,1,100)

#C = FuncNorm(pdf_cos_th_func,-1,1)
#C = NumericNormalization(e_sym,)
#PDF_kn = KN_PDF(kot_plot,E0)
#p1, p2, p3 = C1_num * PDF_kn[0], C2_num * PDF_kn[1], C3_num * PDF_kn[2]

#plt.plot(kot_plot,p1)
#plt.plot(energija_plot,PDF_kn)

#%%   ###  PLOTTING distributions
izbira = 2 # 0-energija, 1 - prostorski kot, 2-theta
   # 1 MI NE RATA DOBIT PORAZDELITVE!!! 

 # histogram
plt.figure(figsize=(10,6))
fig, ax = plt.subplots(figsize=(10,6))

     # porszdelitev po energijah
#plt.hist(ComptonEnergies,bins=100,density=True,label=f'energy generator histogram (N = {N_gama:.0e})')
     # porazdelitev po kotih
hist_plot = [ComptonEnergies,ComptonAngles[0],ComptonAngles[0]][izbira]
hist_label = [f'energy generator histogram (N = {N_gama:.0e})',f'angle generator histogram (N = {N_gama:.0e})',f'angle generator histogram (N = {N_gama:.0e})'][izbira]
plt.hist(hist_plot,bins=100,density=True,label=hist_label)
if izbira > 0:
    plt.xticks([0,np.pi/2,np.pi],[0,r'$\pi/2$',r'$\pi$'])     
#plt.xticks([])

  # plot theoretical pdf
x_axis = [energija_plot,arccos_kot,kot_plot][izbira]
#y_axis = [p1,p2,p3][izbira]
y_axis = [C1_num*pdf_eps_plot,C2_num*pdf_omega_plot,C3_num*pdf_theta_plot][izbira]
plot_label = [r'PDF \hspace dw/d$\epsilon$',r'PDF dw/$d\Omega$',r'$PDF dw/d\theta$'][izbira]
plt.plot(x_axis,y_axis,'r--',label='Probability Density Function')
#plt.xticks([0,np.pi/2,np.pi],[0,r'$\frac{\pi}{2}$',r'$\pi$'])
xlabel = [r'$\epsilon$',r'$\theta_{scatt}$',r'$\theta_{scatt}$'][izbira]
ylabel = [r'PDF dw/d$\epsilon$',r'$dw/d\Omega$',r'$dw/d\theta$'][izbira]
plt.xlabel(xlabel)
plt.ylabel(ylabel)
#plt.ylabel(r'PDF dw/d$\theta$')

plt.text(0.4,.3,s=f'E = {E0} MeV',fontsize=20,transform=fig.transFigure)
plt.legend(fontsize=14)
plt.show()

#plt.plot(eps_plot,C1_num*pdf_plot)




#%% # deluje??
# i) simbolično integriranje
# C1_sym = ... glej zgorej
integral_norm = sym.integrate(pdf_eps.subs({a_sym:alpha})*C1_sym,(e_sym,e_min,1)).evalf()
#integral_norm = sym.integrate(pdf_eps.subs({a_sym:alpha}),(e_sym,e_min.subs({a_sym:alpha}),1)).evalf()
print(integral_norm)
# ii) numerično integriranje
pdf_eps_func = sym.lambdify(e_sym,pdf_eps.subs({a_sym:alpha}),'numpy')
print(pdf_eps_func(.3))

#pdf_eps_func
#integral_num, error = sc.integrate.quad(pdf_eps_func,0.5,1)
integral_num, error = sc.integrate.quad(pdf_eps_func,e_min,1)
print(integral_num)
C1_num = 1  / integral_num
print(f'  C1_sym={C1_sym};\n  C1_num={C1_num}')

 # plot
eps_plot = np.linspace(e_min,1,30)
prob_plot = C1_sym*pdf_eps_func(eps_plot)
plt.plot(eps_plot,prob_plot,'--')



#%%   PATHTOEXIT - Pot do izhoda kvadra in ploskev, skozi katero izstopi
   # 1. prva ideja
def PathToExit1(d,a,b,r0,direction):
    x0,y0,z0 = r0
    #if x0<0 or

    theta, phi = direction    
    omega = Spherical_to_cartesian(direction)
    #omega = (np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi))
    #omega = [max(10**(-5),faktor) for faktor in omega]
    
    lam2, lam1 = 1/omega[0] * np.array([d-x0,-x0])
    lam4,lam3 = 1/omega[1] * np.array([a/2-y0,-a/2-y0]) 
    lam6,lam5 = 1/omega[2] * np.array([b/2-z0,-b/2-z0]) 
    LAM = [lam1,lam2,lam3,lam4,lam5,lam6]
    #LAM = [lam for lam in LAM if lam>0]
    lam_min = min([num for num in LAM if num>0],default=None)
    ind_ploskev = [0,1,2,3,4,5][LAM.index(lam_min)]
    #return omega
    #return LAM,lam_min, ind_ploskev
    return lam_min, ind_ploskev

#
   # 2. tale je dejansk bolj komplicirana ... ampak dela hitrejš!!
def PathToExit2(d,a,b,r0,direction):
    x0,y0,z0 = r0
    if x0<0 or x0>d or y0<-a/2 or y0>a/2 or z0<-b/2 or z0>b/2:
        return 'Wrong R0'
    ploskve = []
    
        # neposrečen poskus  ... ne dela v robnih primerih???
    theta, phi = direction
    min_value = 10**(-12)
    omega = (np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi))
    omega = [np.sign(value)*min_value if abs(value)<min_value else value for value in omega]
    #print(omega)
    #return omega
    omega = Spherical_to_cartesian(direction)

    if theta < np.pi/2 or theta > 3/2 * np.pi:
        lam_x = 1/omega[0] * (d - x0)
        ploskve.append(1)
    else:
        lam_x = 1/omega[0] * (-x0)
        ploskve.append(0)
    if theta < np.pi:
        lam_y = 1/omega[1] * (a/2 - y0)
        ploskve.append(3)
    else:
        lam_y = 1/omega[1] * (-a/2 - y0)
        ploskve.append(2)
    if phi < np.pi/2:
        lam_z = 1/omega[2] * (a/2 - z0)
        ploskve.append(5)
    else:
        lam_z = 1/omega[2] *(-a/2 - z0)
        ploskve.append(4)
    LAM = [lam_x,lam_y,lam_z]
    lam_min = min(LAM)
    ind_lam = LAM.index(lam_min)
    ind_ploskev = ploskve[ind_lam]
    #return ploskve,LAM,lam_min, ind_ploskev
    return lam_min, ind_ploskev
    

# začetni pogoj - vstop fotona v svinec / poljubna smer
  # ta ne dela za 2. fjo
R0 = [0.0000,0,0]
dir = [0,np.pi/2]
   # to dela za obe!!
R0 = [0.0001,0,0]
dir = [np.pi/2+0.001,0.00001]
  # dodatno
#R0 = [2,4,-3]
#dir = [0.6,0.7]
print(f'fja 1: {PathToExit1(10,20,20,R0,dir)}')
print(f'fja 2: {PathToExit2(10,20,20,R0,dir)}')


  # poljubni test ....

#%% Test: ker je bojši: PathToExit al PathToExit2 ??
a = b = 20
d = 10
N_photons = 10**4
np.random.seed(124782)
PhotonLocations = [[d*UNFM(),a/2*UNFM(-1,1),a/2*UNFM(-1,1)] for i in range(N_photons)]
PhotonDirections = [[UNFM(0,2*np.pi),UNFM(0,np.pi)] for i in range(N_photons)]
#PhotonLocations[1]

Ts1 = time.time()
data1 = [PathToExit1(d,a,a,r0,dir) for r0,dir in zip(PhotonLocations,PhotonDirections)]

ind = 7
print(data1[ind])

Ts2 = time.time()
data2 = [PathToExit2(d,a,a,r0,dir) for r0,dir in zip(PhotonLocations,PhotonDirections)]
print(data2[ind])

  # check če dela ? --> definitivno dela!!
print(data1 == data2, ' (,,BISTVEN PODATEK!!!"")')

  # časovna primerjava
print(f'T1 ={Ts2-Ts1} \nT2= {time.time() - Ts2}')
print('\n\n\n-------------dejansk zgleda druga koda boljša!!!----------------')

#%% test enakosti...
ind  = 10
loc = PhotonLocations[ind]
dir = PhotonDirections[ind]

theta, phi = dir
#omega = (np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi))
#print(omega)

print(PathToExit1(10,20,20,loc,dir))
print(PathToExit2(10,20,20,loc,dir))


#%% ChatGPT - ne zna bolš!! (...??)
def PathToExit(d, a, b, r0, direction):
    x0, y0, z0 = r0
    if x0<0 or x0>d or y0<-a/2 or y0>a/2 or z0<-b/2 or z0>b/2:
        return 'Wrong R0'
    theta, phi = direction
    omega = (np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi))
    # Initialize intersection distances with infinite values
    lam_x = lam_y = lam_z = np.inf
    # Calculate intersection distances with each face
    if omega[0] > 0:
        lam_x = (d - x0) / omega[0]
    elif omega[0] < 0:
        lam_x = -x0 / omega[0]
    if omega[1] > 0:
        lam_y = (a / 2 - y0) / omega[1]
    elif omega[1] < 0:
        lam_y = (-a / 2 - y0) / omega[1]
    if omega[2] > 0:
        lam_z = (b / 2 - z0) / omega[2]
    elif omega[2] < 0:
        lam_z = (-b / 2 - z0) / omega[2]
        
    # Find the minimum positive intersection distance and corresponding face index
    lam_min = min(lam_x, lam_y, lam_z)
    ind_ploskev = [1, 2, 3, 4, 5, 6][np.argmin([lam_x, lam_y, lam_z])]

    return lam_min, ind_ploskev


R0 = [4,5,7]
dir = [2.5,0.9]
PathToExit(10,20,20,R0,dir)
#print(1/np.cos(2.5)/np.sin(0.9+))



#%% sprememba smeri - rotacije
# Function to rotate a vector around the z-axis
def rotate_z(vec, angle):
    rot_matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                           [np.sin(angle), np.cos(angle), 0],
                           [0, 0, 1]])
    return np.dot(rot_matrix, vec)
# Function to rotate a vector around the y-axis
def rotate_y(vec, angle):
    rot_matrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                           [0, 1, 0],
                           [-np.sin(angle), 0, np.cos(angle)]])
    return np.dot(rot_matrix, vec)
# Function to rotate a vector around the x-axis
def rotate_x(vec, angle):
    rot_matrix = np.array([[1, 0, 0],
                           [0, np.cos(angle), -np.sin(angle)],
                           [np.sin(angle), np.cos(angle),0]])
    return np.dot(rot_matrix, vec)

def Rotate_zx(vec,az,ax):
    rot_z = np.array([[np.cos(az), -np.sin(az), 0],
                           [np.sin(az), np.cos(az), 0],
                           [0, 0, 1]]) 
    rot_x = np.array([[1, 0, 0],
                           [0, np.cos(ax), -np.sin(ax)],
                           [0,np.sin(ax), np.cos(ax)]])
    rot_matrix = np.dot(rot_x,rot_z)
    return np.dot(rot_matrix,vec)

def Rotate_yz(vec,th,phi):
    rot_y = np.array([[np.cos(th), 0, np.sin(th)],
                           [0, 1, 0],
                           [-np.sin(th), 0, np.cos(th)]])
    rot_z = np.array([[np.cos(phi), -np.sin(phi), 0],
                           [np.sin(phi), np.cos(phi), 0],
                           [0, 0, 1]])
    rot_matrix = np.dot(rot_z,rot_y)
    #return rot_matrix
    return np.dot(rot_matrix,vec)

Rotate_yz((1,0,0),0.9,0.5)

ex = (1,0,0)
vec1 = rotate_z(ex,1)
vec1
vec2=rotate_y(vec1,1)
vec2

 # primer
ex = (1,0,0)
th, phi = 0.7, 1.3
  # rotacija
print(Rotate_zx(ex,az=th,ax=phi))
  # ročno
x_new = np.cos(th)
y_new = np.sin(th)*np.cos(phi)
z_new = np.sin(phi)*np.sin(th)
print([x_new,y_new,z_new])

#array = np.array([[2,2,2],[1,3,1]])
#vec = np.array([-1,3,0]).T
#array
#np.dot(vec.T,array.T)
#np.dot(array,vec)

#%%  Rotacija S -> S'; ex = R*omega oz. pretvoriš v sistem, kjer lahko izvedeš Comptonsko rotacijo
def X_RotationMatrix(vec):
    # returns R = rotation matrix that rotates vector to ex vector: ex = R*vec
    # Normalize the given vector
    vec_norm = vec / np.linalg.norm(vec)
    # Calculate the axis of rotation
    axis = np.cross(vec_norm, [1, 0, 0])
    axis_norm = axis / np.linalg.norm(axis)
    # Calculate the angle of rotation
    angle = np.arccos(np.dot(vec_norm, [1, 0, 0]))
    # Construct the rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x, y, z = axis_norm
    R = np.array([[t*x**2 + c, t*x*y - z*s, t*x*z + y*s],
                  [t*x*y + z*s, t*y**2 + c, t*y*z - x*s],
                  [t*x*z - y*s, t*y*z + x*s, t*z**2 + c]])
    return R

# Example usage:
#vec = np.array([0, 1, 1])  # Example input vector
vec = np.array([-1, 2, 1])  # Example input vector
#vec = np.array([4, 1, 0])  # Example input vector
R = X_RotationMatrix(vec)
#print("Rotation Matrix:")
#print(R)

  # check da dela zadeva; vec = R^-1 * ex ...
print(np.dot(R,vec))
R_inv = np.linalg.inv(R) 
print(np.dot(R_inv,ex))


#%% Direction after Compton
def DirectionAfterCompton(direction,th_cs,phi_cs):
    # returns 
    th, phi = direction
    omega = np.array([np.sin(phi)*np.cos(th),np.sin(phi)*np.sin(th),np.cos(phi)])
    R = X_RotationMatrix(omega)
    R_inv = np.linalg.inv(R)
    omega2 = np.array([np.cos(th_cs),np.sin(th_cs)*np.cos(phi_cs), np.sin(th_cs)*np.sin(phi_cs)])
    new_direction =  np.dot(R_inv,omega2)
    #return new_direction
    return new_direction/np.linalg.norm(new_direction)
    
    ## TEST
DirectionAfterCompton(direction=(0,np.pi/2),th_cs=np.pi/4,phi_cs=0) # ex -> sipanje z kot pi/4=45° 
#DirectionAfterCompton(direction=(np.pi*3/2,np.pi/2),th_cs=np.pi/4,phi_cs=0) 

#%% kako zarotirat, če maš splošen omega?
def Rotate_vector(vector, axis, angle):
    # Normalize the axis vector
    axis = axis / np.linalg.norm(axis)
    
    # Construct the rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    # os vrtenja
    x, y, z = axis
    rotation_matrix = np.array([[t*x*x + c, t*x*y - s*z, t*x*z + s*y],
                                [t*x*y + s*z, t*y*y + c, t*y*z - s*x],
                                [t*x*z - s*y, t*y*z + s*x, t*z*z + c]])
    
    # Apply the rotation to the vector
    return np.dot(rotation_matrix, vector)

# Example usage
vector = np.array([0, 1, 0])  # Example vector
axis = np.array([0, 0, 1])    # Example axis (z-axis)
angle = np.pi / 2             # Example angle (90 degrees)
rotated_vector = Rotate_vector(vector, axis, angle)
print(rotated_vector)


#%%  compton rotation (??)
def ComptonRotation(direction,th_sc,phi_sc):
    th, phi = direction
    omega = np.array([np.sin(phi)*np.cos(th),np.sin(phi)*np.sin(th),np.cos(phi)])
    z_perp = 1
    r_perp = np.array([1,1,z_perp])
    x,y,z = omega
    return omega   

number = [UNFM(-10,10) for _ in range(10**5)]
plt.hist(number,bins=100)

#%%  perpendicular vector 1
def PerpendicularVector(vec):
    #vec2 = [UNFM(-1,1),UNFM(-1,1),UNFM(-1,1)]
    vec2 = np.array([1,1,1])
    return np.cross(vec,vec2)

def VectorAngle(vec):
    return np.arctan(vec[1]/vec[0])

ez = np.array([0,0,1])
PerpVecs = [PerpendicularVector(ez) for i in range(10**3)]
Angles = [VectorAngle(vec) for vec in PerpVecs]

PerpVecs[:20]
plt.hist(Angles,bins=100)

#%% perpendicular vector 2
def PerpendicularVector(vec):
    # Generate a random angle
    phi = np.random.uniform(0.0001, 2*np.pi-0.0001)
    
    # Rotate vec around an axis orthogonal to it
    axis = np.cross(vec, [1, 0, 0])  # Choose any vector not parallel to vec
    axis = axis /  np.linalg.norm(axis)  # Normalize the axis vector
    rot_matrix = np.array([[np.cos(phi), -np.sin(phi), 0],
                           [np.sin(phi), np.cos(phi), 0],
                           [0, 0, 1]])  # Rotation matrix for 2D rotation in xy-plane
    #axis = np.dot(rot_matrix,axis)
    return np.dot(rot_matrix, vec)
def VectorAngle(vec):
    if vec[0] == 0:
        vec[0] += 0.0001
    #else:
     #   vec[0] += 0.0001
    return np.arctan(vec[1] / vec[0])

ez = np.array([0, 0, 1])
PerpVecs = [PerpendicularVector(ez) for _ in range(10**3)]
Angles = [VectorAngle(vec) for vec in PerpVecs]

  #
#PerpVecs[:10]
#plt.hist(Angles,bins=20)


#%%  simulacija propagacije fotonov skozi svinec

def PhotonSimulation(E0):
    r0 = (0,0,0)
    omega = (1,0,0)
    E_out = np.zeros(6)
    E_current = E0
    #ac_fe, ac_compton = 
    while True:
        E_current
        U = [- np.log(UNFM()), - np.log(UNFM())]
        s_ints = min([U[i]/mu[i] for i in range(2)])
        lam, ploskev = PathToExit2(d,a,a,r0,direction=(th,phi))
        if s_int > lam:
            E_out[ploskev] += E_current
        #else:
        #    return None
