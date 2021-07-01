import numpy as np 
from matplotlib import pyplot as plt
from numpy.core.numeric import ones
import numpy as np
import scipy.constants as scic
import scipy.integrate as scint
import scipy.optimize as scis
from scipy.integrate import *
from scipy.constants import *
from scipy.optimize import *
from sympy import symbols, Eq, solve
from sympy.functions.elementary.miscellaneous import root

eV = 1.6e-019; 
kB = 1.38e-023; # Boltzmann constant, J/mol/K
h= 6.626e-034; # Planck's constant, Js
m_0 = 9.1e-031; # Rest mass of an electron, kg
#Parameters

meff_elect = 1.08*m_0; # Effective mass of electron
meff_hole = 0.56*m_0; # Effective mass of holes 
Eg = 1.1*eV; # Energy gap of Si #setting E_v=0
Nd = 1e17
deltaE = 0.045*eV 
Ed=(Eg-deltaE) #Energy of the donor band

Efi= lambda T: Eg/2 + (3/4)*kB*T*np.log(meff_hole/meff_elect) #intrinsic fermi energi
ni= lambda T: 2*(2*np.pi*kB*T/h**2)**(3/2)*(meff_elect*meff_hole)**(3/4)*np.exp(-Eg/(2*kB*T)) / 10**6
Nc = lambda T: 2*(2*np.pi*kB*T/h**2)**(3/2)*(meff_elect)**(3/2)/ 10**6
Nv= lambda T: 2*(2*np.pi*kB*T/h**2)**(3/2)*(meff_hole)**(3/2)/ 10**6

#functions
def n(Ef,T): #carrier distribution in conduction
    return ni(T)*np.exp((Ef-Efi(T))/(kB*T))
def p(Ef,T): #holes distribution in valence
    return ni(T)*np.exp((Efi(T)-Ef)/(kB*T))
def Ndp(Ef,T): # doped atoms distribution
    return Nd*(1-1/(1+0.5*np.exp((Ed-Ef)/(kB*T))))
def g(Ef,T): # neutral charge equation
    return  Ndp(Ef,T) + p(Ef,T) - n(Ef,T)
def nc(Ef,T): 
    return Nc(T)*np.exp((Ef-Eg)/(kB*T))
def pv(Ef,T): 
    return Nv(T)*np.exp((-Ef)/(kB*T))

Tvals = np.linspace(60,1500,1500)
Tvals1 = [1000/T for T in Tvals]
root1=np.ones(len(Tvals))
c=0
for i in Tvals:
    root1[c]=scis.fsolve(g, Eg , args=(i),xtol=eV**2)
    print(root1[c])
    c +=1 
eg= np.linspace(0,1.1,len(Tvals))
# plt.plot(Tvals1,root1/eV, label=r"$E_F$ Fermi energy")
# plt.plot(Tvals1,Efi(Tvals)/eV, label=r"$E_{F_i}$ Intrinsic Fermi energy")
# plt.plot(Tvals1,Eg*np.ones(len(Tvals))/eV,'--',label=r"$E_g$ Energy gap" )
# plt.plot(Tvals1,Ed*np.ones(len(Tvals))/eV,'--',label=r"$E_D$ Donor's energy band")
# plt.plot(Tvals1,0*np.ones(len(Tvals)),'--',label=r"$E_v$ Valence band")
# plt.xlabel(r"$1000/T [K^{-1}]$")
# plt.ylabel(r"$E [eV]$")
# plt.legend()
# plt.grid()
# plt.show()
#plt.semilogy(Tvals1, p(root1,Tvals),'#d62728', label=r"p hole carrier concentration in valence band")
plt.semilogy(Tvals1, n(root1,Tvals), label=r"n carrier concentration in conduction band")
#adding the approximations 

plt.semilogy(Tvals1[900:], ni(Tvals[900:]), '--', label=r"High temperature approximation", linewidth=2)
#plt.semilogy([Tvals1[400],Tvals1[400]], [1e16,np.max(ni(Tvals))],'#ff7f0e', '--')
plt.semilogy(Tvals1[80:900], Nd*np.ones(len(Tvals1[80:900])), linestyle='dashed',color='#d62728',label=r"medium temperature approximation", linewidth=2)
#plt.semilogy([4,4], [np.min(ni(Tvals[500:])),np.max(n(root1,Tvals))],'#ff7f0e', '--')
###
#low temperatures
###
nstar= lambda T: Nc(T)/2*np.exp(-deltaE/(kB*T))
zztop= -nstar(Tvals[2:150])/2 + np.sqrt(0.25*nstar(Tvals[2:150])**2 + nstar(Tvals[2:150])*Nd)
#plt.semilogy(Tvals1[2:150], zztop, '--')
plt.semilogy(Tvals1[1:80], np.sqrt(0.5*Nd*Nc(Tvals[1:80]))*np.exp((Ed-Eg)/(2*kB*Tvals[1:80])), linestyle='dashed',color='#2ca02c',label=r"Low temperature approximation", linewidth=2)
plt.semilogy()
plt.grid(True, which="both", ls="--")
plt.xlabel(r"$1000/T [K^{-1}]$")
plt.ylabel(r"$ n [cm^{-3}]$")
plt.legend()
plt.show()
#####
# #Ef and concentration
# root2=np.ones(len(Tvals))
# c=0
# for rr in range (10,19):
#     for i in Tvals:
#         Nd=10**rr
#         root2[c]=scis.fsolve(g, Eg , args=(i),xtol=eV**2)
#         c +=1
#     plt.plot(Tvals, (root2-Efi(Tvals))/eV,'#d62728')
# plt.plot(Tvals,Eg/2*np.ones(len(Tvals))/eV,'--',label=r"$E_g$ Energy gap" )
# plt.plot(Tvals,-Eg/2*np.ones(len(Tvals))/eV,'--',label=r"$E_D$ Donor's energy band")
# #plt.plot(Tvals1,0*np.ones(len(Tvals)),'--',label=r"$E_v$ Valence band")
# plt.xlabel(r"$1000/T [K^{-1}]$")
# plt.ylabel(r"$E [eV]$")
# plt.show() 


#print(g(root1,Tvals))
#plt.semilogy(eg,n(root1,Tvals)+Ndp(root1,Tvals))



