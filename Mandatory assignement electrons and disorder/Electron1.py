import pandas as pd
from pandas import DataFrame
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import numpy as np
import scipy.constants as scic
import scipy.integrate as scint
import scipy.optimize as scis
from scipy.integrate import *
from scipy.constants import *
from scipy.optimize import *
from sympy.solvers import solve
from sympy import Symbol, diff, exp, Function, lambdify, Eq , rsolve, log , tanh 

#################
#Define constants
#################
#we put kb=1 and Tf=1 since everything is normalised
#V=1 hbar=1 m=1  

#########
#Define functions
#########
def f_FD (energy , beta , mu):
    return 1/(np.exp((energy-mu)*beta)+1)
def Dos (energy):
    return 2*np.sqrt(2)*pow(energy,1/2)/(2*pow(scic.pi,2))
def g(mu, beta):
    domain= np.linspace(0,100,1000) 
    dx=domain[1]-domain[0]
    integral=0
    n0 = 2*np.sqrt(2)/(3*pow(scic.pi,2))
    for en in domain:
        integral+= Dos(en)*f_FD(en, beta ,mu )*dx 
    return integral-n0


T=np.linspace(0.01,1.3,100)
Beta=[1/val for val in T]
######################
#Temperature dependence of mu
######################
#numeric evaluation of integral rectangular:
root=[]
"""
for i in Beta:
    root.append(scis.fsolve(g, 1 , args=(i)))
#print ("radice: ", root)
plt.plot(T,root)
plt.xlabel(r"$\frac{T}{T_f} $")
plt.ylabel(r"$\frac{\mu}{\epsilon_f}$")
plt.grid()
plt.show()
"""

T=[0.01,0.1,0.5,1,1.2]
Beta=[1/val for val in T]
root1=[]
for i in Beta:
    root1.append(scis.fsolve(g, 1 , args=(i)))
#Point b)
x=np.linspace(0,4,1000)
for i in range (0,len(T)):
    plt.plot(x,f_FD(x,Beta[i],1))
    #plt.plot(x,f_FD(x,Beta[i], root[i]),'--')
plt.xlabel(r"$\frac{\epsilon}{\epsilon_f} $")
plt.ylabel(r"$f_{FD}(\frac{\epsilon}{\epsilon_f})$")
plt.plot(x,0.5*np.ones(len(x)),':', color='grey')
plt.legend(['0.01'r"$T_f$",'0.1'r"$T_f$",'0.5'r"$T_f$",'1'r"$T_f$",'1.2'r"$T_f$"])
plt.grid()
plt.show()
#Point d)
x=np.linspace(0,4,1000)
for i in range (0,len(T)):
    #plt.plot(x,f_FD(x,Beta[i],1))
    plt.plot(x,f_FD(x,Beta[i], root1[i]))
plt.xlabel(r"$\frac{\epsilon}{\epsilon_f} $")
plt.ylabel(r"$f_{FD}(\frac{\epsilon}{\epsilon_f})$")
plt.plot(x,0.5*np.ones(len(x)),':', color='grey')
plt.legend(['0.01'r"$T_f$",'0.1'r"$T_f$",'0.5'r"$T_f$",'1'r"$T_f$",'1.2'r"$T_f$"])
plt.grid()
plt.show()

#comparison about point b and d
x=np.linspace(0,4,1000)
for i in range (0,len(T)):
    plt.plot(x,f_FD(x,Beta[i],1), label= str(T[i]) + str(r"$T_f$"))
    plt.plot(x,f_FD(x,Beta[i], root1[i]),'--', label=str(T[i]) + str(r"$T_f$"))
plt.xlabel(r"$\frac{\epsilon}{\epsilon_f} $")
plt.ylabel(r"$f(\epsilon)$")
plt.plot(x,0.5*np.ones(len(x)),':', color='grey')
plt.legend()
plt.grid()
plt.show()