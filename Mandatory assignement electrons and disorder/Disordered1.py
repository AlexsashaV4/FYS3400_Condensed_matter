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

def tau(E,beta,k=1):
    return np.log(k*E**3*pow(np.tanh(beta*E/2),-1))
def C(E,beta):
    return pow(beta*E,2)*np.exp(beta*E)/((np.exp(beta*E)+1)**2)
t = Symbol('t')
f=diff(pow(t,2)*exp(t)/((exp(t)+1)**2))
g1 = diff(log(t**3*pow(tanh(10000*t/2),-1)))
g2 = diff(log(t**3*pow(tanh(100*t/2),-1)))
g3 = diff(log(t**3*pow(tanh(1*t/2),-1)))

fx=lambdify(t,f,modules=['numpy']) 
g1x=lambdify(t,g1) 
g2x=lambdify(t,g2) 
g3x=lambdify(t,g3) 
#t=np.linspace(0,100,1000)
#f1=t**2*exp(t)/(exp(t) + 1)**2 - 2*t**2*exp(2*t)/(exp(t) + 1)**3 + 2*t*exp(t)/(exp(t) + 1)**2
x=np.linspace(0.0001,100,10000000)
y=fx(x)
T=[0.001,0.1,1,10]
######
#Plot C and tau in the kbT range and compare#
######
for i in T: 
    fig=plt.figure()
    x=np.linspace(0.0001,10*i,10000000)
    idx=np.where(abs(C(x,1/i) - max(C(x,1/i))/2)<10**-3)[0][0]
    idxm=np.where(C(x,1/i)==max(C(x,1/i)))[0][0]
    idx1=np.where(abs(C(x[idxm:],1/i) - C(x[idx],1/i))<10**-3)[0][0]
    #plt.plot(x[idxm:],C(x[idxm:],1/i),'--',label=i)
    print("idx ", idx ,"C", (C(x[idx],1/i)),"\n")
    print("idxm", idxm, "\n")
    print("idx1 ", idx1,"\n")
    val=x[idx]
    val1=x[idx1+idxm]
    #plt.plot([val,val],[0,C(val,1/i)],'--',label=i)
    #plt.plot([val1,val1],[0,C(val1,1/i)],'--',label=i)
    #plt.axvspan(val, val1, facecolor='g', alpha=0.5)
    ax1 = fig.add_subplot(111)
    ax1.plot(x,C(x,1/i),label=str(i)+"K")
    ax1.plot([val,val],[0,C(val,1/i)],'--',color='orange')
    ax1.plot([val1,val1],[0,C(val1,1/i)],'--',color='orange')
    ax1.axvspan(val, val1, facecolor='g', alpha=0.3 )
    ax1.set_ylabel(r"$C_{TLS}(E)$")
    ax1.grid()
    ax1.set_xlabel(r"$E$")
    ax2 = ax1.twinx()
    ax2.plot(x, tau(x,1/i), 'r-')
    ax2.set_ylabel(r"$ln(\frac{1}{\tau(E)}) $", color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    #plt.plot([2.4*i,2.4*i],[0,tau(2.4*i,1/i)],':',label=i)
    plt.show()
