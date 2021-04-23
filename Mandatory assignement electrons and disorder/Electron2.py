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

def KPM_Soln(a,b,U0,eps_range):  #input a, b spacing in Angstoms, potential in eV, and total desired range to output
    #Constants
    h_bar = 1.054*1e-34    #J*s
    m = 9.109*1e-31        #kg
    alpha_0 = (2*m*U0*1.602*1e-19/h_bar**2)**(1/2)  #m^-1

    #Kronig-Penny Solution (LHS),  with epsilon = E/U0
    def KPM_p(eps):  #for epsilon > 1
        return (1-2*eps)/(2*(eps*(eps-1))**(1/2))*np.sin(alpha_0*a*1e-10*eps**(1/2))*np.sin(alpha_0*b*1e-10*(eps-1)**(1/2))+np.cos(alpha_0*a*1e-10*eps**(1/2))*np.cos(alpha_0*b*1e-10*(eps-1)**(1/2))
    def KPM_m(eps):  #for epison < 1
        return (1-2*eps)/(2*(eps*(1-eps))**(1/2))*np.sin(alpha_0*a*1e-10*eps**(1/2))*np.sinh(alpha_0*b*1e-10*(1-eps)**(1/2))+np.cos(alpha_0*a*1e-10*eps**(1/2))*np.cosh(alpha_0*b*1e-10*(1-eps)**(1/2))

    #Define epsilon space to plot
    epslist = np.linspace(0,eps_range,200000)
    f_eps = np.piecewise(epslist, [epslist < 1, epslist > 1], [KPM_m, KPM_p])
    
    return epslist, f_eps

def Eband_KP(epslist,f_eps): #outputs energy band data
    k=[]
    bandlist=[]
    Eps=[]
    epsbuildlist=[]
    for i in range(len(f_eps)-1):
       if 1 >= f_eps[i] >= -1:
           bandlist.append(f_eps[i])
           epsbuildlist.append(epslist[i])
           if (1 < f_eps[i+1] or  -1 > f_eps[i+1]):
               k.append(bandlist)
               Eps.append(epsbuildlist)
               bandlist=[]
               epsbuildlist=[]

    for i in range(len(k)):
        k[i]=np.arccos(k[i])/np.pi
        if i % 2 == 0:
            Eps[i]=np.concatenate((Eps[i][::-1],Eps[i][::1]))
            k[i]=np.concatenate((-1*np.array(k[i],dtype=float)[::-1],k[i][::1]))
        else:
            k[i]=np.concatenate((k[i][::1],-1*np.array(k[i],dtype=float)[::-1]))
            Eps[i]=np.concatenate((Eps[i][::1],Eps[i][::-1]))
                                                      
    return k, Eps


xmin=-4*scic.pi
xmax=4*scic.pi
x=np.linspace(xmin, xmax,10000)
V=5
f=V*np.sin(x)/x +np.cos(x)
#fig=plt.figure()
#ax1=fig.add_subplot(111)
plt.plot(x,f)
plt.plot([xmin, xmax],[1,1],color= 'orange')
plt.plot([xmin, xmax],[-1,-1],color= 'orange')
plt.ylabel(r"$f(ka)$")
plt.xlabel(r"$ka$")
plt.xticks(np.arange(xmin, xmax+scic.pi, step=scic.pi),['-4'r"$\pi$", '-3'r"$\pi$", '-2'r"$\pi$",r"$\pi$",'0',r"$\pi$",'2'r"$\pi$",'3'r"$\pi$",'4'r"$\pi$"])
y1=f*0+1
y2=f*0-1
plt.fill_between(x, f, y1, where=f >= y1,
                 facecolor="orange", # The fill color
                 color='blue',       # The outline color
                 alpha=0.2)          # Transparency of the fill
plt.fill_between(x, f, y2, where=f <= y2,
                 facecolor="orange", # The fill color
                 color='blue',       # The outline color
                 alpha=0.2)  

# Show the plot
plt.grid()
plt.show()
def ff(c):
    return V*np.sin(c)/c +np.cos(c)
print(np.arange(xmin, xmax+scic.pi, step=scic.pi))
b=[]
en1=[]
for i in np.arange(xmin, xmax+scic.pi, step=scic.pi):

    k=np.linspace(i+0.0001,i+scic.pi,1000)
    xband=np.where(abs(ff(k))<1)[0]
    f=ff(k[xband])
    energy= ((abs(np.arccos(f))+i)**2)
    plt.plot(k[xband],energy)
#plt.axvspan(x[xband],x[xband], facecolor='g', alpha=0.5)
#plt.plot([xmin,xmax],[1,1],color= 'red')
#plt.plot([xmin,xmax],[-1,-1],color= 'red')
#plt.savefig('David_gorgonzola.png')  
plt.grid()
plt.show()


#Ec, Ev = eps[2], eps[1]
#Eg = (np.min(Ec)-np.max(Ev))*V

