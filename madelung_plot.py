from matplotlib import pyplot as plt # import the plotting tool "pyplot" and use the keyword "plt" to call it
import numpy as np # import numerical library for vectors and use keyword np
import pandas as pd # import pandas, useful library, with keyword pd
import scipy.constants as sc
from sympy import symbols, Eq, solve
# Plot 1
def Pauli(r,lam,ro):
    return  lam*np.exp(-r/ro)
def elPot(r,e):
    return pow(e,2)/(4*np.pi*sc.epsilon_0*r)
N = np.linspace(0, 50, 50) # creates linearly spaced vector from 0 to 50, 50 points
s=0
mad = [0]
q = sc.e
lam =(1.0e3)*q
ro =0.32e-10
a =2.82e-10
eps = sc.epsilon_0
energy = [0]
z = 2
madth=2*np.log(2)
Utot = [0]
Utot_th = [0]
# now calculate the madelung constant no boundary conditions
for i in range (1,50,1):
    s=s+(pow(-1,i-1)/i)
    mad.append(2*s)
    energy.append(z*Pauli(a,lam,ro)-2*s*elPot(a,q))
    Utot.append(i*energy[i])
    Utot_th.append(i*(z*Pauli(a,lam,ro)-2*np.log(2)*elPot(a,q)))
#print(Utot)
#plot madelung
plt.plot(N, mad, color='blue', linewidth=2.2)
plt.title('Madelung constant as a function of N')
plt.xlabel('N')
plt.ylabel('Madelung constant')
y= np.ones(50)
plt.plot(N,2*np.log(2)*y)
plt.grid(True)
plt.show()


#plot energy as a function of N 
plt.plot(N,Utot_th,color='orange')
plt.plot(N,Utot)

plt.title('U as a function of N')
plt.xlabel('N')
plt.ylabel('U [J]')
plt.grid(True)
plt.show()

#plot energy as a function of radius
x1=np.linspace(1e-10,8e-10,100)
for i in range(10,len(N)+1,10):
    plt.plot(x1,N[i-1]*(z*Pauli(x1,lam,ro)-madth*elPot(x1,q)),label='N=%d' %i)
plt.grid(True)
plt.title('U as a function of r and N')
plt.xlabel('r [m]')
plt.ylabel('U [J]')
plt.legend()
plt.show()

#To find the minimum I'm going to implement newton
r=0.5e-10
for l in range(1,100,1):
    der1=-(z/ro)*Pauli(r,lam,ro) + (madth*elPot(r,q)/r)
    der2=(z/pow(ro,2))*Pauli(r,lam,ro) - (2*madth*elPot(r,q)/pow(r,2))
    r=r-der1/der2
print("r:",r)