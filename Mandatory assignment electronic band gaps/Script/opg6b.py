import numpy as np
import matplotlib.pyplot as plt
#from plot_set import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d


def E(kx, ky):
    # notice that it should really be +- t
    E = t*np.sqrt(1 + 4*np.cos(np.sqrt(3)*ky*a/2)*np.cos(3*kx*a/2) + 4*np.cos(np.sqrt(3)*ky*a/2)**2)
    return E

def Edx(kx, ky): #dE/dk_x
    sqrt3 = np.sqrt(3)
    return t**2/(2*E(kx,ky))*(-6*a*np.cos(sqrt3*ky*a/2)*np.sin(3*kx*a/2))

def Edy(kx, ky): #dE/dk_y
    sqrt3 = np.sqrt(3)
    return t**2/(2*E(kx,ky))*(-2*sqrt3*a*np.sin(sqrt3*ky*a/2)*np.cos(3*kx*a/2) - 4*sqrt3*a*np.cos(sqrt3*ky*a/2)*np.sin(sqrt3*ky*a/2))

def m_star(kx, ky): #effective mass in [01] direction (b1,b2 basis)
    d = np.array([1/3, -1/np.sqrt(3)]) # in kx,ky basis
    return (d[0]*kx+d[1]*ky)/(d[0]*Edx(kx,ky) + d[1]*Edy(kx,ky)) #[hbar^2]

def plot3D_energy():
    fig=plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111,projection='3d') 
    N = int(3e2+1)
    span = np.pi/a
    #span = 1/a
    x = np.linspace(-span,span, N)
    y = np.linspace(-span,span, N)
    X, Y = np.meshgrid(x, y)
    Z = E(X, Y)
    #ax = plt.axes(projection='3d')
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    #ax.axes.set_xlim3d(left=-1.0, right=1.0) 
    #ax.axes.set_ylim3d(bottom=-1.0, top=1.0)     
    ax.set_zlabel('E(kx,ky)')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.plot_surface(X, Y, -Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_title('Energy 3D')
    ax.grid()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    plt.show()


def heatmap(X,Y,Z):
    plt.pcolor(X,Y,Z)
    plt.colorbar()

def plot_energy(X,Y,Z):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surf=ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    ax.plot_surface(X, Y, -Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    ax.set_title('Graphene energy band')


def plot_01_direction(): 
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k') # [01]
    N = int(1e4+1)
    span  = 1
    k2 = np.linspace(-span, span, N)
    m = m_star(b2[0]*k2, b2[1]*k2)
    tol = 1000 # remove high peaks
    m_stripped = np.where(np.abs(m) > tol, np.nan, m) # remove high peaks

    plt.title("Effective mass along direction [01]")
    plt.plot(k2, m_stripped )#,label=r'Effective mass')
    plt.xlabel(r"$k$", fontsize=10)
    plt.ylabel(r"$m^*$", fontsize=10)
    plt.ylim([-250, 250])
    plt.plot(0,m_star(0.001,0.001),'ro', label=r"$m^*_{\Gamma} = 0.84 m_e$")
    plt.plot([0.5,0.5],[-1000,+1000],'#ff7f0e',linestyle='dashed')
    plt.plot([-0.5,-0.5],[-1000,+1000],'#ff7f0e',linestyle='dashed')
    plt.grid()
    plt.legend(loc=8)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2);

def m_star_gamma(): #find effective mass in limit of kx -> 0 with y = 0
    kx = 1e-20
    m = m_star(kx, 0) #[hbar^2 1/(nm^2*eV)]
    m_e = 9.1093837015e-31
    hbar = 6.582119569e-16 #[eV*s]

    #Convert
    m *= (1e9)**2   # --> [hbar^2 1/(m^2*eV)]
    m *= hbar**2    # --> [eV s**2/m**2]

    convert = 1.60217662e-19 # [eV s**2/m**2] --> [kg]
    m *= convert    # --> [kg]
    m /= m_e        # --> [m_e]
    print(m)


def plot_Bzone():
    plt.figure(num=1, dpi=80, facecolor="w", edgecolor="k") #[B-zone]
    N = int(1e2+1)

    span = 1/2*np.linalg.norm(b2)*3

    kx = np.linspace(-span,span, N) #B-zone?
    ky = np.linspace(-span,span, N) #B-zone?
    KX, KY = np.meshgrid(kx,ky)
    m3D = m_star(KX,KY)
    tol = 1000
    m3D_stipped = np.where(np.abs(m3D) > tol, np.nan, m3D) # remove high peaks

    plt.title("...")
    # plot3D(KX,KY,m3D_stipped)
    heatmap(KX, KY, m3D_stipped)
    plt.xlabel(r"$k_x$", fontsize=14)
    plt.ylabel(r"$k_y$", fontsize=14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)




if __name__ == "__main__":
    t = -3 # [eV]
    a = 0.142 # [nm]
    # a = 0.142e-9 #[m]


    b1 = 2*np.pi/a*np.array([1/3, +1/np.sqrt(3)])
    b2 = 2*np.pi/a*np.array([1/3, -1/np.sqrt(3)])
    N = int(1e2+1)
    span = 1/2*np.linalg.norm(b2)*3

    kx = np.linspace(-span,span, N) #B-zone?
    ky = np.linspace(-span,span, N) #B-zone?
    span = 2*np.pi
    x = np.linspace(-span/a,span/a, N)
    y = np.linspace(-span/a,span/a, N)
    X, Y = np.meshgrid(x, y)
    Z = E(X, Y)
    #plot_energy(X,Y,Z)
    plot_01_direction()
    plot3D_energy()
    m_star_gamma()
    #plot_Bzone()


    plt.show()
