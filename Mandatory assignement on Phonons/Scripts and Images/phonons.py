from matplotlib import pyplot as plt # import the plotting tool "pyplot" and use the keyword "plt" to call it
import numpy as np # import numerical library for vectors and use keyword np
import pandas as pd # import pandas, useful library, with keyword pd
import scipy.constants as sc
import numpy as np 
from matplotlib import pyplot as plt 
import matplotlib
from sympy import symbols, Eq, solve

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

w = np.linspace(0, 0.95, 10000)
plt.plot(w, np.full(10000, 1.5), linewidth=1.8, color='violet', label='DOS(k) = $L/2\pi$ ')
plt.grid()
ax = plt.gca()
ax.axes.xaxis.set_ticklabels([], fontsize=14)
ax.axes.yaxis.set_ticklabels([], fontsize=14)
plt.title('Density of states', fontsize=14)
plt.legend()
plt.xlabel('k')
plt.ylabel('DOS(k)')
plt.show()
plt.plot(w, 1/np.sqrt(1-w**2), linewidth=1.8, color='orange', label=r'DOS($\omega$)')
plt.title('Density of states', fontsize=14)
plt.legend()
plt.xlabel('$\omega$')
plt.ylabel('DOS($\omega$)')
plt.grid()
plt.show()