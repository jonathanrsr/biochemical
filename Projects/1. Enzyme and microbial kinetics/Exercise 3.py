# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:33:46 2024

@author: pauli
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from scipy.optimize import curve_fit

def model(s, vmax, km):
    return vmax*(s/(km+s))

x= np.array([0.68, 0.16, 0.006])
y= np.array([0.16,0.13,0.0385])
s = np.linspace(0,1,1000)

popt_1, pcov_1 = curve_fit(model, x, y)
vmax, km = popt_1
print(popt_1)

plt.plot(s, model(s, popt_1[0], popt_1[1]), label="Optimized model", c='dodgerblue')
plt.scatter(x,y,label='Experimental data',c='black',s=20)
plt.axhline(y=vmax, linestyle=':', label = '$v_{max}$ = 0.158')
plt.plot([-0.05, km ], [vmax/2,vmax/2], linestyle='--', c='darkmagenta')
plt.plot([km, km ], [0,vmax/2], linestyle='--', c='darkmagenta')
plt.scatter(km,vmax/2, c='darkmagenta', s=20)
plt.annotate('$K_M$= 0.022', xy=(km,0), xytext=(0.1, 0.01),
             arrowprops=dict(color="violet",
                             headwidth=5,
                             headlength=3,
                             width=0.7))

plt.xlim([-0.02,1])
plt.ylim([0,0.17])



plt.xlabel('P [mmol/L]')
plt.ylabel('$r_s$')
plt.grid()
plt.legend()
#plt.savefig('bonus_graphfinal.jpg')
