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
    return vmax*(s/(km+s)) # Michaelis-Menten model 

t = np.array([0, 2, 6, 10]) # h
S = np.array([1, 0.68, 0.16, 0.006]) # Sucrose concentration, mmol/L
S_average = np.zeros(len(S) - 1) # Average between samples
rs = np.zeros(len(S) - 1) # reaction rate, mmol/L/h

#Calculate the reaction rate between 2 measurements (Euler's method) and average sucrose concentration
for x in range(len(S) - 1):
    rs[x] = -(S[x + 1] - S[x])/(t[x + 1] - t[x])
    S_average[x] = (S[x + 1] + S[x])/2

s = np.linspace(0, 2, 1000)
popt_1, pcov_1 = curve_fit(model, S_average, rs)  #Call of the function that will fit the model to the experimental data given
vmax, km = popt_1 #Parameters we are trying to optimize 

print(popt_1)
print(S_average)
print(rs)

plt.scatter(S_average, rs, label = 'data', color = 'black', s = 4.0)
plt.plot(s, model(s, popt_1[0], popt_1[1]), label = "model", color = 'black', linewidth = 1.0)
plt.axhline(y = vmax, color = 'black', linestyle=':', linewidth = 1.0)
plt.plot([0, km], [vmax/2, vmax/2], linestyle=':', linewidth = 1.0, color = 'black')
plt.plot([km, km], [0, vmax/2], linestyle=':', linewidth = 1.0, color = 'black')
plt.annotate(r'$\text{K}_{m}$', xy = (km, 0), xytext = (km*1.2, 0.02),
             arrowprops = dict(arrowstyle = '-|>', color = 'black'))
plt.annotate(r'$\nu_{max}$', xy = (0, vmax), xytext = (s[-1]/5, vmax*0.95))
plt.annotate(r'$\frac{\nu_{max}}{2}$', xy = (0, vmax/2), xytext = (km/5, vmax/2*1.05))
plt.xlabel(r'$\text{[S] (mmol.L}^{-1}\text{)}$')
plt.ylabel(r'$\text{r}_{s}$')
plt.xlim(0, s[-1])
plt.grid(linewidth = 0.25)
plt.legend()
plt.title("Michaelis-Menten plot")
plt.show()
#plt.savefig('bonus_graphfinal.jpg')