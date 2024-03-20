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
    return vmax*(s/(km+s)) #Michaelis-Menten model 

x= np.array([0.84, 0.42, 0.083]) #Sucrose concentration  
y= np.array([0.16,0.13,0.0385]) #Rates derived from the product concentration
s = np.linspace(0,1,1000)

popt_1, pcov_1 = curve_fit(model, x, y)  #Call of the function that will fit the model to the experimental data given
vmax, km = popt_1 #Parameters we are trying to optimize 
print(popt_1)

plt.rcParams['font.family'] = 'serif'
plt.plot(s, model(s, popt_1[0], popt_1[1]), label="Optimized model", c='dodgerblue')
plt.scatter(x,y,label='Experimental data',c='black',s=20)
plt.axhline(y=vmax, linestyle=':', label = '$v_{max}$ = 0.158 mmol/L/hr')
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

#--------------- LINEWEAVER BURK PLOT -----------------------------

Substr= np.array([0.84, 0.42, 0.083])

rs= np.array([0.16,0.13,0.0385])

s = np.linspace(0,1,1000)

plt.scatter(1/Substr, 1/rs)
plt.ylim([0,30])

def model_lin_LB(x,a,b):
    return a*x+b

popt_2, pcov_2 = curve_fit(model_lin_LB, (1/Substr), (1/rs))
A, B = popt_2
print(popt_2)

overS = np.linspace(-80,200, 500)
overRS= np.linspace(-5, 30, 500)

fig, ax=plt.subplots(figsize=(9,7))
ax.axhline(y=0, c='black')
ax.axvline(x=0, c='black')
ax.plot(overS, model_lin_LB(overS, A, B), linestyle = '--')
ax.scatter(1/Substr, 1/rs)
ax.set_ylabel('$1/r_S$')
ax.set_xlabel('$1/S$')
ax.scatter(-B/A, 0, label=f'-1/$K_M$ = {-B/A:.2f}')
ax.scatter(0, B, label='$1/μ_{max}$= 6.51', c='r')
ax.legend()
ax.grid()
#fig.savefig('LBP_Bonus.jpg')
