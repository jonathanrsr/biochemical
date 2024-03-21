# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:33:46 2024

@author: pauli
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

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

#--------------- LINEWEAVER BURK PLOT -----------------------------

# Linear regression of 1/nu versus 1/S
x = (1/S_average).reshape((-1, 1))
y = 1/rs
model2 = LinearRegression().fit(x, y)
y_intercept = float(model2.intercept_)
slope = model2.coef_[0]
x_intercept = -y_intercept/slope
det_coefficient = model2.score(x, y)

print(f'Slope: {round(slope, 2)}\nIntercept (x): {round(x_intercept, 2)}\nIntercept (y): {round(y_intercept, 2)}')

# Find kinetics parameter from the linear regression
Km = -1/x_intercept
nu_max = 1/y_intercept

# Create 2 points for plotting the model
x_graph = np.array([x_intercept, 1/S_average[-1]*1.2]) # Multiply by 1.2 to make the line go further than the last data point
y_graph = model2.predict(x_graph.reshape((-1, 1)))

# Lineweaver-Burke type plot
plt.plot(1/S_average, 1/rs, color = 'black', marker = 'o', markersize = '4', linestyle = 'None', label = 'data')
plt.plot(x_graph, y_graph, color = 'black', marker = 'None', linestyle = '--', linewidth=1.0, label = r'$\text{model, R}^{2} = ' + r'\text{' + f'{round(det_coefficient, 3):.3f}' + r'}$')
plt.annotate(r'$\frac{1}{\nu_{\mathrm{max}}}$', xy = (0, y_intercept), xycoords = 'data', xytext = (50, 0), textcoords = 'offset points', arrowprops = dict(arrowstyle = '-|>', color = 'black'))
plt.annotate(r'-$\frac{1}{K_{\mathrm{M}}}$', xy = (x_intercept, 0), xycoords = 'data', xytext = (-6, 50), textcoords = 'offset points', arrowprops = dict(arrowstyle = '-|>', color = 'black'))
plt.xlabel(r'$\text{1/S (L.mmol}^{-1}\text{)}$')
plt.ylabel(r'$\text{1/r}_{s}\text{ (L.s.mmol}^{-1}\text{)}$')
plt.vlines(0, 0, np.max(y_graph)*1.05, color = 'black', linestyle = '-', linewidth = 1.0) # Multiply by 1.05 for better design
plt.ylim([0, np.max(y_graph)*1.05]) # Multiply by 1.05 for better design
plt.grid(which = 'both', linewidth = 0.25)
plt.legend(loc = 'lower right')
plt.title('Lineweaver-Burke type plot')
plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\lineweaver_burke_type_plot_3.eps', format = 'eps')
plt.show()

print(f'Km = {round(Km, 2)}\nmax = {round(nu_max, 2)}')

#--------------- MICHAELIS-MENTEN PLOT -----------------------------
S_model = np.linspace(0, 2, 1000)
rs_model = nu_max*(S_model/(Km + S_model))

plt.scatter(S_average, rs, label = 'data', color = 'black', s = 10.0)
plt.plot(S_model, rs_model, label = "model", color = 'black', linestyle = '--', linewidth = 1.0)
plt.axhline(y = nu_max, color = 'black', linestyle = ':', linewidth = 1.0)
plt.plot([0, Km], [nu_max/2, nu_max/2], linestyle = ':', linewidth = 1.0, color = 'black')
plt.plot([Km, Km], [0, nu_max/2], linestyle = ':', linewidth = 1.0, color = 'black')
plt.annotate(r'$\text{K}_{M}$', xy = (Km, 0), xytext = (Km*1.2, 0.02),
             arrowprops = dict(arrowstyle = '-|>', color = 'black'))
plt.annotate(r'$\nu_{max}$', xy = (0, nu_max), xytext = (S_model[-1]/5, nu_max*0.95))
plt.annotate(r'$\frac{\nu_{max}}{2}$', xy = (0, nu_max/2), xytext = (Km/5, nu_max/2*1.05))
plt.xlabel(r'$\text{S (mmol.L}^{-1}\text{)}$')
plt.ylabel(r'$\text{1/r}_{s}\text{ (mmol.L}^{-1}\text{.s}^{-1}\text{)}$')
plt.xlim(0, S_model[-1])
plt.ylim(0, nu_max*1.1)
plt.grid(linewidth = 0.25)
plt.legend()
plt.title("Michaelis-Menten plot")
plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\michaelis_menten_plot.eps', format = 'eps')
plt.show()