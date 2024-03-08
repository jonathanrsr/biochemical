import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

time_cell = np.array([0.0, 0.54, 0.9, 1.23, 1.58, 1.95, 2.33, 2.7])
cell_concentration = np.array([15.5, 23.0, 30.0, 38.8, 48.5, 58.3, 61.3, 62.5])
time_lactose = np.array([0.0, 0.54, 0.54, 0.9, 0.9, 1.23, 1.23, 1.58, 1.58, 1.95, 1.95, 2.33, 2.33, 2.7])
lactose_concentration = np.array([137.0, 137.0, 114.0, 114.0, 90.0, 90.0, 43.0, 43.0, 29.0, 29.0, 9.0, 9.0, 2.0, 2.0])

plt.plot(time_cell, cell_concentration, label = 'Cell concentration (g/L)')
#plt.plot(time_lactose, lactose_concentration, label = 'Lactose concentration (g/L)')
plt.legend()
plt.grid(axis = 'both')
plt.show()

plt.plot(time_cell, np.log(cell_concentration))
plt.show()

mu = np.zeros(7)
for x in range(len(mu)):
    mu[x] = 1/cell_concentration[x]*(cell_concentration[x + 1] - cell_concentration[x])/(time_cell[x + 1] - time_cell[x])


lactose_plot = np.array([137, 114, 90, 43, 29, 9, 2])
plt.plot(1/lactose_plot, 1/mu)
plt.show()

mu_exp = mu[0:5]
lactose_plot = np.array([137, 114, 90, 43, 29])

slope = Polynomial.fit(1/lactose_plot, 1/mu_exp, 1).convert().coef
print(slope)

y = 1/lactose_plot*slope[1] + slope[0]

#Faire pente jusqua abscisse

plt.plot(1/lactose_plot, 1/mu_exp, 1/lactose_plot, y)
plt.ylim([0, 2])
plt.xlim([0, 0.035])
plt.show()

mu_max = 1/slope[0]
Ks = slope[1]*mu_max

print(mu_max, Ks)

dt = np.log(2)/mu