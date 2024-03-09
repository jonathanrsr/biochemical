import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker
from sklearn.linear_model import LinearRegression

# Empirical data
sampling_time = np.array([0, 0.54, 0.9, 1.23, 1.58, 1.95, 2.33, 2.7]) # h
cell_concentration = np.array([15.5, 23.0, 30.0, 38.8, 48.5, 58.3, 61.3, 62.5]) # g/L
lactose_concentration = np.array([137.0, 114.0, 90.0, 43.0, 29.0, 9.0, 2.0]) # g/L

# Plot cell and lactose concentration as a function of time
fig, ax1 = plt.subplots()
ax1.stairs(lactose_concentration, sampling_time, baseline = None, color = 'red', linestyle = '--', linewidth = 1.0, label = 'Lactose concentration')
#ax1.plot(sampling_time[:-1] + np.diff(sampling_time)/2, lactose_concentration, color = 'r', marker = 'o', linestyle = 'None')
ax1.set_ylabel('Lactose concentration, $\mathregular{g.L^{-1}}$')
ax1.set_yticks(np.linspace(0, 160, 9))
ax1.legend(loc = 'upper left')
ax1.grid(linewidth = 0.25)

ax2 = ax1.twinx()
ax2.plot(sampling_time, cell_concentration, color = 'blue', marker = 'o', markersize = '4', linestyle = '--', linewidth = 1.0, label = 'Cell concentration')
ax2.set_xlabel('Time, h')
ax2.set_ylabel('Cell concentration, $\mathregular{g.L^{-1}}$')
ax2.grid(linewidth = 0.25)
ax2.set_yticks(np.linspace(0, 80, 9))
ax2.legend(loc = 'upper right')

plt.title('Cell and lactose concentration versus time')
plt.show()

# Plot semi-log graph of cell concentrationas a function of time
plt.semilogy(sampling_time, cell_concentration, color = 'black', marker = 'o', markersize = '4', linestyle = '--', linewidth = 1.0)
plt.xlabel('Time, h')
plt.ylabel('Cell concentration, ln')
plt.grid(which = 'both', linewidth = 0.25)
plt.gca().yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
plt.title('Cell concentration versus time, semi-log')
plt.show()

# Calculate growth rate mu = 1/Xmin*(dX/dt) using Euler's formula for the derivative
mu = np.zeros(len(cell_concentration) - 1)

for x in range(len(mu)):
    mu[x] = 1/cell_concentration[x]*(cell_concentration[x + 1] - cell_concentration[x])/(sampling_time[x + 1] - sampling_time[x])



# Lineweaver-Burke type plot
plt.plot(1/lactose_concentration, 1/mu, color = 'black', marker = 'o', markersize = '4', linestyle = 'None', label = 'Cell concentration')
plt.title('Lineweaver-Burke type plot')
plt.show()