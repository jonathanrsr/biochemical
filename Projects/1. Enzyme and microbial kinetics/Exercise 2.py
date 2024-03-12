import numpy as np
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
ax1.set_ylabel('Lactose concentration, g.L⁻¹')
ax1.set_yticks(np.linspace(0, 160, 9))
ax1.grid(linewidth = 0.25)
ax1.legend(loc = 'upper left')

ax2 = ax1.twinx()
ax2.plot(sampling_time, cell_concentration, color = 'blue', marker = 'o', markersize = '4', linestyle = '--', linewidth = 1.0, label = 'Cell concentration')
ax2.set_xlabel('Time, h')
ax2.set_ylabel('Cell concentration, g.L⁻¹')
ax2.set_yticks(np.linspace(0, 80, 9))
ax2.grid(linewidth = 0.25)
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

# Slicing data to only take the growth phase part
index_growth_phase_start = 0
index_growth_phase_finish = 5
cell_concentration_growth_phase = cell_concentration[index_growth_phase_start:index_growth_phase_finish]
lactose_concentration_growth_phase = lactose_concentration[index_growth_phase_start:index_growth_phase_finish - 1]

# Calculate growth rate mu = 1/Xmin*(dX/dt) using Euler's formula for the derivative
mu = np.zeros(len(lactose_concentration_growth_phase))

for x in range(len(mu)):
    mu[x] = 1/cell_concentration[x]*(cell_concentration[x + 1] - cell_concentration[x])/(sampling_time[x + 1] - sampling_time[x])

# Linear regression of 1/mu versus 1/lactose concentration
x = (1/lactose_concentration_growth_phase).reshape((-1, 1))
y = 1/mu
model = LinearRegression().fit(x, y)
y_intercept = float(model.intercept_)
slope = model.coef_[0]
x_intercept = -y_intercept/slope
det_coefficient = model.score(x, y)

# Find kinetics parameter from the linear regression
Ks = -1/x_intercept
mu_max = 1/y_intercept

# Create 2 points for plotting the model
x_graph = np.array([x_intercept, 1/lactose_concentration_growth_phase[-1]*1.2])
y_graph = model.predict(x_graph.reshape((-1, 1)))

# Lineweaver-Burke type plot
plt.plot(1/lactose_concentration_growth_phase, 1/mu, color = 'black', marker = 'o', markersize = '4', linestyle = 'None', label = 'Cell concentration⁻¹')
plt.plot(x_graph, y_graph, color='black', marker='None', linestyle='--', linewidth=1.0, label=f'Cell concentration⁻¹, predicted (R² = {round(det_coefficient, 2)})')
plt.xlabel('1/X')
plt.ylabel('1/μ')
plt.vlines(0, 0, np.max(y_graph)*1.05, color = 'black', linestyle = '-', linewidth = 1.0)
plt.xlim([-0.07, 0.04])
plt.ylim([0, np.max(y_graph)*1.05])
plt.grid(which = 'both', linewidth = 0.25)
plt.legend(loc = 'upper left')
plt.title('Lineweaver-Burke type plot')
plt.show()

# Find doubling-time, using eq. X = exp(mu*t)*X0 and X = 2*X0 ==> t = ln(2)/u
mu_average = np.average(mu)
doubling_time = np.log(2)/mu_average

# Print kinetics parameters
print(f'Ks = {round(Ks, 2)}\nμₘₐₓ = {round(mu_max, 2)}\ndoubling time = {round(doubling_time, 2)} h')