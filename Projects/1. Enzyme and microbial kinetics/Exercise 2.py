import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# Empirical data
sampling_time = np.array([0, 0.54, 0.9, 1.23, 1.58, 1.95, 2.33, 2.7]) # h
cell_concentration = np.array([15.5, 23.0, 30.0, 38.8, 48.5, 58.3, 61.3, 62.5]) # g/L
lactose_concentration = np.array([137.0, 114.0, 90.0, 43.0, 29.0, 9.0, 2.0]) # g/L


# Plot cell and lactose concentration as a function of time
fig, ax1 = plt.subplots()
ax1.stairs(lactose_concentration, sampling_time, baseline = None, color = 'red', linestyle = '--', linewidth = 1.0, label = 'Lactose concentration')
ax1.set_xlabel('Time, h')
ax1.set_ylabel(r'$\text{Lactose concentration, g.L}^{-1}$')
ax1.set_yticks(np.linspace(0, 160, 9))
ax1.grid(linewidth = 0.25)
ax1.legend(loc = 'upper left')

ax2 = ax1.twinx()
ax2.plot(sampling_time, cell_concentration, color = 'blue', marker = 'o', markersize = '4', linestyle = '--', linewidth = 1.0, label = 'Cell concentration')
ax2.set_ylabel(r'$\text{Cell concentration, g.L}^{-1}$')
ax2.set_yticks(np.linspace(0, 80, 9))
ax2.grid(linewidth = 0.25)
ax2.legend(loc = 'upper right')

plt.title('Cell and lactose concentration versus time')
plt.show()


# Calculate growth rate mu = 1/Xmin*(dX/dt) using Euler's formula for the derivative
mu = np.zeros(len(lactose_concentration))
for x in range(len(mu)):
    mu[x] = 1/cell_concentration[x]*(cell_concentration[x + 1] - cell_concentration[x])/(sampling_time[x + 1] - sampling_time[x])


# Plot cell concentrations and mu as a function of time
fig, ax1 = plt.subplots()
ax1.stairs(mu, sampling_time, baseline = None, color = 'red', linestyle = '--', linewidth = 1.0, label = r'Growth rate')
ax1.set_ylabel(r'$\text{Growth rate, g.L}^{-1}\text{.s}^{-1}$')
ax1.set_yticks(np.linspace(0, 1, 11))
ax1.grid(linewidth = 0.25)
ax1.legend(loc = 'upper left')

ax2 = ax1.twinx()
ax2.plot(sampling_time, np.log(cell_concentration) + np.log(cell_concentration[0]), color = 'blue', marker = 'o', markersize = '4', linestyle = '--', linewidth = 1.0, label = 'Cell concentration')
ax2.set_xlabel('Time, h')
ax2.set_ylabel('Cell concentration, ln')
ax2.set_yticks(np.linspace(5.4, 7.1, 11))
ax2.set_ylim(5.4, 7.1)
ax2.legend(loc = 'upper right')

plt.title('Cell concentration and growth rate versus time')
plt.show()


# Slicing data to only take the growth phase part
index_growth_phase_start = 0
index_growth_phase_finish = 6
cell_concentration_growth_phase = cell_concentration[index_growth_phase_start:index_growth_phase_finish]
lactose_concentration_growth_phase = lactose_concentration[index_growth_phase_start:index_growth_phase_finish - 1]
mu = mu[index_growth_phase_start:index_growth_phase_finish - 1]


# Linear regression of 1/mu versus 1/lactose concentration
x = (1/lactose_concentration_growth_phase).reshape((-1, 1))
y = 1/mu
model = LinearRegression().fit(x, y)
y_intercept = float(model.intercept_)
slope = model.coef_[0]
x_intercept = -y_intercept/slope
det_coefficient = model.score(x, y)

print(f'Slope: {round(slope, 2)}\nIntercept (x): {round(x_intercept, 2)}\nIntercept (y): {round(y_intercept, 2)}')


# Find kinetics parameter from the linear regression
Ks = -1/x_intercept
mu_max = 1/y_intercept


# Create 2 points for plotting the model
x_graph = np.array([x_intercept, 1/lactose_concentration_growth_phase[-1]*1.2]) # Multiply by 1.2 to make the line go further than the last data point
y_graph = model.predict(x_graph.reshape((-1, 1)))


# Lineweaver-Burke type plot
plt.plot(1/lactose_concentration_growth_phase, 1/mu, color = 'black', marker = 'o', markersize = '4', linestyle = 'None', label = 'data')
plt.plot(x_graph, y_graph, color='black', marker='None', linestyle='--', linewidth=1.0, label = r'$\text{model, R}{^2} = ' + r'\text{' + str(round(det_coefficient, 2)) + r'}$')
plt.annotate(r'$\frac{1}{\mu_{\mathrm{max}}}$', xy = (0, y_intercept), xycoords = 'data', xytext = (-50, 25), textcoords = 'offset points', arrowprops = dict(arrowstyle = '-|>', color = 'black'))
plt.annotate(r'-$\frac{1}{K_{\mathrm{s}}}$', xy = (x_intercept, 0), xycoords = 'data', xytext = (-6, 50), textcoords = 'offset points', arrowprops = dict(arrowstyle = '-|>', color = 'black'))
plt.xlabel('1/S')
plt.ylabel('1/μ')
plt.vlines(0, 0, np.max(y_graph)*1.05, color = 'black', linestyle = '-', linewidth = 1.0) # Multiply by 1.05 for better design
plt.ylim([0, np.max(y_graph)*1.05]) # Multiply by 1.05 for better design
plt.grid(which = 'both', linewidth = 0.25)
plt.legend(loc = 'upper left')
plt.title('Lineweaver-Burke type plot')
plt.show()


# Find doubling-time, using eq. X = exp(mu*t)*X0 and X = 2*X0 ==> t = ln(2)/u
mu_average = np.average(mu[:3]) # Only the first 3 growth rates are used, before its decrease
doubling_time = np.log(2)/mu_average
max_doubling_time = np.log(2)/mu_max


# Print kinetics parameters
print(f'Ks = {round(Ks, 2)}\nμmax = {round(mu_max, 2)}\nμaverage = {round(mu_average, 2)}\ndoubling time = {round(doubling_time, 2)} h\nMax doubling time = {round(max_doubling_time, 2)}')