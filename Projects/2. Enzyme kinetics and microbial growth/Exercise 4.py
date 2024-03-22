import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

######## - DATA - ########

# Dilution rate, 1/h
D = np.array([0.05, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.84])

# Cell concentration, g/L
X = np.array([3.2, 3.7, 4, 4.4, 4.75, 4.9, 4.5, 0.5])

# Carbon substrate concentration, g/L
S = np.array([0.012, 0.028, 0.05, 0.1, 0.15, 0.176, 0.8, 9])

# Carbon substrate concentration in the feed, g/L
S_feed = 10



######## - Estimate (assuming Monod kinetics) maximum growth rate 
# (mu_max) and saturation/Monod constant (Ks) - ########

# Specific growth rate mu is equal to dilutation rate D = mu

# Lineweaver-Burke plot can be used to determine mumax and Ks by doing
# a linear regression with y = 1/D and x = 1/S

x = (1/S).reshape((-1, 1))
y = 1/D
model = LinearRegression().fit(x, y)
y_intercept = float(model.intercept_) # ==> 1/mumax
slope = model.coef_[0] # ==> Ks/mumax
x_intercept = -y_intercept/slope # ==> -1/Ks
det_coefficient = model.score(x, y) #R^2

mumax = np.round(1/y_intercept, 2) # 1/h
Ks = np.round(-1/x_intercept, 2) # g/L

print(f'Question 1\nmumax = {mumax} 1/h\nKs = {Ks} g/L\n')

plt.scatter(x, y, color = 'black', s = 10.0, label = 'data')
plt.plot(x, model.predict(x), color = 'black', linewidth = 1.0,
         label = f'model, R^')
plt.show()