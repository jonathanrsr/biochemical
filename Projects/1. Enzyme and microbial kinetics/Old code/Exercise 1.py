import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

# We set the initial concentrations and rate constants for the reaction system.
S0 = 1 #M
E0 = 0.1 #M
P0 = ES0 = 0 #M

kf = 10e4 #M-1.s-1
kb = 20 #s-1
kcat = 10 #s-1

# We define the system of differential equations representing the rate of change of reactants and products in the system.
def equations(variables,t):
    E,S,ES,P = variables
    
    dE = (kb+kcat)*ES-kf*E*S
    dS = -kf*E*S+kb*ES
    dES = -(kb+kcat)*ES+kf*E*S
    dP = kcat*ES
    
    return dE, dS, dES, dP

# We specify the initial conditions for the system and create a time array from 0 to 5 seconds.
# with a step size of 0.0001 seconds.
initial_conditions = [E0, S0, ES0, P0]
t = np.arange(start=0, stop=5, step=0.0001)

# We use the odeint function to solve the system of ordinary differential equations defined by the equations function,
# with the specified initial_conditions and over the time interval defined by t.
results = odeint(equations, initial_conditions, t)

# We unpack the results into separate arrays for enzyme (E), substrate (S), enzyme-substrate complex (ES),
# and product (P).
E, S, ES, P = results.T

# We plot the species evolution over time.
plt.plot(t, results)
plt.xlabel("Time [s]")
plt.ylabel("Molar concentration [mol/L]")
plt.legend("E,S,ES,P".split(","))

plt.grid()
plt.title("Species' evolution over time", y=-0.19, fontweight='bold', fontsize=11)

plt.xlim(left=0, right=5) 
plt.ylim(bottom=0, top=1.01) 

plt.xticks(np.arange(0, 5.1, 0.5)) 
plt.yticks(np.arange(0, 1.1, 0.1))

plt.show()

# We define the system of differential equations, representing the rate of change of reactants and products in the
# system, obtained upon making the Quasi-Steady-State approximation.
def equations_approx(variables_a, t):
    E_a, S_a, ES_a, P_a = variables_a
    
    ES_a = (kf*S_a*E0)/(kb+kcat+kf*S_a)
    
    dE_a = 0
    dS_a = -kf*(E0 - ES_a)*S_a+kb*ES_a
    dES_a = 0
    dP_a = kcat * ES_a
    
    return dE_a, dS_a, dES_a, dP_a

# We specify the initial conditions for the system and create a time array from 0 to 5 seconds.
# with a step size of 0.0001 seconds.
ini_conditions = [E0, S0, ES0, P0]
t = np.arange(start=0, stop=5, step=0.0001)

# We use the odeint function to solve the system of ordinary differential equations defined by the equations function,
# with the specified initial_conditions and over the time interval defined by t.
results_a = odeint(equations_approx, ini_conditions, t)

# We unpack the results into separate arrays for enzyme (E), substrate (S), enzyme-substrate complex (ES),
# and product (P).
E_a, S_a, ES_a, P_a = results_a.T

# We plot the product's evolution over time, resulting from solving both the QSSA-system of ODEs and the exact one.
plt.plot(t, P_a, label='Product concentration (with approximation)')
plt.plot(t, P, label='Product concentration (no approximation)')

plt.xlabel("Time [s]")
plt.ylabel("Product concentration [mol/L]")
plt.title("Product concentration over time", y=-0.19, fontweight='bold', fontsize=11)
plt.legend()

plt.xlim(left=0, right=5) 
plt.ylim(bottom=0, top=2.01) 

plt.xticks(np.arange(0, 5.1, 0.5)) 
plt.yticks(np.arange(0, 2.01, 0.1))

plt.grid()
plt.show()