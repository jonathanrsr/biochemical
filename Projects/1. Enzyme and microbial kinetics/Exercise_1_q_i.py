import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
from sklearn.linear_model import LinearRegression

#which parameters would you change to increase selectivity?

def equation (variables, t):
    CE, CS, CES, CP, CP_1 = variables / V
    
    r1 = kf * CE * CS
    r2 = kb *CES
    r3 = kcat1 * CES
    r4 = kcat2 *CES
    
    RE = -r1 +r2 +r3 + r4
    RS = -r1 +r2
    RES= r1 -r2 -r3 - r4
    RP= r3
    RP_1= r4
    
    dNE_dt = RE*V
    dNS_dt = RS*V
    dNES_dt = RES * V
    dNP_dt = RP * V
    dNP_1_dt = RP_1*V
    
    return dNE_dt, dNS_dt, dNES_dt, dNP_dt, dNP_1_dt


kf = 10000 #l/mol/s

kb = 200 #1/s

kcat1 = 10 #1/s

kcat2 = 5 #1/s

E0 = 1 #mol/l

initial_conditions = [E0*V, S0*V, ES0*V, P0*V, P0_1*V]

results = odeint(equation, initial_conditions, t)
CE, CS, CES, CP, CP_1 = results.T/V

plt.figure(1, figsize=(8,7))
plt.title("Non Specific Enzymatic Reaction",  size=20)
plt.plot(t, CE, t, CS, t, CES, t, CP, t, CP_1)
plt.xlabel("time [s]", size=15)
plt.ylabel("Ci [mol/l]", size=15)
plt.legend(["E", "S", "ES", "P", "P'"])

plt.grid()
plt.savefig('non1_enzymatic_reaction.eps', format='eps')
plt.show()

print(CP_1[-1])
print(CP[-1])
