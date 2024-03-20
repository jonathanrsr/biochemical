import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression

# We set the initial concentrations and rate constants for the reaction system.
S0 = 1 #M
E0 = 0.1 #M
P0 = ES0 = 0 #M

kf = 10e4 #M-1.s-1
kb = 20 #s-1
kcat = 10 #s-1

# We define the system of differential equations representing the rate of change of reactants and products in the system.
def equations(variables,t):
    E, S, ES, P = variables
    
    dE = (kb + kcat)*ES - kf*E*S
    dS = -kf*E*S + kb*ES
    dES = -(kb + kcat)*ES + kf*E*S
    dP = kcat*ES
    
    return dE, dS, dES, dP

# We specify the initial conditions for the system and create a time array from 0 to 5 seconds.
# with a step size of 0.0001 seconds.
initial_conditions = [E0, S0, ES0, P0]
t = np.arange(start = 0, stop = 5, step = 0.0001)

# We use the odeint function to solve the system of ordinary differential equations defined by the equations function,
# with the specified initial_conditions and over the time interval defined by t.
results = odeint(equations, initial_conditions, t)

# We unpack the results into separate arrays for enzyme (E), substrate (S), enzyme-substrate complex (ES),
# and product (P).
E, S, ES, P = results.T

# We plot the species evolution over time.
plt.plot(t, results, linewidth = 1.0)
plt.xlabel("Time (s)")
plt.ylabel(r"$\text{Species' concentration (mol.L}^{-1}\text{)}$")
plt.legend("E,S,ES,P".split(","))

plt.grid(linewidth = 0.25)
plt.title("Species' evolution over time")

plt.xlim(left = 0, right = 5) 
plt.ylim(bottom = 0, top = 1.1) 

plt.xticks(np.arange(0, 5.1, 0.5)) 
plt.yticks(np.arange(0, 1.2, 0.1))

plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\fig1.eps', format='eps')
plt.show()

# We define the system of differential equations, representing the rate of change of reactants and products in the
# system, obtained upon making the Quasi-Steady-State approximation.
def equations_approx(variables_a, t):
    E_a, S_a, ES_a, P_a = variables_a
    
    ES_a = (kf*S_a*E0)/(kb +kcat +kf*S_a)
    
    dE_a = 0
    dS_a = -kf*(E0 - ES_a)*S_a + kb*ES_a
    dES_a = 0
    dP_a = kcat*ES_a
    
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
plt.plot(t, P_a, label='Product concentration (with approximation)', linewidth = 1.0, color = 'blue')
plt.plot(t, P, label='Product concentration (no approximation)', linewidth = 1.0, color = 'red')

plt.xlabel("Time (s)")
plt.ylabel(r"$\text{Product concentration (mol.L}^{-1}\text{)}$")
plt.title("Product concentration over time")
plt.legend()

plt.xlim(left=0, right=5) 
plt.ylim(bottom=0, top=1.1) 

plt.xticks(np.arange(0, 5.1, 0.5)) 
plt.yticks(np.arange(0, 1.2, 0.1))

plt.grid(linewidth = 0.25)
plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\fig2.eps', format='eps')
plt.show()



#Exercise 1, Q) e
S0 = 1 #mol/l
E0= 0.1 #mol/l
P0 = 0 #mol/l
ES0=P0 #mol/l
kf = 10**4 #l/(mols)
kb = 20 #1/s
kcat = 10 #1/s

# Creating ranges for the parameters one wants to modify
kf_range = np.logspace(1, 8) #L/mol
kb_range = np.logspace(1, 4) #1/s
kcat_range = np.logspace(1, 4) #1/s
E0_range = np.linspace(0,1) #mol/L

v_0_kf_vary = kcat*E0*S0/(((kb + kcat)/kf_range) + S0)
v_0_kb_vary = kcat*E0*S0/(((kb_range + kcat)/kf) + S0)
v_0_kcat_vary = kcat_range*E0*S0/(((kb + kcat_range)/kf) +S0)
v_0_E0_vary = kcat*E0_range*S0/(((kb + kcat)/kf) + S0)

#Plots 
fig, axs = plt.subplots(2, 2, figsize=(12,13))
axs[0,0].set_title(r'$\text{(a) }\text{Initial reaction rate vs k}_{f}$')
axs[0,0].semilogx(kf_range, v_0_kf_vary, color = 'black', linewidth = 1.0)
axs[0,0].set_xlabel(r'$\text{k}_{f}\text{ (L.mol}^{-1}\text{.s}^{-1}\text{)}$')
axs[0,0].set_ylabel(r'$\nu_{init}\text{ (mol.L}^{-1}\text{.s}^{-1}\text{)}$')
axs[0,0].grid(True, which="both", linewidth = 0.25)

axs[0,1].set_title(r'$\text{(b) }\text{Initial reaction rate vs k}_{b}$')
axs[0,1].semilogx(kb_range, v_0_kb_vary, color = 'black', linewidth = 1.0)
axs[0,1].set_xlabel(r'$\text{k}_{b}\text{ (s}^{-1}\text{)}$')
axs[0,1].set_ylabel(r'$\nu_{init}\text{ (mol.L}^{-1}\text{.s}^{-1}\text{)}$')
axs[0,1].grid(True, which="both", linewidth = 0.25)

axs[1,0].set_title(r'$\text{(c) }\text{Initial reaction rate vs k}_{cat}$')
axs[1,0].semilogx(kcat_range, v_0_kcat_vary, color = 'black', linewidth = 1.0)
axs[1,0].set_xlabel(r'$\text{k}_{cat}\text{ (s}^{-1}\text{)}$')
axs[1,0].set_ylabel(r'$\nu_{init}\text{ (mol.L}^{-1}\text{.s}^{-1}\text{)}$')
axs[1,0].grid(True, which="both", linewidth = 0.25)

axs[1,1].set_title('(c) Initial reaction rate vs initial enzyme concentration')
axs[1,1].plot(E0_range, v_0_E0_vary, color = 'black', linewidth = 1.0)
axs[1,1].set_xlabel(r'$\text{E}_{0}\text{ (mol.L}^{-1}\text{)}$')
axs[1,1].set_ylabel(r'$\nu_{init}\text{ (mol.L}^{-1}\text{.s}^{-1}\text{)}$')
axs[1,1].grid(True, which="both", linewidth = 0.25)

plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\Initial_velocity.eps', format='eps')
plt.show()


#Exercise 1, Q) f
v0_kf_slope = v_0_kf_vary[np.log(kf_range) <= 4]
kf_range_slope = kf_range[np.log(kf_range) <= 4]
v0_kb_slope = v_0_kb_vary[np.log(kb_range) >= 7.5]
kb_range_slope = kb_range[np.log(kb_range) >= 7.5]
v0_kcat_slope = v_0_kcat_vary[np.log(kcat_range) >= 7]
kcat_range_slope = kcat_range[np.log(kcat_range) >= 7]
v0_E0_slope = v_0_E0_vary[E0_range > 0]
E0_range_slope = E0_range[E0_range > 0]

######### Fits #########
z_kf = np.polyfit(np.log(kf_range_slope), np.log(v0_kf_slope), 1)
x1_kf=np.linspace(2, 4, 100)
y1_kf = np.polyval(z_kf,x1_kf);
R_kf = np.corrcoef(np.log(kf_range_slope), np.log(v0_kf_slope))[0, 1]

######### Fits #########
z_kb = np.polyfit(np.log(kb_range_slope), np.log(v0_kb_slope), 1)
x1_kb=np.linspace(7.5, 9.5, 100)
y1_kb = np.polyval(z_kb,x1_kb);
R_kb = np.corrcoef(np.log(kb_range_slope), np.log(v0_kb_slope))[0, 1]

######### Fits #########
z_kcat = np.polyfit(np.log(kcat_range_slope), np.log(v0_kcat_slope), 1)
x1_kcat=np.linspace(7, 9.2, 100)
y1_kcat = np.polyval(z_kcat,x1_kcat);
R_kcat = np.corrcoef(np.log(kcat_range_slope), np.log(v0_kcat_slope))[0, 1]

######### Fits #########
z_E0 = np.polyfit(np.log(E0_range_slope), np.log(v0_E0_slope), 1)
x1_E0=np.linspace(-4, 0, 100)
y1_E0 = np.polyval(z_E0,x1_E0);
R_E0 = np.corrcoef(np.log(E0_range_slope), np.log(v0_E0_slope))[0, 1]

######### Plot #########
fig, axs = plt.subplots(2, 2, figsize=(12, 13))
axs[0,0].set_title(r'$\text{Initial reaction rate (ln) vs k}_{f}\text{ (ln)}$')
axs[0,0].scatter(np.log(kf_range_slope), np.log(v0_kf_slope), color = 'black', s = 10.0, label = "data")
axs[0,0].plot(x1_kf, y1_kf, color = 'black', linestyle = '--', linewidth = 1.0, label = r'$\text{model, R}^{2} = ' + r'\text{' + f'{round(R_kf**2, 3):.3f}' + r'}$')
axs[0,0].set_xlabel(r'$\text{ln(k}_{f}\text{)}$')
axs[0,0].set_ylabel(r'$\text{ln(}\nu_{init}\text{)}$')
axs[0,0].legend()
axs[0,0].grid(True, which="both", linewidth = 0.25)

axs[0,1].set_title(r'$\text{Initial reaction rate (ln) vs k}_{b}\text{ (ln)}$')
axs[0,1].scatter(np.log(kb_range_slope), np.log(v0_kb_slope), color = 'black', s = 10.0, label = "data")
axs[0,1].plot(x1_kb, y1_kb, color = 'black', linestyle = '--', linewidth = 1.0, label = r'$\text{model, R}^{2} = ' + r'\text{' + f'{round(R_kb**2, 3):.3f}' + r'}$')
axs[0,1].set_xlabel(r'$\text{ln(k}_{b}\text{)}$')
axs[0,1].set_ylabel(r'$\text{ln(}\nu_{init}\text{)}$')
axs[0,1].legend()
axs[0,1].grid(True, which="both", linewidth = 0.25)

axs[1,0].set_title(r'$\text{Initial reaction rate (ln) vs k}_{cat}\text{ (ln)}$')
axs[1,0].scatter(np.log(kcat_range_slope), np.log(v0_kcat_slope), color = 'black', s = 10.0, label = "data")
axs[1,0].plot(x1_kcat, y1_kcat, color = 'black', linestyle = '--', linewidth = 1.0, label = r'$\text{model, R}^{2} = ' + r'\text{' + f'{round(R_kcat**2, 3):.3f}' + r'}$')
axs[1,0].set_xlabel(r'$\text{ln(k}_{cat}\text{)}$')
axs[1,0].set_ylabel(r'$\text{ln(}\nu_{init}\text{)}$')
axs[1,0].legend()
axs[1,0].grid(True, which="both", linewidth = 0.25)

axs[1,1].set_title(r'$\text{Initial reaction rate (ln) vs initial enzyme concentration (ln)}$')
axs[1,1].scatter(np.log(E0_range_slope), np.log(v0_E0_slope), color = 'black', s = 10.0, label = "data")
axs[1,1].plot(x1_E0, y1_E0, color = 'black', linestyle = '--', linewidth = 1.0, label = r'$\text{model, R}^{2} = ' + r'\text{' + f'{round(R_E0**2, 3):.3f}' + r'}$')
axs[1,1].set_xlabel(r'$\text{ln(E}_{0}\text{)}$')
axs[1,1].set_ylabel(r'$\text{ln(}\nu_{init}\text{)}$')
axs[1,1].legend()
axs[1,1].grid(True, which="both", linewidth = 0.25)

plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\logloginitialvelocity.eps', format='eps')
plt.show()

#Exercise 1, Q) h
t = np.linspace(start = 0, stop= 5, num = 5001) #s

S0 = 1 #mol/l
E0 = 0.1 #mol/l
P0 = 0 #mol/l
ES0 = P0 #mol/l
P0_1 = P0 #mol/l
kf = 10000 #l/mol/s
kb = 20 #1/s
kcat1 = 10 #1/s
kcat2 = 5 #1/s
V = 100 #l
initial_conditions = [E0*V, S0*V, ES0*V, P0*V, P0_1*V]

#Function that allows to solve the mass balance 
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

results = odeint(equation, initial_conditions, t)
CE, CS, CES, CP, CP_1 = results.T/V

#Plot of the system
plt.plot(t, CE, t, CS, t, CES, t, CP, t, CP_1, linewidth = 1.0)
plt.title('Non Specific Enzymatic Reaction')
plt.xlabel('time (s)')
plt.ylabel(r'$\text{C}_{i}\text{ (mol.L}^{-1}\text{)}$')
plt.legend(['E', 'S', 'ES', 'P', 'P'])
plt.grid(linewidth = 0.25)
plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\non_enzymatic_reaction.eps', format='eps')
plt.show()

print(CP_1[-1]) #Undesired product concentration at steady state
print(CP[-1]) #Desired productconcentration at steady state
print(t[CP>=0.66][0]) # Time to steady state

#Exercise 1, Q) i
#which parameters would you change to increase selectivity?
#increase kcat1 will allow better selectivity. Example: 
kf = 10000 #l/mol/s
kb = 200 #1/s
kcat1 = 100 #1/s
kcat2 = 5 #1/s
E0 = 1 #mol/l

initial_conditions = [E0*V, S0*V, ES0*V, P0*V, P0_1*V]
results = odeint(equation, initial_conditions, t)
CE, CS, CES, CP, CP_1 = results.T/V

#Plot figure of example change in selectivity
plt.plot(t, CE, t, CS, t, CES, t, CP, t, CP_1, linewidth = 1.0)
plt.title("Non Specific Enzymatic Reaction")
plt.xlabel("time (s)")
plt.ylabel(r'$\text{C}_{i}\text{ (mol.L}^{-1}\text{)}$')
plt.xlim([0, 0.1])
plt.legend(["E", "S", "ES", "P", "P'"])

plt.grid(linewidth = 0.25)
plt.savefig('Projects\\1. Enzyme and microbial kinetics\\Images\\non1_enzymatic_reaction.eps', format='eps')
plt.show()

print(CP_1[-1]) #Undesired product concentration at steady state
print(CP[-1]) #Desired productconcentration at steady state