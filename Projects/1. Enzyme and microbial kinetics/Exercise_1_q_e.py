import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
from sklearn.linear_model import LinearRegression


S0 = 1 #mol/l

E0= 0.1 #mol/l

P0 = 0 #mol/l

ES0=P0 #mol/l

kf = 10**4 #l/(mols)

kb = 20 #1/s

kcat = 10 #1/s

kf_range = np.logspace(1, 8) #l/(mols)
kb_range = np.logspace(1, 4) #1/s
kcat_range = np.logspace(1, 4) #1/s
E0_range = np.linspace(0,1) #mol/l

v_0_kf_vary = kcat*E0*S0/(((kb+kcat)/kf_range) +S0)
v_0_kb_vary = kcat*E0*S0/(((kb_range+kcat)/kf) +S0)
v_0_kcat_vary = kcat_range*E0*S0/(((kb+kcat_range)/kf) +S0)
v_0_E0_vary = kcat*E0_range*S0/(((kb+kcat)/kf) +S0)

#plt.semilogx( kb_range, v_0_kb_vary)
fig, axs = plt.subplots(2, 2, figsize=(12,13))
axs[0,0].set_title('(a) v(t=0) vs kf')
axs[0,0].semilogx(kf_range,v_0_kf_vary,'b', linewidth=2)
axs[0,0].set_xlabel('kf [l/mol/s]', size=10)
axs[0,0].set_ylabel('Initial reaction rate, v(t=0), [mol/l/s]',  size=10)
axs[0,0].legend(["kf variation"], loc='lower right', frameon=False)
axs[0,0].grid(True, which="both")

axs[0,1].set_title('(b) v(t=0) vs kb')
axs[0,1].semilogx(kb_range,v_0_kb_vary,'r', linewidth=2)
axs[0,1].set_xlabel('kb [1/s]', size=10)
axs[0,1].set_ylabel('Initial reaction rate, v(t=0), [mol/l/s]',  size=10)
axs[0,1].legend(["kb variation"], loc='lower left', frameon=False)
axs[0,1].grid(True, which="both")

axs[1,0].set_title('(c) v(t=0) vs kcat')
axs[1,0].semilogx(kcat_range,v_0_kcat_vary,'g', linewidth=2)
axs[1,0].set_xlabel('kcat [1/s]', size=10)
axs[1,0].set_ylabel('Initial reaction rate, v(t=0), [mol/l/s]',  size=10)
axs[1,0].legend(["kcat variation"], loc='upper left', frameon=False)
axs[1,0].grid(True, which="both")

axs[1,1].set_title('(d) v(t=0) vs E0')
axs[1,1].plot(E0_range,v_0_E0_vary,'y', linewidth=2)
axs[1,1].set_xlabel('E0 [mol/l]', size=10)
axs[1,1].set_ylabel('Initial reaction rate, v(t=0), [mol/l/s]',  size=10)
axs[1,1].legend(["E0 variation"], loc='upper left', frameon=False)
axs[1,1].grid(True, which="both")

plt.savefig('Initial_velocity.eps', format='eps')
