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


v0_kf_slope=v_0_kf_vary[np.log(kf_range)<=4]
kf_range_slope=kf_range[np.log(kf_range)<=4]
v0_kb_slope=v_0_kb_vary[np.log(kb_range)>=7.5]
kb_range_slope=kb_range[np.log(kb_range)>=7.5]
v0_kcat_slope=v_0_kcat_vary[np.log(kcat_range)>=7]
kcat_range_slope=kcat_range[np.log(kcat_range)>=7]
v0_E0_slope=v_0_E0_vary[E0_range>0]
E0_range_slope=E0_range[E0_range>0]

######### Fits #######
z_kf = np.polyfit(np.log(kf_range_slope), np.log(v0_kf_slope), 1)
x1_kf=np.linspace(2,4,100)
y1_kf = np.polyval(z_kf,x1_kf);

######### Fits ############

z_kb = np.polyfit(np.log(kb_range_slope), np.log(v0_kb_slope), 1)
x1_kb=np.linspace(7.5,9.5,100)
y1_kb = np.polyval(z_kb,x1_kb);

########## Fits #########

z_kcat = np.polyfit(np.log(kcat_range_slope), np.log(v0_kcat_slope), 1)
x1_kcat=np.linspace(7,9.2,100)
y1_kcat = np.polyval(z_kcat,x1_kcat);


########## Fits ###########

z_E0 = np.polyfit(np.log(E0_range_slope), np.log(v0_E0_slope), 1)
x1_E0=np.linspace(-4,0,100)
y1_E0 = np.polyval(z_E0,x1_E0);



############################
fig, axs = plt.subplots(2, 2, figsize=(12,13))
axs[0,0].set_title('(a) log(v(t=0)) vs log(kf)')
axs[0,0].scatter(np.log(kf_range_slope), np.log(v0_kf_slope), color = '#02666e', label = "data")
axs[0,0].plot(x1_kf, y1_kf, 'b', linewidth=2, label=f'fit y={z_kf[0]:.2}*x+{z_kf[1]:.2}')
axs[0,0].set_xlabel('log(kf) [log(l/mol/s)]', size=10)
axs[0,0].set_ylabel('log(Initial reaction rate), log(v(t=0)), [log(mol/l/s)]',  size=10)
axs[0,0].legend()
axs[0,0].grid(True, which="both")


axs[0,1].set_title('(b) log(v(t=0)) vs log(kb)')
axs[0,1].scatter(np.log(kb_range_slope), np.log(v0_kb_slope), color = '#9e0505', label = "data")
axs[0,1].plot(x1_kb, y1_kb, 'r', linewidth=2, label=f'fit y={z_kb[0]:.2}*x+{z_kb[1]:.2}')
axs[0,1].set_xlabel('log(kb) [log(1/s)]', size=10)
axs[0,1].set_ylabel('log(Initial reaction rate), log(v(t=0)), [log(mol/l/s)]',  size=10)
axs[0,1].legend()
axs[0,1].grid(True, which="both")

axs[1,0].set_title('(c) log(v(t=0)) vs log(kcat)')
axs[1,0].scatter(np.log(kcat_range_slope), np.log(v0_kcat_slope), color = '#165e28', label = "data")
axs[1,0].plot(x1_kcat, y1_kcat, 'g', linewidth=2, label=f'fit y={z_kcat[0]:.2}*x+{z_kcat[1]:.2}')
axs[1,0].set_xlabel('log(kcat) [log(1/s)]', size=10)
axs[1,0].set_ylabel('log(Initial reaction rate), log(v(t=0)), [log(mol/l/s)]',  size=10)
axs[1,0].legend()
axs[1,0].grid(True, which="both")

axs[1,1].set_title('(d) log(v(t=0)) vs log(E0)')
axs[1,1].scatter(np.log(E0_range_slope), np.log(v0_E0_slope), color = '#f3ba3f', label="data")
axs[1,1].plot(x1_E0, y1_E0, 'y', linewidth=2, label=f'fit y={z_E0[0]:.2}*x+{z_E0[1]:.2}')
axs[1,1].set_xlabel('log(E0) [log(mol/l)]', size=10)
axs[1,1].set_ylabel('log(Initial reaction rate), log(v(t=0)), [log(mol/l/s)]',  size=10)
axs[1,1].legend()
axs[1,1].grid(True, which="both")

plt.savefig('logloginitialvelocity.eps', format='eps')
