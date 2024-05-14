#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


##############################################################
############        code chunk 1             #################
########### Setting the parameters and IC    #################
##############################################################


"""
   Initial conditions
"""


N = 2200000
Hw0 = 5
Aw0 = 5
Hm0 = 0
Am0 = 0 
S0 = N - Hw0 -Aw0 - Hm0 - Am0
R0 = 0

"""
   model parameters
"""

gamma = 0.05  # recovery rate

beta_w = 1.5*gamma # transmission rate wild type
phi_w = 1.6 # change of infectivity Aw versus the infectivity of Hw

beta_m = 1.5*gamma # transmission rate mutant
phi_m = 2.0 # change of infectivity Am versus the infectivity of Hm

eta_w =0.45 # fraction that develops severe disease after wild-type infection 
eta_m =0.2 # fraction that develops severe disease after mutant infection

eps = 0.0001 # fraction of infections that mutate
alpha = 1/180 # 1/alpha = mean duration of natural immunity

rho = 1e-5 # spillover transmission rate 



#############################################################
############# Model equations ################################
##############################################################



def model(y, t):
    S, Hw, Aw, Hm, Am, R=y
    # Equations
    dSdt = -beta_w*(S/N)*(Hw + phi_w*Aw + rho*N) - beta_m*(S/N)*(Hm + phi_m*Am) + alpha*R
    dHwdt = eta_w*(1.-eps)*beta_w*(S/N)*(Hw + phi_w*Aw + rho*N) - gamma*Hw
    dAwdt = (1.-eta_w)*(1.-eps)*beta_w*(S/N)*(Hw + phi_w*Aw + rho*N) - gamma*Aw
    dHmdt = eta_m*beta_m*(S/N)*(Hm + phi_m*Am) + eta_w*eps*beta_w*(S/N)*(Hw + phi_w*Aw + rho*N) - gamma*Hm
    dAmdt = (1.-eta_m)*beta_m*(S/N)*(Hm + phi_m*Am) + (1.-eta_w)*eps*beta_w*(S/N)*(Hw + phi_w*Aw + rho*N) - gamma*Am
    dRdt = gamma*(Hw + Aw + Hm + Am) - alpha*R
    return dSdt, dHwdt, dAwdt, dHmdt, dAmdt, dRdt


################################################################
############# Integration and plots  ###########################
################################################################

# A grid of time points (in days)
t = np.linspace(0, 1000, 1001)

# Initial conditions vector
y0 = S0, Hw0, Aw0, Hm0, Am0, R0

# Solution of the model
sol = odeint(model, y0, t)
S, Hw, Aw, Hm, Am, R = sol.T

plt.style.use('seaborn-paper')

############# infected class  ##################################
fig, ax = plt.subplots(figsize=(11,7), tight_layout=True)

ax.plot(t, Hw, linewidth=3, label=r'$H_w$')
ax.plot(t, Aw, linewidth=3, label=r'$A_w$')
ax.plot(t, Hm, linewidth=3, label=r'$H_m$')
ax.plot(t, Am, linewidth=3, label=r'$A_m$')
#ax.plot(t, S, linewidth=3, label=r'$S$')

# Customise some display properties
ax.set_xlabel('Time (days)', fontsize=20)
ax.set_ylabel('Hospitalized individuals', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
#ax.set_title('Hospitalized class', labelsize=30)
ax.legend(loc='best', fontsize=20)
#plt.savefig("endemic_reservoir.pdf", bbox_inches = 'tight')
plt.show()



