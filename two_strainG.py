#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator



##############################################################
############        code chunk 1             #################
########### Setting the parameters and IC    #################
##############################################################


"""
   Initial conditions
"""


N = 100000
Hw0 = 1000 # IC 10 # EE 52000
Aw0 = 2000 # IC 10 # EE 63000
Hm0 = 30
Am0 = 50 
S0 = N - Hw0 -Aw0 - Hm0 - Am0
R0 = 0

"""
   model parameters
"""

recoveryW = 1/7  # base gamma_w=1/20 =0.05
recoveryM = 1/8  # base gamma_m=0.05

beta_w = 0.8*(1/7)/N # transmission rate wild type 1.2
phi_w = 1.6 # change of infectivity Aw versus the infectivity of Hw

beta_m = 0.8*(1/7)/N # transmission rate mutant 1.4
phi_m = 1.6 # change of infectivity Am versus the infectivity of Hm

eta_w =0.4 # fraction that develops severe disease after wild-type infection 
eta_m =0.4 # fraction that develops severe disease after mutant infection

eps = 1e-3 # fraction of infections that mutate
loss_immunity = 1/60 # a/alpha = mean duration of natural immunity

rho = N*1e-6 # spillover transmission rate 




##############################################################
############        code chunk 2             #################
########### Creating a GillesPy2 Model #######################
##############################################################

# Import the types that'll be needed to define your Model.
from gillespy2.core import (
    Model,
    Species,
    Reaction,
    Parameter
)



class two_strain(Model):
     def __init__(self, parameter_values=None):

            
            Model.__init__(self, name="two_strain")
            
            """
            Parameters
            """
            
            # spillover
            p1 = Parameter(name="p1", expression=rho*beta_w*eta_w*(1-eps))
            p2 = Parameter(name="p2", expression=rho*beta_w*(1-eta_w)*(1-eps))
            p3 = Parameter(name="p3", expression=rho*beta_w*eta_w*eps)
            p4 = Parameter(name="p4", expression=rho*beta_w*(1-eta_w)*eps)
            
            # Infection (Hospitalized wild-type)
            p5 = Parameter(name="p5", expression=beta_w*eta_w*(1-eps))
            p6 = Parameter(name="p6", expression=beta_w*eta_w*eps)
            p7 = Parameter(name="p7", expression=beta_w*(1-eta_w)*(1-eps))
            p8 = Parameter(name="p8", expression=beta_w*(1-eta_w)*eps)
            
            # Infection (Asymptomatic wild-type)
            p9 = Parameter(name="p9", expression=beta_w*phi_w*eta_w*(1-eps))
            p10 = Parameter(name="p10", expression=beta_w*phi_w*eta_w*eps)
            p11 = Parameter(name="p11", expression=beta_w*phi_w*(1-eta_w)*(1-eps))
            p12 = Parameter(name="p12", expression=beta_w*phi_w*(1-eta_w)*eps)
            
            # Infection (hospitalized mutant)
            p13 = Parameter(name="p13", expression=beta_m*eta_m)
            p14 = Parameter(name="p14", expression=beta_m*(1-eta_m))
            
            # Infection (hospitalized mutant)
            p15 = Parameter(name="p15", expression=beta_m*phi_m*eta_m)
            p16 = Parameter(name="p16", expression=beta_m*phi_m*(1-eta_m))
            
            gammaW = Parameter(name="gammaW", expression=recoveryW)
            gammaM = Parameter(name="gammaM", expression=recoveryM)
            alpha = Parameter(name="alpha", expression=loss_immunity)
            
            # Add the Parameters to the Model.
            self.add_parameter([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,gammaW,gammaM,alpha])
            
            """
            Species 
            """
            
            S = Species(name="S", initial_value=S0)
            Hw = Species(name="Hw", initial_value=Hw0)
            Aw = Species(name="Aw", initial_value=Aw0)
            Hm = Species(name="Hm", initial_value=Hm0)
            Am = Species(name="Am", initial_value=Am0)
            R = Species(name="R", initial_value=R0)
            
            # Add the Species to the Model.
            self.add_species([S, Hw, Aw, Hm, Am, R])
            
            """
            Reactions
            """
            
            # Spillover reactions
            r1 = Reaction(
                    name="spilloverHw",
                    reactants={S: 1}, 
                    products={S: 0, Hw: 1},
                    rate=p1
                )
            
            r2 = Reaction(
                    name="spilloverAw",
                    reactants={S:1}, 
                    products={S: 0, Aw: 1},
                    rate=p2
                )
            
            r3 = Reaction(
                    name="spilloverHm",
                    reactants={S: 1}, 
                    products={S: 0, Hm: 1},
                    rate=p3
                )
            
            r4 = Reaction(
                    name="spilloverAm",
                    reactants={S:1,}, 
                    products={S: 0, Am: 1},
                    rate=p4
                )
            
            # Infection (hospitalized wild-type)
            r5 = Reaction(
                    name="infectionHwHw",
                    reactants={S: 1, Hw: 1}, 
                    products= {S: 0, Hw: 2},
                    rate=p5
                )
            
            r6 = Reaction(
                    name="infectionHwHm",
                    reactants={S: 1, Hw:1}, 
                    products= {S: 0, Hw: 1, Hm: 1},
                    rate=p6
                )
            
            
            r7 = Reaction(
                    name="infectionHwAw",
                    reactants={S:1, Hw: 1}, 
                    products={S: 0, Hw: 1, Aw: 1},
                    rate=p7
                )
            
            r8 = Reaction(
                    name="infectionHwAm",
                    reactants={S:1, Hw: 1}, 
                    products={S: 0, Hw: 1, Am: 1},
                    rate=p8
                )
            
            # Infection (asymptomatic wild-type)
            r9 = Reaction(
                    name="infectionAwHw",
                    reactants={S: 1, Aw: 1},
                    products={S: 0, Aw: 1, Hw: 1},
                    rate=p9
                )
            
            
            r10 = Reaction(
                    name="infectionAwHm",
                    reactants={S: 1, Aw: 1},
                    products={S: 0, Aw: 1, Hm: 1},
                    rate=p10
                )
            
            r11 = Reaction(
                    name="infectionAwAw",
                    reactants={S: 1, Aw: 1},
                    products={S: 0, Aw: 2},
                    rate=p11
                )
            
            r12 = Reaction(
                    name="infectionAwAm",
                    reactants={S: 1, Aw: 1},
                    products={S: 0, Aw: 1, Am: 1},
                    rate=p12
                )
            
            
            # Infection (hospitalized mutant)
            r13 = Reaction(
                    name="infectionHmHm",
                    reactants={S: 1, Hm: 1}, 
                    products={S: 0, Hm: 2},
                    rate=p13
                )
            
            
            r14 = Reaction(
                    name="infectionHmAm",
                    reactants={S: 1, Hm: 1}, 
                    products={S: 0, Hm: 1, Am: 1},
                    rate=p14
                )
            
            # Infection (asymptomatic mutant)
            r15 = Reaction(
                    name="infectionAmHm",
                    reactants={S: 1, Am: 1}, 
                    products={S: 0, Am: 1, Hm: 1},
                    rate=p15
                )
            
            
            r16 = Reaction(
                    name="infectionAmAm",
                    reactants={S: 1, Am: 1}, 
                    products={S: 0, Am: 2},
                    rate=p16
                )
            
            # recovery
            
            r17 = Reaction(
                    name="recoveryHw",
                    reactants={Hw: 1}, 
                    products={Hw: 0, R: 1},
                    rate=gammaW
                )
            
            r18 = Reaction(
                    name="recoveryAw",
                    reactants={Aw: 1}, 
                    products={Aw: 0, R: 1},
                    rate=gammaW
                )
            
            r19 = Reaction(
                    name="recoveryHm",
                    reactants={Hm: 1}, 
                    products={Hm: 0, R: 1},
                    rate=gammaM
                )
            
            r20 = Reaction(
                    name="recoveryAm",
                    reactants={Am: 1}, 
                    products={Am: 0, R: 1},
                    rate=gammaM
                )
            
            # loss of immunity
            r21 = Reaction(
                    name="loss_immunity",
                    reactants={R: 1}, 
                    products={R: 0, S: 1},
                    rate=alpha
                )
            
            self.add_reaction([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21])
            
            # Use NumPy to set the timespan of the Model.
            self.timespan(np.linspace(0, 300, 301))
            
            
            
model = two_strain()



##############################################################
################   code chunk  3          ####################
#################### Solving the Model #######################
##############################################################



# ---------------------------------------Stochastic solution of the model
from gillespy2.solvers.numpy import TauLeapingSolver

num_traj = 500
stochastic_sol = model.run(solver=TauLeapingSolver, number_of_trajectories=num_traj, seed=10)


# ---------------------------------------Deterministic solution

    
from gillespy2.solvers.numpy import ODESolver

deterministic_sol = model.run(solver=ODESolver)



##############################################################
################   code chunk  4          ####################
############  Plot the results of the simulations ############
##############################################################



# ------------------------------------------------

plt.style.use('seaborn')


#-----------------------------------------------------------------------------
# Hospitalized
figH, (ax1, ax2, ax3) = plt.subplots(figsize=(20,5), ncols=3, nrows=1, constrained_layout=True)
# Asymptomatic
#figA, ax2 = plt.subplots(figsize=(8,5), tight_layout=True)


# stochastic realizations
#for index in range(0, num_traj-1):
for index in range(0, num_traj):
    trajectory = stochastic_sol[index]
    ax1.plot(trajectory['time'], trajectory['Hw'], 'seagreen', alpha=0.2)
    ax1.plot(trajectory['time'], trajectory['Hm'], 'orange', alpha=0.2)
    ax2.plot(trajectory['time'], trajectory['Aw'], 'firebrick', alpha=0.2)
    ax2.plot(trajectory['time'], trajectory['Am'], 'dodgerblue', alpha=0.2)
    ax3.plot(trajectory['time'], trajectory['Hw']+trajectory['Aw'], 'silver', alpha=0.2)
    ax3.plot(trajectory['time'], trajectory['Hm']+trajectory['Am'], 'pink', alpha=0.2)


# deterministic solution (mean field approximation)
ax1.plot(deterministic_sol['time'],deterministic_sol['Hw'],'green', linewidth=2, label=r'$H_w$')
ax1.plot(deterministic_sol['time'],deterministic_sol['Hm'], 'darkorange', linewidth=2, label=r'$H_m$')
ax2.plot(deterministic_sol['time'],deterministic_sol['Aw'],'firebrick', linewidth=2, label=r'$A_w$')
ax2.plot(deterministic_sol['time'],deterministic_sol['Am'], 'blue', linewidth=2, label=r'$A_m$')
ax3.plot(deterministic_sol['time'],deterministic_sol['Hw']+deterministic_sol['Aw'], 'gray', linewidth=2, label=r'$H_{w}+A_{w}$')
ax3.plot(deterministic_sol['time'],deterministic_sol['Hm']+deterministic_sol['Am'], 'orchid', linewidth=2, label=r'$H_{m}+A_{m}$')


# Customise some display properties

ax1.set_xlabel('Time (days)', fontsize=18)
ax1.set_ylabel('Individuals', fontsize=18)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', labelsize=18)
#ax1.xaxis.set_major_locator(MaxNLocator(nbins=9))
ax1.yaxis.set_major_locator(MaxNLocator(nbins=6))
#ax1.set_title(r'$\beta = \gamma$, $\phi=1.5$, $\tau=1e-4$')
ax1.legend(loc='best', fontsize=18)
#ax1.set_ylim([0, 5000])

ax2.set_xlabel('Time (days)', fontsize=18)
#ax2.set_ylabel('Individuals', fontsize=18)
ax2.tick_params(axis='x', labelsize=18)
ax2.tick_params(axis='y', labelsize=18)
#ax1.xaxis.set_major_locator(MaxNLocator(nbins=9))
ax2.yaxis.set_major_locator(MaxNLocator(nbins=6))
#ax1.set_title(r'$\beta = \gamma$, $\phi=1.5$, $\tau=1e-4$')
ax2.legend(loc='best', fontsize=18)
#ax2.set_ylim([0, 5000])

ax3.set_xlabel('Time (days)', fontsize=18)
#ax2.set_ylabel('Individuals', fontsize=18)
ax3.tick_params(axis='x', labelsize=18)
ax3.tick_params(axis='y', labelsize=18)
#ax1.xaxis.set_major_locator(MaxNLocator(nbins=9))
ax3.yaxis.set_major_locator(MaxNLocator(nbins=6))
#ax1.set_title(r'$\beta = \gamma$, $\phi=1.5$, $\tau=1e-4$')
ax3.legend(loc='best', fontsize=18)
#ax3.set_ylim([0, 5000])



#plt.savefig("GammaM8days.pdf", bbox_inches = 'tight')
plt.show()







