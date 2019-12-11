# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 08:19:23 2019

@author: caisa
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

#andrew.lifson@thep.lu.se

#For double well, this will look very similar but change a few parameters

# Parameters
ximin=-16
ximax=-ximin
Nsteps=2001    #increasing from 1001
nu0=4.0
xi0 = 1.1

#Defining our vectors
xi = np.linspace(ximin, ximax, Nsteps)
phi = np.zeros(Nsteps)
h = xi[1]-xi[0]
nu = -nu0*(np.exp(-(xi-xi0)**2) + np.exp(-(xi+xi0)**2))  #adjusted to two potential wells

#energy guess, bisection method
    
def phii(eps):
    k = np.sqrt(-eps)
    phi[1] = np.exp(k*h)
    for i in np.arange(2, Nsteps):
        phi[i] = (2+h**2*(nu[i-1]-eps))*phi[i-1]-phi[i-2]
    return phi.tolist()

def bisec(emin=-0.01, emax=-0):    
    while (np.abs((emin-emax)/(emin+emax))>1e-8):   #relative error
        e0=(emin+emax)/2
        if ((phii(e0)[-1]*phii(emax)[-1])>0):
            emax = e0
        else:
            emin = e0
    return (print('epsilon is ', e0), e0)

plt.figure()
plt.plot(xi, phii(bisec()[1]))
plt.xlabel('xi')
plt.ylabel('phi')
plt.grid()
plt.show()
    
"""
Questions 2.1:
1. when xi0 = 1000, we find none. Since potential decrease with decreasing xi0, the energy states will be next to zero
and thus finding no valid ground state.
2. 2;   ground:  -1.1221 eV
        1st excite: -0.7814 eV
3. since states are splitting when wells are brought together, we would expect double the states from single well, but since
the barrier height was halfed, we got half the states, 2.

Questions 2.2
1. three (fourth disappeared)
2. skriv in värden här
3. ??
4. higher states give more peaks to interrct and overlap,
thus the ratio och overlaping should be larger at excited states
5. 

"""





