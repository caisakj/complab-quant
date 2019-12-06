# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 13:32:35 2019

@author: caisa
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

#andrew.lifson@thep.lu.se

# Parameters
ximin=-8
ximax=-ximin
Nsteps=1001
nu0=4.0

#Defining our vectors
xi = np.linspace(ximin, ximax, Nsteps)
phi = np.zeros(Nsteps)
h = xi[1]-xi[0]
nu = -nu0*np.exp(-xi**2)

#Energy guess, nu0<eps<0, by hand
#eps = -2.37547051  #crisp

#define more parameters
k = np.sqrt(-eps)
phi[0] = 1
phi[1] = np.exp(k*h)

#eps = -2.375470

#For loop time!
for i in np.arange(2, Nsteps):
    phi[i] = (2+h**2*(nu[i-1]-eps))*phi[i-1]-phi[i-2]

#energy guess, bisection method
    
def phii(eps):
    k = np.sqrt(-eps)
    phi[1] = np.exp(k*h)
    for i in np.arange(2, Nsteps):
        phi[i] = (2+h**2*(nu[i-1]-eps))*phi[i-1]-phi[i-2]
    return phi.tolist()

def bisec(emin=-3.5, emax=-2):
    while (np.abs((emin-emax)/(emin+emax))>1e-8):   #relative error
        e0=(emin+emax)/2
        if ((phii(e0)[-1]*phii(emax)[-1])>0):
            emax = e0
        else:
            emin = e0
    return (print('epsilon is ', e0), e0)

#print(bisec()[1])
#plt.figure()
#plt.plot(xi, phii(bisec()[1]))
#plt.xlabel('xi')
#plt.ylabel('phi')
#plt.grid()
#plt.show()

#To check excited states, change emin and emax cleverely
#bisec(-2.3754704552404284, 0)
#plt.figure()
#plt.plot(xi, phii(bisec(-2.3754704552404284, 0)[1]))
#plt.xlabel('xi')
#plt.ylabel('phi')
#plt.grid()
#plt.show()

"""
Questions:
1. 0 nodes
2. 1 node
3. larger energy --> larger eps --> smaller k --> smaller exponent that kills wavefunction
4. doesn't go above 0, make k imaginary - only 1 excited state
"""

#bisec(-2.3754704552404284, 0)
#bisec(-0.1498260235979687, 0)














