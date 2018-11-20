#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:29:26 2018

@author: martin.luttmann
"""

from Cranknicolson import *
import matplotlib.pyplot as plt 
import numpy as np
from math import *
import time



#initial condition
def gaussian(x, x0, sigma, U_0):
    return U_0 * np.exp(-(x-x0)**2/(sigma**2))/(2*sigma*sqrt(2*pi))

# Those will not change
x_min, x_max, nx = -20, 20., 1000
t_min, t_max, nt = 0, 10., 100

# Construct the general Crank-Nicolson solver
crank = CrankNicolson2()
crank.set_grid(x_min, x_max, nx, t_min, t_max, nt)

E0=0.01
omega=1.
#potentiel V(x,t)


def Vc(x): #soft coulomb potential
    return -1/sqrt(2+x**2)
    

def E(t):
    return E0*cos(omega*t)
def V(t,x):
    return Vc(x)-x*E(t)

delta_x= (x_max-x_min)/nx
delta_t=(t_max-t_min)/nt

#potentiel discretise
    #tableau des potentiels aux differents instants t
V_matrix=np.zeros([nt, nx])
for i in range(nt):
    for p in range(nx):
        V_matrix[i][p]=V(t_min+i*delta_t, x_min+p*delta_x)
    



#f = lambda u: mu*u - u**3

def f(psi,n):
    result=psi
    for i in range(nx):
        result[i]=psi[i]*V_matrix[n][i]
    return result
    
    
crank.set_parameters(f)

s, x0, u0 = 1., 0., 1.
u = gaussian(crank.x_pts, x0, s, u0)
crank.solve(u, boundary_conditions=['neumann', 'neumann'])

plt.figure()
plt.plot(crank.x_pts, crank.get_first_psi().real, '-b', lw=2)
plt.plot(crank.x_pts, crank.get_first_psi().imag, '-r', lw=2)
plt.title('$t=0$')


plt.figure()
plt.plot(crank.x_pts, crank.get_second_psi().real, '-b', lw=2)
plt.plot(crank.x_pts, crank.get_second_psi().imag, '-r', lw=2)
plt.title('$t=dt$')

plt.figure()
plt.plot(crank.x_pts, crank.get_final_psi().real, '-b', lw=2)
plt.plot(crank.x_pts, crank.get_final_psi().imag, '-r', lw=2)
#plt.plot(crank.x_pts, Vc(), '--k', lw=1)
plt.title('$t=t_{max}$')
#plt.plot(crank.x_pts, abs(crank.get_final_psi())**2, '-k', lw=2)

plt.show()

Norms=[]
for i in range(nt):
    Norms.append(np.linalg.norm(abs(crank.psi_matrix[i,:])**2))
    
    

plt.figure()
plt.plot(crank.t_pts,Norms)
plt.show()


#plt.pcolormesh(crank.x_pts, crank.t_pts, crank.psi_matrix.real)