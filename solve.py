#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:29:26 2018

@author: martin.luttmann
"""

from Cranknicolsoneasy import *
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import numpy as np
from scipy.constants import epsilon_0, c
from math import *
import time




    
x_min, x_max, nx = -50, 50., 1024
t_min, t_max, nt = 0, 1000., 10000


# Construct the general Crank-Nicolson solver
crank = CrankNicolson2()
crank.set_grid(x_min, x_max, nx, t_min, t_max, nt)
tau=800.
E0=10**(-2)
omega=0.05
#potentiel V(x,t)

Uc=1.


def Vc(x): #soft coulomb potential
    return -Uc/sqrt(2+x**2)
    #return -x**2
    
#field should be zero and not varying at t=0
def E(t):
    #return E0*np.sin(omega*t)
    if t<tau:
        return E0*np.sin(t/tau*pi)**2*np.sin(omega*t)
    else:
        return np.zeros_like(t)
    
    
def V(t,x):
    return Vc(x)-x*E(t)
    

#delta_x= (x_max-x_min)/nx
    
delta_t=crank.delta_t
delta_x=crank.delta_x
#delta_t=(t_max-t_min)/nt

#potentiel discretise
    #tableau des potentiels aux differents instants t
V_matrix=np.zeros([nt, nx])
for i in range(nt):
    for p in range(nx):
        V_matrix[i][p]=V(t_min+i*delta_t, x_min+p*delta_x)
    



#f = lambda u: mu*u - u**3

def f(psi,n):
    return psi*V_matrix[n]
    #result=psi
    #for i in range(nx):
     #   result[i]=psi[i]*V_matrix[n][i]
    #return result
    
    
crank.set_parameters(f)

s, x0, u0 = 1., 0., 1.
u = np.load('Groundstate.npy')
#N=np.trapz(np.trapz(abs(u)**2,crank.x_pts))
#print(N)
#for elt in u:
#    elt=elt/N


crank.solve(u, V_matrix, boundary_conditions=['dirichlet', 'dirichlet'])

plt.figure()
plt.plot(crank.x_pts, crank.get_first_psi().real, '-b', lw=2)
plt.plot(crank.x_pts, crank.get_first_psi().imag, '-r', lw=2)
plt.title('$t=0$')


#plt.figure()
#plt.plot(crank.x_pts, crank.get_second_psi().real, '-b', lw=2)
#plt.plot(crank.x_pts, crank.get_second_psi().imag, '-r', lw=2)
#plt.title('$t=dt$')

plt.figure()
plt.plot(crank.x_pts, abs(crank.get_first_psi())**2, '-b', lw=2)
plt.plot(crank.x_pts, abs(crank.get_final_psi())**2, '-r', lw=2)
#plt.plot(crank.x_pts, Vc(), '--k', lw=1)
plt.title('$t=t_{max}$')
#plt.plot(crank.x_pts, abs(crank.get_final_psi())**2, '-k', lw=2)

plt.show()

Norms=[]
logNorms=[]
for i in range(nt):
    Norms.append(np.trapz(abs(crank.psi_matrix[i])**2,crank.x_pts))
    #logNorms.append(log(np.linalg.norm(crank.psi_matrix[i,:])**2))
    

plt.figure()
plt.semilogy(crank.t_pts,Norms)
plt.xlabel('$t \  (a.u.)$', fontsize='large')
plt.ylabel('$\psi^2 $',fontsize='large')
plt.title('Norm of the wavefunction')

plt.figure()
plt.semilogy(crank.t_pts,crank.position)
plt.xlabel('$t \  (a.u.)$', fontsize='large')
plt.ylabel('$\psi^2 $',fontsize='large')
plt.title('Average position')

plt.figure()
plt.semilogy(2*np.pi*np.fft.fftfreq(crank.n_t, crank.delta_t), np.fft.fft(crank.dipole_acceleration))
plt.xlabel('$t \  (a.u.)$', fontsize='large')
plt.ylabel('$\psi^2 $',fontsize='large')
plt.title('spectrum')


plt.figure()
#plt.plot(crank.t_pts,logNorms)
plt.pcolormesh(crank.x_pts, crank.t_pts, abs(crank.psi_matrix)**2, norm=LogNorm(vmin=1e-3,vmax=1))
plt.xlabel('$x \  (a.u.)$', fontsize='large')
plt.ylabel('$t \ (a.u.)$', fontsize='large')
plt.title('Squared modulus of the wavefunction')
#plt.figure()
#plt.pcolormesh(crank.x_pts, crank.t_pts, crank.psi_matrix.imag)

plt.figure()
plt.plot(crank.t_pts, [E(t) for t in crank.t_pts], '-k', lw=2)
plt.xlabel('$t \ (a.u.)$', fontsize='large')
plt.ylabel('$E(t) \ (a.u.)$', fontsize='large')
plt.title('Electric field pulse')

