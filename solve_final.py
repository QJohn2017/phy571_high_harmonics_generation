#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:29:26 2018

@author: martin.luttmann
"""

from Cranknicolson_final import *
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import numpy as np





########### Parameters of the time-space grid ###########
x_min, x_max, nx = -50, 50., 1024
t_min, t_max, nt = 0, 1000., 10000


########### Construct the general Crank-Nicolson solver ###########
crank = CrankNicolson()
crank.set_grid(x_min, x_max, nx, t_min, t_max, nt)
tau=800.
E0=10**(-2)
omega=0.05





########### V(x,t) potential (atomic and electric potential) ###########

#soft coulomb potential
Uc=1.
def Vc(x): #soft coulomb potential
    return -Uc/np.sqrt(2+x**2)
    
    
#Elctric field pulse
def E(t):
    
    if t<tau:
        return E0*np.sin(t/tau*np.pi)**2*np.sin(omega*t)
    else:
        return np.zeros_like(t)
    
#global effective potential on the electron   
def V(t,x):
    return Vc(x)-x*E(t)
    


    
delta_t=crank.delta_t
delta_x=crank.delta_x


#potentiel discretise
    #tableau des potentiels aux differents instants t
V_matrix=np.zeros([nt, nx])
for i in range(nt):
    for p in range(nx):
        V_matrix[i][p]=V(t_min+i*delta_t, x_min+p*delta_x)
    
########### Second term f in the linear system ###########
def f(psi,n):
    return psi*V_matrix[n]
    
    
    
crank.set_parameters(f)

########### Solving Schrodinger equation ###########

#loading initial wavefunction
u = np.load('Groundstate.npy')
crank.solve(u, boundary_conditions=['dirichlet', 'dirichlet'])



#Plotting initial wavefunction
plt.figure()
plt.plot(crank.x_pts, crank.get_first_psi().real, '-b', lw=2)
plt.plot(crank.x_pts, crank.get_first_psi().imag, '-r', lw=2)
plt.title('Wavefunction at $t=0$')


#Plotting initial and final squared modulus of the wavefunction
plt.figure()
plt.plot(crank.x_pts, abs(crank.get_first_psi())**2, '-b', lw=2)
plt.plot(crank.x_pts, abs(crank.get_final_psi())**2, '-r', lw=2)
plt.title('$t=t_{max}$')
plt.show()



#Recording the norm of the wavefunction in time
Norms=[]
logNorms=[]
for i in range(nt):
    Norms.append(np.trapz(abs(crank.psi_matrix[i])**2,crank.x_pts))
    
plt.figure()
plt.plot(crank.t_pts,Norms)
plt.xlabel('$t \  (a.u.)$', fontsize='large')
plt.ylabel('$\psi^2 $',fontsize='large')
plt.title('Norm of the wavefunction')



#Plotting the squared modulus of the wavefunction with respect to t
plt.figure()
plt.pcolormesh(crank.x_pts, crank.t_pts, abs(crank.psi_matrix)**2, norm=LogNorm(vmin=1e-3,vmax=1))
plt.xlabel('$x \  (a.u.)$', fontsize='large')
plt.ylabel('$t \ (a.u.)$', fontsize='large')
plt.title('Squared modulus of the wavefunction')



#Plotting the shape of the laser pulse
plt.figure()
plt.plot(crank.t_pts, [E(t) for t in crank.t_pts], '-k', lw=2)
plt.xlabel('$t \ (a.u.)$', fontsize='large')
plt.ylabel('$E(t) \ (a.u.)$', fontsize='large')
plt.title('Electric field pulse')