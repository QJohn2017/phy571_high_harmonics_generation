# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:42:21 2018

@author: marti
"""

import numpy as np
import scipy.integrate as scp
import matplotlib.pyplot as plt

class Groundstate:
    """This class allows you to calculate the groundstate wavefunction of a particle inside a potential V"""
    def __init__(self, N, dt, xmin, xmax, V, phi):    
        self.N = N 
        self.dt = dt #temporal resolution"""
        self.xmin = xmin #inferior spatial boundary
        self.xmax = xmax #superior spatial boundary
        self.V = V #An arbitrary potential with at least one bounded state, written as a size N array
        self.X, self.dx = np.linspace(self.xmin, self.xmax, N, retstep = True)
        self.K = 2*np.pi*np.fft.fftfreq(N, self.dx)
        self.Vk = np.fft.fft(self.V)
        self.phi = phi
    
    def propagate(self):
        """Takes an arbitrary wavefunction as a size N array, fourrier-transforms it, applies a wick rotation then inverse-fourrier-transforms it"""
        phik = np.fft.fft(self.phi)
        phiksec = np.exp(-1j*(1/2*self.K**2)*self.dt)*phik #here is the wick rotation
        phisec = np.fft.ifft(phiksec)
        phisec = phisec*np.exp(-1j*self.V*self.dt)
        phisec = self.renormalize(phisec)
        self.phi = phisec

    def renormalize(self, phi):
        return phi/scp.simps(np.abs(phi)**2, dx = self.dx)**(1/2)
    
    def get_laplacian(self):
        """terurn the 1D laplacian of the wavefunction"""
        L = np.zeros(self.N)
        for i in range(self.N-1):
            L[i] = (self.phi[i+1]+self.phi[i-1]-2*self.phi[i])
        return L/self.dx**2
        #return (self.phi[-2:]+self.phi[:-2]-2*self.phi[1:-1])/self.dx**2
    
    def get_energy(self):
        """return the average energy over time""" 
        energy = scp.simps(self.phi.conjugate()*(-1/2*self.get_laplacian()+self.V*self.phi), dx = self.dx)
        return energy
    

Xtest = np.linspace(-20,20, 4096)
Vtest = 1/2*Xtest**2
phitest = np.exp(-0.01*Xtest**2)+0j
Test = Groundstate(4096, 1e-3*(1-0.7j),-20, 20, Vtest, phitest)
#phitest = np.cos(np.pi/10*Xtest)**2

S = 1000
E = np.zeros(S)
for i in range(S):
    Test.propagate()
    E[i] = Test.get_energy()

    
plt.plot(Xtest, 100*np.abs(Test.phi)**2)
plt.plot(Xtest, Vtest)
plt.show()
plt.plot([k for k in range(S)], E)
plt.show()
print(E[-2])
