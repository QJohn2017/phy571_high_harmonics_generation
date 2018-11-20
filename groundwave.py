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
    def __init__(self, N, dt, xmin, xmax, V):    
        self.N = N 
        self.dt = dt #temporal resolution"""
        self.xmin = xmin #inferior spatial boundary
        self.xmax = xmax #superior spatial boundary
        self.V = V #An arbitrary potential with at least one bounded state, written as a size N array
        self.X, self.dx = np.linspace(self.xmin, self.xmax, N, retstep = True)
        self.K = np.fft.fftfreq(N, self.dx)
        self.Vk = np.fft.fft(self.V)
    
    def propagate(self, phi):
        """Takes an arbitrary wavefunction as a size N array, fourrier-transforms it, applies a wick rotation then inverse-fourrier-transforms it"""
        phik = np.fft.fft(phi)
        phiksec = np.exp(-1j*(1/2*self.K**2)*self.dt)*phik #here is the wick rotation
        phisec = np.fft.ifft(phiksec)
        phisec = phisec*np.exp(-1j*self.V*self.dt)
        phisec = self.renormalize(phisec)
        return phisec

    def renormalize(self, phi):
        return phi/scp.simps(np.abs(phi)**2, dx = self.dx)**(1/2)
    

Xtest = np.linspace(-5, 5, 1024)
Vtest = 1/4*Xtest**4
Test = Groundstate(1024, (1-0.1j)*0.00001,-5, 5, Vtest)
phitest = np.cos(np.pi/10*Xtest)**2
#phitest = np.exp(-10*Xtest**2)

for i in range(1000):
    phitest = Test.propagate(phitest)

plt.plot(Xtest, 100*np.abs(phitest)**2)
plt.plot(Xtest, Vtest)
plt.show()
