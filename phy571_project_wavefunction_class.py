#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 14:16:26 2018

@author: martin.guillot
"""

import numpy as np
import scipy.integrate as scp

class wavefunction:
    """create the wave function class"""
    def __init__(self, N, T, xmin, xmax, t0, tfinal, psi0):
        """create a wavefunction object with spatial boundaries [xmin,xmax], temporal boundaries [t0, tfinal] and a spatial/temporal precision of N/T. At t = t0, the wavefunction is equal to psi0"""
        self.N = N
        self.T = T
        self.xmin = xmin
        self.xmax = xmax
        self.t0 = t0
        self.tfinal = tfinal
        self.X, self.dx = np.linspace(xmin, xmax, N, retstep = True)
        self.K = np.fft.fftfreq(N, self.dx)
        psi = np.zeros(T, N)
        psi[0] = psi0
        self.psi = psi """a T x N array representing the wavefunction during the experiment"""
    
    
    def get_position(self):
        """return the average position over time"""
        position = np.zeros(self.T)
        X = np.linspace(self.xmin, self.xmax, self.N)
        for i in range(self.T):
            position[i] = scp.simps(self.phi[i].conjugate()*X*self.phi[i], dx = self.dx)
        return position
    
    def get_laplacian(self):
        """terurn the 1D laplacian of the wavefunction"""
        return (self.phi[:][:2]+self.phi[:][:-2]-2*self.phi[:][1:-1])/self.dx**2
    
    def get_energy(self, V, E):
        """return the average energy over time""" 
        energy = np.zeros(self.T)
        X = np.array([self.xmin+i*self.dx for i in range(N)])
        for i in range(self.T):
                energy[i] = scp.simps(self.phi[i].conjugate()*(-1/2*self.get_laplacian()[i]+V*self.psi[i]-X*E[i]*self.psi[i]), dx = self.dx)
        return energy
    
    def get_dipole_acceleration(self, V, E):
        """return the dipole acceleratetion over time for a time-independant potential V and a slowly variating uniform E field """
        A = np.zeros(self.T)
        dx = (self.xmax-self.xmin)/self.N
        for i in range(self.T):
                A[i] = scp.simps(self.psi[i].conjugate()*(-np.gradient(V, dx))*self.psi[i], dx = self.dx)+E[i]
        return A
    
    def get_spectrum(self, V, E):
        """return the spectrum of the radiation emitted by the electron by fourrier-transforming the dipole acceleration"""
        S = np.fft.fft(self.get_dipole_acceleration(V, E))
        return S
                       
        
