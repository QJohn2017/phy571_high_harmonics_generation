#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:09:38 2018
@author: martin.luttmann
"""

from matplotlib import pyplot as plt
import scipy
import numpy as np
import scipy.sparse
import scipy.linalg as la



########### A class that solves Schrodingers equation ###########
class CrankNicolson:
    
    
    def set_grid(self, x_min, x_max, n_x, t_min, t_max, n_t):

        self.x_min, self.x_max, self.n_x = x_min, x_max, n_x
        self.t_min, self.t_max, self.n_t = t_min, t_max, n_t
        self.x_pts, self.delta_x = np.linspace(x_min, x_max, n_x, retstep=True, endpoint=False)
        self.t_pts, self.delta_t = np.linspace(t_min, t_max, n_t, retstep=True, endpoint=False)
        self.position = np.zeros(self.n_t)*1j #This array contains the average position over time
        self.dipole_acceleration = np.zeros(self.n_t)*1j #This array contains the average dipole acceleration over time
    def set_parameters(self, f):
        
        self.f =  f
    
    
    # Returns the 1D laplacian of the wavefunction 
    def get_laplacian(self, psi):
        
        L = np.zeros(self.n_x)*1j
        for i in range(1,self.n_x-1):
            L[i] = (psi[i+1]+psi[i-1]-2*psi[i])
        return L/self.delta_x**2


    def solve(self, psi_init, V_matrix, E_matrix, boundary_conditions):
            
        sig = 1j* self.delta_t / (4.* self.delta_x**2)
        
      
        data_type = type(sig*psi_init[0])
        
        
        # Matrix of the wavefunctions at differents t: each column of psi_matrix is a wavefunction at fixed t
        self.psi_matrix = np.zeros([self.n_t, self.n_x], dtype=data_type)  
  

        # Calculating the matrices of the linear system 
        A = self._fillA_sp(sig, self.n_x, data_type)
        B = self._fillB_sp(sig, self.n_x, data_type)
            
        
        
        # Setting boundary conditions
        for b in [0,1]:
            if boundary_conditions[b] == 'dirichlet':
                # u(x,t) = 0
                A[1,-b] = 1.0
                A[2*b,1-3*b] = 0.0
                B[-b,-b] = 0.0
                B[-b,1-3*b] = 0.0
            elif boundary_conditions[b] == 'neumann':
                # u'(x,t) = 0
                A[2*b,1-3*b] = -2*sig
                B[-b,1-3*b] = 2*sig
                
        ##### Propagating the solution #####
        psi = psi_init
        fpsi_1=self.f(psi,0)
        for n in range(self.n_t):
            self.psi_matrix[n,:] = psi
            Laplacian_psi = self.get_laplacian(psi)
            self.position[n] = scipy.integrate.simps(psi.conjugate()*self.x_pts*psi, dx = self.delta_x) #compute the average position at time n*delta_t
            self.dipole_acceleration[n] = scipy.integrate.simps(psi.conjugate()*-1*np.gradient(V_matrix[n], self.delta_x)*psi+E_matrix[n], dx = self.delta_x) #compute the average dipole acceleration at time n*delta_t
            
            fpsi_2=self.f(psi,n)
            
            psi = la.solve_banded((1,1),A, B.dot(psi) - 1j*self.delta_t * (3/2*fpsi_2-1/2*fpsi_1),
                                check_finite=False)
            fpsi_1=fpsi_2
            
    def get_final_psi(self):
        
        return self.psi_matrix[-1,:].copy()
    
    
    def get_first_psi(self):
       
        return self.psi_matrix[0,:].copy()
    
    def get_second_psi(self):
        return self.psi_matrix[1,:].copy()
        

    ####### Returns a tridiagonal matrix in compact form ab[1+i-j,j]=a[i,j] #######
    def _fillA_sp(self, sig, n, data_type):
        
        A = np.zeros([3,n], dtype=data_type) # A has three diagonals and size n
        A[0] = -sig # superdiagonal
        A[1] = 1+2*sig # diagonal
        A[2] = -sig # subdiagonal
        return A

    ####### Returns a tridiagonal sparse matrix in csr-form #######
    def _fillB_sp(self, sig, n, data_type):
        
        _o = np.ones(n, dtype=data_type)
        supdiag = (sig)*_o[:-1]
        diag = (1-2*sig)*_o
        subdiag = (sig)*_o[:-1]
        return scipy.sparse.diags([supdiag, diag, subdiag], [1,0,-1], (n,n), format="csr")