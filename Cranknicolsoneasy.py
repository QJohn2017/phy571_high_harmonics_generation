#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:09:38 2018

@author: martin.luttmann
"""

from matplotlib import pyplot as plt
from math import *



import numpy as np
import scipy.sparse
import scipy.linalg as la




class CrankNicolson2:
    """A class that solves Schrodingers equation"""
    
    def set_grid(self, x_min, x_max, n_x, t_min, t_max, n_t):

        self.x_min, self.x_max, self.n_x = x_min, x_max, n_x
        self.t_min, self.t_max, self.n_t = t_min, t_max, n_t
        self.x_pts, self.delta_x = np.linspace(x_min, x_max, n_x, retstep=True, endpoint=False)
        self.t_pts, self.delta_t = np.linspace(t_min, t_max, n_t, retstep=True, endpoint=False)
        
    def set_parameters(self, f):
        
        self.f =  f

#get next psi
    def solve(self, psi_init, boundary_conditions):
            
        sig = 1j* self.delta_t / (4.* self.delta_x**2)
        
        #norm_psi=np.trapz(abs(psi_init)**2,self.x_pts)
        #psi_init/=np.sqrt(norm_psi)
        
        
        # Figure the data type
        data_type = type(sig*psi_init[0])
        
        self.psi_matrix = np.zeros([self.n_t, self.n_x], dtype=data_type)  #matrice des fonctions d'onde
  


        '''
            
        A = self._make_tridiag(sig, self.n_x, data_type)
        B = self._make_tridiag(-sig,  self.n_x, data_type)
        # Set boundary conditions
        for b in [0,1]:
            if boundary_conditions[b] == 'dirichlet':
                    # u(x,t) = 0
                A[-b,-b] = 1.0
                A[-b,1-3*b] = 0.0
                B[-b,-b] = 0.0
                B[-b,1-3*b] = 0.0
            elif boundary_conditions[b] == 'neumann':
                    # u'(x,t) = 0
                A[-b,1-3*b] = -2*sig
                B[-b,1-3*b] = 2*sig

            # Propagate
        psi = psi_init
        #print(B)
        #print(sig)
        inv_A=np.linalg.inv(A)
        fpsi_1=self.f(psi,0)
        
            
        for n in range(self.n_t):
            self.psi_matrix[n,:] = psi
            #fpsi = self.f(psi,n)
            fpsi_2=self.f(psi,n)
            
            psi = la.solve(A, B.dot(psi) - 1j*self.delta_t * (3/2*fpsi_2-1/2*fpsi_1))
            #psi=inv_A.dot(B.dot(psi) - 1j*self.delta_t * (3/2*fpsi_2-1/2*fpsi_1))
            fpsi_1=fpsi_2
            
            #psi = np.linalg.solve(A, B.dot(psi))
            
        '''
  
        A = self._fillA_sp(sig, self.n_x, data_type)
        B = self._fillB_sp(sig, self.n_x, data_type)
            
            # Set boundary conditions
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
                
        # Propagate
        psi = psi_init
        fpsi_1=self.f(psi,0)
        for n in range(self.n_t):
            self.psi_matrix[n,:] = psi
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
        
    '''
    
    def _make_tridiag(self, sig, n, data_type):
    
        M = np.diagflat(np.full(n, (1+2*sig), dtype=data_type)) + \
            np.diagflat(np.full(n-1, -(sig), dtype=data_type), 1) + \
            np.diagflat(np.full(n-1, -(sig), dtype=data_type), -1)

        return M
    '''  
    
    def _fillA_sp(self, sig, n, data_type):
        """Returns a tridiagonal matrix in compact form ab[1+i-j,j]=a[i,j]"""
        
        A = np.zeros([3,n], dtype=data_type) # A has three diagonals and size n
        A[0] = -sig # superdiagonal
        A[1] = 1+2*sig # diagonal
        A[2] = -sig # subdiagonal
        return A

    def _fillB_sp(self, sig, n, data_type):
        """Returns a tridiagonal sparse matrix in csr-form"""
        
        _o = np.ones(n, dtype=data_type)
        supdiag = (sig)*_o[:-1]
        diag = (1-2*sig)*_o
        subdiag = (sig)*_o[:-1]
        return scipy.sparse.diags([supdiag, diag, subdiag], [1,0,-1], (n,n), format="csr")