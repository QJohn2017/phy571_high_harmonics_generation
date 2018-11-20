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


j=complex(0,1)

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
    def solve(self, psi_init, boundary_conditions=('neumann','neumann')):
            
        sig = j* self.delta_t / 4. / self.delta_x**2
        
        
        # Figure the data type
        data_type = type(sig*psi_init[0])
        
        self.psi_matrix = np.zeros([self.n_t, self.n_x], dtype=data_type)  #matrice des fonctions d'onde
  


   
            
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
        for n in range(self.n_t):
            self.psi_matrix[n,:] = psi
            fpsi = f(psi,n)
            if n==0: fpsi_old = fpsi
            psi = la.solve(A, B.dot(psi) + self.delta_t * (fpsi))
            fpsi_old = fpsi
  
            
    def get_final_psi(self):
        
        return self.psi_matrix[-1,:].copy()
    
    
    def get_first_psi(self):
       
        return self.psi_matrix[0,:].copy()
    
    def get_second_psi(self):
        return self.psi_matrix[1,:].copy()
        

    
    def _make_tridiag(self, sig, n, data_type):
    
        M = np.diagflat(np.full(n, (1+2*sig), dtype=data_type)) + \
            np.diagflat(np.full(n-1, -(sig), dtype=data_type), 1) + \
            np.diagflat(np.full(n-1, -(sig), dtype=data_type), -1)

        return M