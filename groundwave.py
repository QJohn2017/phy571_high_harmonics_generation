# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:42:21 2018

@author: marti
"""

import numpy as np

class Groundstate:
    """This class allows you to calculate the groundstate wavefunction of a particle inside a potential V"""
    def __init__(self, N, dt, xmin, xmax, V):    
        self.N = N 
        self.dx = (xmax-xmin)/N """spatial resolution"""
        self.dt = dt """temporal resolution"""
        self.xmin = xmin """inferior spatial boundary"""
        self.xmax = xmax """superior spatial boundary"""
        self.V = V
    
    