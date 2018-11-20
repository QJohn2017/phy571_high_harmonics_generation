#This file contains all the functions which will be usefull for the project

import numpy as np
import numpy.random as rnd
import numpy.linalg
import matplotlib.pyplot as plt
from Physical_quantities import * 

#We present here the possible function used for the laser source : Monochromatic, 2*omega and 3*omega

"""The Simplest Laser source (monochromatic)"""
def E_1(t):
    return E1*np.sin(omega*t+phi1)

"""2 Frequencies"""
def E_2(t):
    return E1*np.sin(omega*t+phi1)+E2*np.sin(2*omega*t+phi2)

"""3 Frequencies"""
def E_3(t):
    return E1*np.sin(omega*t+phi1)+E2*np.sin(2*omega*t+phi2)+E3*np.sin(3*omega*t+phi3)

#Here are the Theoretical expressions for the movement when the source is monochromatic. t_prime is the ionization time.

def X_1_theorical(t):
    Vd = q*E1/(m*omega)*np.cos(omega*t_prime+phi1)
    X0 = q*E1/(m*omega)*(1/omega*np.sin(omega*t_prime+phi1)-np.cos(omega*t_prime+phi1))
    return -q*E_1(t)/(m*omega**2) + Vd*t + X0

def V_1_theorical(t):
    Vd = q*E1/(m*omega)*np.cos(omega*t_prime+phi1)
    return -q*E1/(m*omega)*np.cos(omega*t)+Vd

def Ekin_1_theorical(t):
    return 0.5*m*(V_1_theorical(t))**2

