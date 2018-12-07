#This file contains all the implemented functions which will be useful for the project

import numpy as np
import numpy.random as rnd
import numpy.linalg
import matplotlib.pyplot as plt
from Physical_quantities import * 
from scipy.integrate import odeint
from scipy.optimize import bisect
from color import *

################################################################################################################################
#We present here the possible function used for the laser source : Monochromatic, 2*omega and 3*omega and a more realistic pulse with an enveloppe
################################################################################################################################

"""The Simplest Laser source (monochromatic)"""
def E_1(t):
    return E0*np.sin(omega*t+phi1)

"""2 Frequencies"""
def E_2(t):
    return E0*np.sin(omega*t+phi1)+E2*np.sin(2*omega*t+phi2)

"""3 Frequencies"""
def E_3(t):
    return E0*np.sin(omega*t+phi1)+E2*np.sin(2*omega*t+phi2)+E3*np.sin(3*omega*t+phi3)

"""Enveloppe"""
def E_env(t):
    return E0_env*(np.sin(np.pi*t/(tau))**2)*np.sin(omega_env*t)

#######################################################################################################################
#Here are the Theoretical expressions for the movement when the source is monochromatic. t_prime is the ionization time.
########################################################################################################################
def X_1_theoretical(t,t_prime):
    Vd = q*E0/(m*omega)*np.cos(omega*t_prime+phi1)
    X0 = q*E0/(m*omega)*(1/omega*np.sin(omega*t_prime+phi1)-t_prime*np.cos(omega*t_prime+phi1))
    return -q*E_1(t)/(m*omega**2) + Vd*(t) + X0

#X_1_theorical=np.vectorize(X_1_theorical)

def V_1_theoretical(t,t_prime):
    Vd = q*E0/(m*omega)*np.cos(omega*t_prime+phi1)
    return -q*E0/(m*omega)*np.cos(omega*t)+Vd

def Ekin_1_theoretical(t,t_prime):
    return 0.5*m*(V_1_theoretical(t,t_prime))**2


#Calculate the theoretical return time as a function of t_ionization

def find_Treturn_theoretical(t_prime):
    L=10000
    Temps=np.linspace(t_prime,2*T,L)
    X=X_1_theoretical(Temps,t_prime)
    idx=idx_change_of_sign(X)
    t_before_crossing=Temps[idx-1]
    t_after_crossing=Temps[idx]
    if np.abs(X[idx])<np.abs(X[idx-1]):
        return t_after_crossing
    else :
        return t_before_crossing

find_Treturn_theoretical=np.vectorize(find_Treturn_theoretical)


#Calculate the theoretical kinetic energy

def Kinetic_energy_theoretical(t_prime):
     return 16*(np.pi)**2*Ekin_1_theoretical(find_Treturn_theoretical(t_prime),t_prime)

Kinetic_energy_theoretical=np.vectorize(Kinetic_energy_theoretical)




################################################################################################################################
####### We integrate the movement equation for any laser source represented by a function E(t) using the odeint method. 
################################################################################################################################


"""Definition of the vector which will be integrated by the odeint function depending on the electric field use"""
def Laser_source_1(y,t):
    x, x_dot = y
    dydt = [x_dot, q*E_1(t)/m]
    return dydt

def Laser_source_2(y,t):
    x, x_dot = y
    dydt = [x_dot, q*E_2(t)/m]
    return dydt

def Laser_source_3(y,t):
    x, x_dot = y
    dydt = [x_dot, q*E_3(t)/m]
    return dydt

def Laser_source_env(y,t):
    x, x_dot = y
    dydt = [x_dot, q*E_env(t)/m]
    return dydt

def Integration(LASER_SOURCE,y0,time_array):
    return odeint(LASER_SOURCE,y0,time_array)


"""TEST des fonctions
def x(t):
    return t**2

print(Laser_source_1([0,0],2))
print(Integration(Laser_source_1, [0,0],[0,1,2])[0])"""


################################################################################################################################
#We define the function trajectory which will convey the position and velocity of the electron at each time between t_prime and 3*T_prime. The inputs are a time of ionization and the LASER_SOURCE function.
################################################################################################################################

"""The first argument of the ouput is a table of time starting at t_prime and ending at 3*T. Dim of the matrix = (1,N)
The second argument of the output are the position and velocity for every time. Dim of the matrix = (N,2)  """
def trajectory(t_prime,LASER_SOURCE,periode):
    y0=[0,0]
    TRAJ=[np.zeros(Number_of_points),np.zeros([2,Number_of_points])]
    T_PRIME = np.linspace(t_prime,periode,Number_of_points)
    #Choice of the laser source
    SOL = Integration(LASER_SOURCE,y0,T_PRIME)   
    TRAJ[0]=T_PRIME
    TRAJ[1]=SOL
    return TRAJ



################################################################################################################################
# The function idx_change_of_sign will give the first index where the an array change of sign. It will be tested for arrays where the first term is always equal to zero (x(0)=0). Therefore, the change is saught starting at the second element.
################################################################################################################################

def idx_change_of_sign(Array):
    #We try to find out if the trajectory crosses the x axis (i.e if the sign changes)
    N=len(Array)
    T=[]
    for i in range(2,N):
        #We start at 2 because the first term will always be 0 (the initial position is x(0)=0)
        if Array[i]*Array[i-1]<0:
            T.append(i)
    if T!=[]:
        return (T[0])
    else :
        print('!!!!!!!! The electron never returns to the atom !!!!!!!!')
        return (-1)


################################################################################################################################
# The function idx_change_of_sign will give the first index where the an array change of sign. It will be tested for arrays where the first term is always equal to zero (x(0)=0). Therefore, the change is saught starting at the second element.
################################################################################################################################

def find_Treturn(t_prime,LASER_SOURCE,periode):
    #Find the return time if it exists. 
    #In fact, we'll chosse t_return such as the trajectory is the nearest from 0. It could either be the time before the change or the time after.
    indice=idx_change_of_sign(trajectory(t_prime,LASER_SOURCE,periode)[1][:,0]) #The function is built in the Useful_functions.py file
    t_after_crossing=trajectory(t_prime,LASER_SOURCE,periode)[0][indice]
    t_before_crossing=trajectory(t_prime,LASER_SOURCE,periode)[0][indice-1]
    #Chose which t is correct (the one for which the value of the trajectory is the minimum)
    if np.abs(trajectory(t_prime,LASER_SOURCE,periode)[1][indice,0])<np.abs(trajectory(t_prime,LASER_SOURCE,periode)[1][indice-1,0]):
        return t_after_crossing
    else :
        return t_before_crossing


def Kinetic_energy(t_prime,LASER_SOURCE,periode):
    #Calcultate the kinetic energy 
    indice=idx_change_of_sign(trajectory(t_prime,LASER_SOURCE,periode)[1][:,0])
    V2_max=max((trajectory(t_prime,LASER_SOURCE,periode)[1][indice,1])**2,(trajectory(t_prime,LASER_SOURCE,periode)[1][indice-1,0])**2)
    return 0.5*m*V2_max

Kinetic_energy=np.vectorize(Kinetic_energy)







