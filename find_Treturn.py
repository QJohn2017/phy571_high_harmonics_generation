#This file is used to calculate the time of return when it exists using the function trajectory calculated in an other file.
#t_return must be greater than t_ionization.
#We then plot t_return as a function of t_prime.
#The time of flight and the E_kin can also be plotted as functions of the return time.

import numpy as np
import numpy.random as rnd
import numpy.linalg
import matplotlib.pyplot as plt
import scipy.optimize 
from Integration import *
from trajectory_evolution import *


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


def find_Treturn(t_prime):
    #Find the return time if it exists. 
    #In fact, we'll chosse t_return such as the trajectory is the nearest from 0. It could either be the time before the change or the time after.
    indice=idx_change_of_sign(trajectory(t_prime)[1][:,0])
    t_after_crossing=trajectory(t_prime)[0][indice]
    t_before_crossing=trajectory(t_prime)[0][indice-1]
    #Chose which t is correct
    if np.abs(trajectory(t_prime)[1][indice,0])<np.abs(trajectory(t_prime)[1][indice-1,0]):
        return t_after_crossing
    else :
        return t_before_crossing


def Kinetic_energy(t_prime):
    #Calcultate the kinetic energy 
    indice=idx_change_of_sign(trajectory(t_prime)[1][:,0])
    V2_max=max((trajectory(t_prime)[1][indice,1])**2,(trajectory(t_prime)[1][indice-1,0])**2)
    return 0.5*m*V2_max


TPRIME=np.linspace(0.25001,0.499,1000)
TPRIME_bis=np.linspace(0.75001,0.999,1000)
# We plot here the return time as a function of the ionization time. It should be a strictly decreasing
# function in the case of a monochromatic source which
# suggests that there's only one return time for one ionization time.
#TEST for t' in [0.25,0.5]
"""
for i in range(len(TPRIME)):
    plt.scatter(TPRIME[i],find_Treturn(TPRIME[i]))
    plt.scatter(TPRIME_bis[i],find_Treturn(TPRIME_bis[i]))"""


# Here is plotted the time of flight (t_return-t_ionization) as a function of the return time.
#TEST for t' in [0.25,0.5]
"""
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i]),find_Treturn(TPRIME[i])-TPRIME[i])
    plt.scatter(find_Treturn(TPRIME_bis[i]),find_Treturn(TPRIME_bis[i])-TPRIME_bis[i])"""


# Here is plotted the time of flight (t_return-t_ionization) as a function of the return time.
#TEST for t' in [0.25,0.5]

for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i]),16*(np.pi)**2*Kinetic_energy(TPRIME[i]))
    #plt.scatter(find_Treturn(TPRIME_bis[i]),find_Treturn(TPRIME_bis[i])-TPRIME_bis[i])"""




plt.show()

