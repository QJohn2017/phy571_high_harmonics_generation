#This file is used to plot the evolution of the trajectory with respect to the time. The indidence of the inization time t_prime on the result is studied.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Physical_quantities import * 
from Usefull_functions import * 
from Integration import * 
from color import *

Number_of_points=1000

"""t_prime is the time when the electron leaves the atomical orbit. It's the beginning of the classical solution. We define here a vector of t_prime oh length N_tprime."""
N_tprime=6
T_Prime_converge=np.linspace(0.25*T,0.5*T,N_tprime)
T_Prime_diverge=np.linspace(0*T,0.25*T,N_tprime)

"""Utilisation of the echelle function implemented in color.py"""
Echelle_c = echelle(jaune,rouge,N_tprime)  #gradient of color from yellow to red
Echelle_d = echelle(bleu,jaune,N_tprime)   #gradient of color from blue to yellow


def trajectory(t_prime):
    y0=[0,0]
    TRAJ=np.zeros([2,Number_of_points])
    T_PRIME = np.linspace(t_prime,t_prime+2*T,Number_of_points)
    SOL = Integration(Laser_source_1,y0,T_PRIME)
    TRAJ[0]=T_PRIME
    TRAJ[1]=SOL
    return TRAJ

Solution=trajectory(0.25*T)

print(Solution[0])


"""Exemple when t_prime is taken between 0 and T/2. 
for i in range(N_tprime-1):
  
    y0 = [0, 0]                            #Initial Conditions

    t_prime_c=T_Prime_converge[i]          #List of t_prime between 0
    t_prime_d=T_Prime_diverge[i]
    
    N=1000
    t_c = np.linspace(t_prime_c, t_prime_c+2*T, N)
    solE_c = odeint(Simplified_Laser_source, y0, t_c, args=(q,m,E1,omega,phi1))
    t_d = np.linspace(t_prime_d, t_prime_d+2*T, N)
    solE_d = odeint(Simplified_Laser_source, y0, t_d, args=(q,m,E1,omega,phi1))
    plt.plot(t_c,solE_c[:,0],color=rgb_to_hex(Echelle_c[i]))#,label='T_ionization = %0.3f'%(t_prime_c))
    plt.plot(t_d,solE_d[:,0],color=rgb_to_hex(Echelle_d[i]))#,label='T_ionization = %0.3f'%(t_prime_d)
    #plt.legend()
plt.plot(np.linspace(0,2.5*T,10),np.linspace(0,0,10),color='k')
plt.grid()
plt.show()"""
