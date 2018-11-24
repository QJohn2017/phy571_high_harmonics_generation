#This file is used to plot the evolution of the trajectory with respect to the time. The indidence of the ionization time t_prime on the result is studied.
#t_prime is the time when the electron leaves the atomical orbit. It's the beginning of the classical solution. We define here a vector of t_prime of length Number_of_t_prime

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Physical_quantities import * 
from Usefull_functions import * 
from Integration import * 
from color import *

#Number of points considered

Number_of_points=100
Number_of_t_prime=11


"""Utilisation of the echelle function implemented in color.py"""
Echelle_c = echelle(jaune,rouge,Number_of_t_prime)  #gradient of color from yellow to red
Echelle_d = echelle(bleu,jaune,Number_of_t_prime)   #gradient of color from blue to yellow


#We define the function trajectory which will convey the position and velocity of the electron at each time between t_prime and t_prime+2*T_prime
"""The first argument of the ouput is a table of time starting at t_prime and ending at t_prime+2*T. Dim of the matrix = (1,N)
The second argument of the output are the position and velocity for every time. Dim of the matrix = (N,2)  """
def trajectory(t_prime):
    y0=[0,0]
    TRAJ=[np.zeros(Number_of_points),np.zeros([2,Number_of_points])]
    T_PRIME = np.linspace(t_prime,3*T,Number_of_points)
    SOL = Integration(Laser_source_1,y0,T_PRIME)
    TRAJ[0]=T_PRIME
    TRAJ[1]=SOL
    return TRAJ


#We plot here the trajectories for an ionization time in [0,T]
T_prime=np.linspace(0,T,Number_of_t_prime)
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_prime[i])[0],trajectory(T_prime[i])[1][:,0],color=rgb_to_hex(Echelle_c[i]))
plt.plot(np.linspace(0,3*T,10),np.linspace(0,0,10),color='k')  
plt.grid() 
plt.show()


#In reality there are to cases which need to be distinguished.
#On the one hand, the case where the electron never goes back to the atom : It occurs for an ionization time in [0,T/4] or in [T/2,3T/4]
#On the other hand, the case where the electron goes back to the atom : It occurs for an ionization time in [T/4,T/2] or in [3T/4,T]. We can then define a return time t_return
#This can clearly be seen on the precedent plot. When the trajectory crosses the horizontal axe, the electron comes back to the atom. 

#We are interested here only in the case where t_prime is in [0,T/2].
plt.figure()
T_Prime_converge=np.linspace(0.25*T,0.5*T,Number_of_t_prime)
T_Prime_diverge=np.linspace(0*T,0.25*T,Number_of_t_prime)
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_Prime_diverge[i])[0],trajectory(T_Prime_diverge[i])[1][:,0],color=rgb_to_hex(Echelle_d[i]))
    plt.plot(trajectory(T_Prime_converge[i])[0],trajectory(T_Prime_converge[i])[1][:,0],color=rgb_to_hex(Echelle_c[i]))
    plt.scatter(T_Prime_diverge[i],E_1(T_Prime_diverge[i]),color=rgb_to_hex(Echelle_d[i]))
    plt.scatter(T_Prime_converge[i],E_1(T_Prime_converge[i]),color=rgb_to_hex(Echelle_c[i]))
s=np.linspace(0,3*T,Number_of_points)
plt.plot(s,np.linspace(0,0,Number_of_points),color='k')
plt.plot(s,E_1(s),color='k')
plt.grid() 
plt.show()

#We are interested here only in the case where t_prime is in [T/2,T].
plt.figure()
T_Prime_converge=np.linspace(0.75*T,1*T,Number_of_t_prime)
T_Prime_diverge=np.linspace(0.5*T,0.75*T,Number_of_t_prime)
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_Prime_diverge[i])[0],trajectory(T_Prime_diverge[i])[1][:,0],color=rgb_to_hex(Echelle_d[i]))
    plt.plot(trajectory(T_Prime_converge[i])[0],trajectory(T_Prime_converge[i])[1][:,0],color=rgb_to_hex(Echelle_c[i]))
    plt.scatter(T_Prime_diverge[i],E_1(T_Prime_diverge[i]),color=rgb_to_hex(Echelle_d[i]))
    plt.scatter(T_Prime_converge[i],E_1(T_Prime_converge[i]),color=rgb_to_hex(Echelle_c[i]))
s=np.linspace(0,3*T,Number_of_points)
plt.plot(s,np.linspace(0,0,Number_of_points),color='k')
plt.plot(s,E_1(s),color='k')
plt.grid()
plt.show()






