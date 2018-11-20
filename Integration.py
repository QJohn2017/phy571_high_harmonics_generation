#This document will be used to derive the solution from the movement equation
#It will use the python function : odeint
# The system needs to be written under the vectorial notation y=[x,x'] where one is looking for y'=[x',x''].


import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Physical_quantities import * 
from Usefull_functions import *
from color import *

#Here is the exemple of the pendulum where one can use of the function odeint :

"""def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt

b = 0.25
c = 5.0

y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 10, 101)


#TEST de la fonction
sol = odeint(pend, y0, t, args=(b, c))
plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()"""


####### We apply the same method but to the Simple Man's model where the laser source is represented by E(t)=E1*sin(omega*t+phi).

"""Definition of the vector which will be integrated by the odeint function in the monochromatic case."""
""""def Simplified_Laser_source(y, t, q, m, E1, omega, phi1):
    x, x_dot = y
    dydt = [x_dot, q*E1*np.sin(omega*t+phi1)/m]
    return dydt"""


####### We integrate the movement equation for any laser source represented by a function E(t).

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

def Integration(Source_function,y0,time_array):
    return odeint(Source_function,y0,time_array)


"""TEST des fonctions
def x(t):
    return t**2

print(Laser_source_1([0,0],2))
print(Integration(Laser_source_1, [0,0],[0,1,2])[0])"""


