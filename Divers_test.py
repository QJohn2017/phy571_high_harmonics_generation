# -*- coding: utf-8 -*-

"""
Représentation temporelle d'une onde
"""

from __future__ import division
from scipy import *
from pylab import *
import scipy.optimize 
from matplotlib import animation
from Physical_quantities import *
#from Useful_functions import *



x = linspace(0,10,500)

# Animation

def y(x,t):
    return cos(2*pi*t-pi*x)				# Fonction à animer

omega=2*pi

def g(x):
    return sin(omega*x)

def h(x):
    return tan(x)-x

def d(x):
    return omega*cos(omega*t_prime)*x+sin(omega*t_prime)-omega*t_prime*cos(omega*t_prime)

S=scipy.optimize.brenth(h,2*pi/3,3*pi/2)


Tprime=linspace(S/(2*pi),S/(2*pi)+1/100,10)
X=linspace(0,1,100)
for i in range(len(Tprime)):
    t_prime=Tprime[i]
    plot(X,d(X))
    plot(X,g(X))

show()
"""
X=linspace(0,0.1,100)

plot(X,f(X))


fig1 = figure()
ax = fig1.add_subplot(111)

xlim(0, 10)                                   			# Limites de l'axe des abscisses
xlabel(ur"$x \, (m)$", fontsize=16)        	# Label de l'axe des abscisses

ylim(-1.5, 1.5)                                			# Limites de l'axe des ordonnées
ylabel(ur"$y \, (m)$", fontsize=16)        	# Label de l'axe des ordonnées

line, = plot([], [], '-', color='b', lw=2)          	# Mise en forme extérieure à la boucle pour effacement entre chaque image

def init():          						# Initialisation avec du vide
    line.set_data([], [])
    return line,

def animate(i):
    t=i/25
    line.set_data(x, y(x,t))
    print("frame num %d" %i)
    return line,						# Ne pas oublier la "," pour le tupple


anim = animation.FuncAnimation(fig1, animate, frames=26, interval=2, init_func=init)

anim.save('representation_temporelle_1D.mp4', fps=25)	# Sauvegarde du film

show()"""

#Plot the theretical solution for a monochromatic source

"""TPrime=np.linspace(0,T,15)

for i in range(len(TPrime)):
    t_prime=TPrime[i]
    Temps=np.linspace(t_prime,3*T,100)
    plt.plot(Temps,X_1_theorical(Temps))"""

def E_1(t):
    return E1*np.sin(omega*t+phi1)


def X_1_theoretical(t):
    Vd = q*E1/(m*omega)*np.cos(omega*t_prime+phi1)
    X0 = q*E1/(m*omega)*(1/omega*np.sin(omega*t_prime+phi1)-t_prime*np.cos(omega*t_prime+phi1))
    return -q*E_1(t)/(m*omega**2) + Vd*(t) + X0

def find_Treturn_theoretical(t):
    return scipy.optimize.brenth(X_1_theoretical,t_prime,T+t_prime)
    
TPRIME=np.linspace(0.25001,0.499,10)
for i in range(len(TPRIME)):
    t_prime=TPRIME[i]
    print(scipy.optimize.brenth(X_1_theoretical,t_prime,T+t_prime))
    #plt.scatter(TPRIME[i],find_Treturn_theoretical(TPRIME[i]))
  
show()
