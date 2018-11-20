import numpy as np
import numpy.random as rnd
import numpy.linalg
import matplotlib.pyplot as plt
import scipy.optimize 
from Integration import *

## Theoretical calculation of the return times

# In[206]:
alpha=q*E1/(m*omega)

#Theoretical Position 
def X(t):
    return alpha*(1/omega*(np.sin(omega*t_prime)-np.sin(omega*t))+np.cos(omega*t_prime)*(t-t_prime))

#Theoretical speed
def V(t):
    return alpha*(np.cos(omega*t_prime)-np.cos(omega*t))

#Kinetic Energy
def Ek(t):
    return 0.5*m*V(t)**2

# Les seuls temps d'ionization interessants sont contenus dans la reunin d'intervalle 
TPRIME=[]#Liste contenant les temps d'ionization
TSECONDE=[]
Energy=[]
N=10000
for i in range(N):
    #On se place sur l'intervalle [T/4,T/2]
    t_prime=0.25*T+0.25*T*rnd.random()
    t_seconde=scipy.optimize.brenth(X,0,10000*T)
    #Energy.append(Ek(t_seconde))
    TPRIME.append(t_prime)
    TSECONDE.append(t_seconde)

#for j in range(N):
    #On se place sur l'intervalle [3T/4,T]
    #t_prime=0.75+0.25*rnd.random()
    #t_seconde=scipy.optimize.brenth(X,0.01,1000*T)
    #Energy.append(Ek(t_seconde))
    #TPRIME.append(t_prime)
    #TSECONDE.append(t_seconde)


plt.figure()
TPRIME_=sorted(TPRIME)
TSECONDE_=sorted(TSECONDE,reverse=True)
plt.plot(TPRIME_,(np.log(TSECONDE_)))
plt.xlabel('T_IONIZATION')  
plt.ylabel('T_RETURN') 
plt.show()


def find_Treturn(function):
    return scipy.optimize.brenth(function,0,10000*T)
    

t_prime=0.26*T
print(find_Treturn(X))
