
# coding: utf-8

# In[43]:



import numpy as np
import numpy.random as rnd
import numpy.linalg
import matplotlib.pyplot as plt
import scipy.optimize 


# ## DESCRIPTION DU FICHIER / CONSTANTES PHYSIQUES

# In[185]:


#Physical Constants
q=1
E0=1
m=1
omega=2*np.pi
T=1
alpha=q*E0/(m*omega)
  


# ## CALCUL DES FONCTIONS UTILES

# In[186]:


#Solution theorique
def X(t):
    return alpha*(1/omega*(np.sin(omega*t_prime)-np.sin(omega*t))+np.cos(omega*t_prime)*(t-t_prime))

#Vitesse theorique
def V(t):
    return alpha*(np.cos(omega*t_prime)-np.cos(omega*t))

#Energie cinetique theorique
def Ek(t):
    return 0.5*m*V(t)**2

#Champ Electrique
def E(t):
    return E0*np.sin(omega*t)



# ## TRACE DU CHAMP ET DU MOUVEMENT DES ELECTRONS

# In[187]:


#Traces du champ et du mouvements des electrons
#On peut reperer les temps d'ionisation pour lesquels l'electron s'echappe
t=np.linspace(0,2*T,100)
T_prime=np.linspace(0.28*T,0.32*T,2)

plt.figure()
for i in range(len(T_prime)):
    t_prime=T_prime[i]
    plt.plot(t+t_prime,X(t+t_prime),label='t_prime = %s T' %(t_prime))
    
    
plt.plot(t,0.1*E(t),color='black',label='Champ E')
plt.plot(1.5*t,np.zeros(100),color='black')
plt.legend()
plt.show()

# ## CALCUL DES TEMPS DE RETOUR 

# In[206]:


# Les seuls temps d'ionization interessants sont contenus dans la reunin d'intervalle 
plt.figure()
TPRIME=[]#Liste contenant les temps d'ionization
TSECONDE=[]
Energy=[]
N=10000
for i in range(N):
    #On se place sur l'intervalle [T/4,T/2]
    t_prime=0.25*T+0.25*T*rnd.random()
    t_seconde=scipy.optimize.brenth(X,0,10000*T)
    Energy.append(Ek(t_seconde))
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



# ## CALCUL DE L'ENERGIE CINETIQUE DE RETOUR
# L'energie cinetique de retour est obtenue en calculant la valeur de la fonction Ek en fonction de T_ionization en prenant t=T_RETOUR

# In[193]:


plt.figure()
Energy=[e*16*(np.pi)**2 for e in Energy] #Normalisation de l'energie Cinetique
plt.scatter(TPRIME,Energy,marker='o')
plt.xlabel('T_IONIZATION')  
plt.ylabel('Ek/Up')
plt.show() 
M=max(Energy)
M

