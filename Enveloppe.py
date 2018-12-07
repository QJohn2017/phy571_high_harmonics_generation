#####################################################################################################################
###################################DESCRIPTION OF THE FILE###########################################################
#####################################################################################################################
# This file is used to fully study the case of a laser pulse <-> E(t)=E0*sin(pi*t/tau)^2*sin(omega*t)
# The trajectories will first be plotted for several ionization times between 0 and tau.
# These trajectories help us to understand that the electron doesn't always come back to the atom. It only does for some value of the ionization time. 
# The ionization times for which the electron goes back to the atom will be found empirically using the function idx_change_of_sign.
# We then find the return time and study the time of flight and the return kinetic energy as functions of t_return.
# Every function is under quotes, you need to unquote one by one and run the file litle by litle in order to get the good results
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################


from Useful_functions import * 


#####################################################################################################################
# The Electric field and the laser source needs to be chosen here.
# In the MONOCHROMATIC CASE we take E_3 = E0*sin(wt+phi1)+E2*sin(2wt+phi2)+E3*sin(3wt+phi3) and Laser_source_3
#####################################################################################################################

LASER_SOURCE=Laser_source_env
E=E_env



""" An array of ionization times """
Number_of_t_prime=11
T_prime=np.linspace(0,tau,Number_of_t_prime)

####################################################################################################################
# Use of the echelle function implemented in color.py
####################################################################################################################
Echelle_c = echelle(jaune,rouge,Number_of_t_prime)  #gradient of color from yellow to red
Echelle_d = echelle(bleu,jaune,Number_of_t_prime)   #gradient of color from blue to yellow





####################################################################################################################
# Visualisation of the trajectories for a time of ionization in [0,T]
####################################################################################################################

#We plot here the trajectories for an ionization time in [0,T]
"""
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_prime[i],LASER_SOURCE,3*tau)[0],trajectory(T_prime[i],LASER_SOURCE,3*tau)[1][:,0],color=rgb_to_hex(Echelle_c[i]),label='t_i=%0.0f' %T_prime[i],linewidth=3)
plt.plot(np.linspace(0,tau,10),np.linspace(0,0,10),color='k')  
plt.grid() 
plt.xlabel('Time $t$ (a.u)', fontsize=30)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.ylabel('Trajectories', fontsize=30)
plt.tight_layout()
"""


#Like in the monochromatic case, the electron doesn't always come back to the atom but the times of ionization for which it does are more complicated to find with a bichromatic source. We will address this issue later.
#For now, we plot the trajectories of the electron for an ionization time in  [0,T/2] or in [T/2,T].


########################################################################################################################
# Visualisation of the trajectories for a time of ionization in [0,T/2], the moment when the electron leaves the atom is highlighted with a colored bullet.
########################################################################################################################

#We are interested here only in the case where t_prime is in [0,T/2].
"""
plt.figure()
T_Prime_converge=np.linspace(0.25*tau,0.5*tau,Number_of_t_prime) #The electron will go back to the atom
T_Prime_diverge=np.linspace(0*tau,0.25*tau,Number_of_t_prime)    #The electron will never go back to the atom
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_Prime_diverge[i],LASER_SOURCE,tau)[0],trajectory(T_Prime_diverge[i],LASER_SOURCE,tau)[1][:,0],color=rgb_to_hex(Echelle_d[i]),linewidth=2)
    plt.plot(trajectory(T_Prime_converge[i],LASER_SOURCE,tau)[0],trajectory(T_Prime_converge[i],LASER_SOURCE,tau)[1][:,0],color=rgb_to_hex(Echelle_c[i]),linewidth=2)
    plt.scatter(T_Prime_diverge[i],E(T_Prime_diverge[i])/E0_env,color=rgb_to_hex(Echelle_d[i]),linewidth=3)
    plt.scatter(T_Prime_converge[i],E(T_Prime_converge[i])/E0_env,color=rgb_to_hex(Echelle_c[i]),linewidth=3)
s=np.linspace(0,1*tau,Number_of_points)
plt.plot(s,np.linspace(0,0,Number_of_points),color='k')
plt.plot(s,E_env(s)/E0_env,color='k',linewidth=2)
plt.xlabel('Ionization Time $t_i$ (u.a)' , fontsize=30)
plt.ylabel('Trajectories and Electric Field', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks( fontsize=20)
plt.grid()
plt.tight_layout() 
"""

########################################################################################################################
# Visualisation of the trajectories for a time of ionization in [T/2,T], the moment when the electron leaves the atom is highlighted with a colored bullet
########################################################################################################################

#We are interested here only in the case where t_prime is in [T/2,T].
"""
plt.figure()
T_Prime_converge=np.linspace(0.75*tau,1*tau,Number_of_t_prime)     
T_Prime_diverge=np.linspace(0.5*tau,0.75*tau,Number_of_t_prime)    

for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_Prime_diverge[i],LASER_SOURCE,1.5*tau)[0],trajectory(T_Prime_diverge[i],LASER_SOURCE,1.5*tau)[1][:,0],color=rgb_to_hex(Echelle_d[i]),linewidth=2)
    plt.plot(trajectory(T_Prime_converge[i],LASER_SOURCE,1.5*tau)[0],trajectory(T_Prime_converge[i],LASER_SOURCE,1.5*tau)[1][:,0],color=rgb_to_hex(Echelle_c[i]),linewidth=2)
    plt.scatter(T_Prime_diverge[i],E(T_Prime_diverge[i])/E0_env,color=rgb_to_hex(Echelle_d[i]),linewidth=3)
    plt.scatter(T_Prime_converge[i],E(T_Prime_converge[i])/E0_env,color=rgb_to_hex(Echelle_c[i]),linewidth=3)

s=np.linspace(0,1*tau,Number_of_points)
plt.plot(s,np.linspace(0,0,Number_of_points),color='k')
plt.plot(s,E_env(s)/E0_env,color='k',linewidth=2)
plt.xlabel('Ionization Time $t_i$ (u.a)' , fontsize=30)
plt.ylabel('Trajectories and Electric Field', fontsize=20)
plt.xticks( fontsize=20)
plt.yticks( fontsize=20)
plt.grid()
plt.tight_layout()
"""

###############################################################################################################################
# What follows is used to calculate the time of return when it exists using the function trajectory calculated in the Useful_functions.py file
#t_return must be greater than t_ionization. 
#We then plot t_return as a function of t_prime.
#The time of flight and the E_kin can also be plotted as functions of the return time.
##############################################################################################################################


#############################################################################################
# CHOICES OF THE INTERVALLS FOR THE IONIZATION TIME
# We first use the function idx_change_of sign to determine empirically the ionization times which give a returning electron.
# the output is a list of ionization times for which the electron comes back to the atom
#############################################################################################
"""
T_prime=np.linspace(0,tau,100)
List_of_good_ionization_times=[]
for i in range(len(T_prime)):
    indice=idx_change_of_sign(trajectory(T_prime[i],LASER_SOURCE,3*tau)[1][:,0])
    if indice>0 :
        List_of_good_ionization_times.append(0.01*i*T)
print(List_of_good_ionization_times)
"""
#the list we get is : [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.19, 0.2, 0.21, 0.22, 0.23, 0.27, 0.28, 0.29, 0.3, 0.31, 0.35000000000000003, 0.36, 0.37, 0.38, 0.43, 0.44, 0.45, 0.46, 0.51, 0.52, 0.53, 0.54, 0.59, 0.6, 0.61, 0.62, 0.67, 0.68, 0.6900000000000001, 0.75, 0.76, 0.77, 0.8300000000000001, 0.84, 0.85, 0.92, 0.93, 0.9400000000000001, 0.9500000000000001, 0.96, 0.97, 0.98, 0.99]




"""Initialization of vectors of ionization time using the precedent result"""

#We can choose to work either on this array taking out the bad numbers.

TPRIME=tau*np.array([ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,0.15, 0.19, 0.2, 0.21, 0.22, 0.23,0.27, 0.28, 0.29, 0.3, 0.31,0.36, 0.37,0.38, 0.44, 0.45, 0.46,  0.52, 0.53, 0.54,  0.6, 0.61, 0.67, 0.68,0.69, 0.75, 0.76, 0.77, 0.83, 0.84, 0.85, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98,0.99])


##############################################################################################################################
# RETURN TIME AS A FUNCTION OF THE IONIZATION TIME. 
# It should be a strictly decreasing function in the case of a monochromatic source which suggests that there's only one return time for one ionization time. (bijection)
# The theoretical solution is compared to the numerical solution
##############################################################################################################################
"""

for i in range(len(TPRIME)):
    plt.scatter(TPRIME[i],find_Treturn(TPRIME[i],LASER_SOURCE,3*tau),marker='+',lw=3)
   

plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Ionization Time $t_i$ (u.a)', fontsize=30)
plt.ylabel('Return time $t_r$ (u.a)', fontsize=30)
plt.xticks( fontsize=20)
plt.yticks( fontsize=20)
plt.tight_layout()
"""

#################################################################################
# TIME OF FLIGHT AS A FUNCTION OF THE RETURN TIME
#Here is plotted the time of flight (t_return-t_ionization) as a function of the return time.
#################################################################################
"""
plt.figure()
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i],LASER_SOURCE,3*tau),find_Treturn(TPRIME[i],LASER_SOURCE,3*tau)-TPRIME[i])

plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Return Time $t_r$ (u.a)', fontsize=30)
plt.ylabel('Time of flight $\\tau$ (u.a)', fontsize=30)
plt.xticks( fontsize=20)
plt.yticks( fontsize=20)
plt.tight_layout()

"""

#################################################################################
# KINETIC RETURN ENERGY AS A FUNCTION OF THE RETURN TIME
# Here is plotted the Kinetic Return Energy as a function of the Return Time.
##############################################################################


plt.figure()
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i],LASER_SOURCE,3*tau),1000*Kinetic_energy(TPRIME[i],LASER_SOURCE,3*tau))
   
plt.grid() 
plt.xlabel('Return Time $t_r$ (u.a)', fontsize=30)
plt.ylabel('Kinetic Return Energy $(1000 \\times E_k)$', fontsize=25)
plt.xticks( fontsize=20)
plt.yticks( fontsize=20)
plt.tight_layout()






plt.legend()
plt.show()
