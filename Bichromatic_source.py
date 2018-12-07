#####################################################################################################################
###################################DESCRIPTION OF THE FILE###########################################################
#####################################################################################################################
# This file is used to fully study the case of a bichromatic laser source <-> E(t)=E0*sin(wt+phi1)+E2*sin(2wt+phi2)
# The trajectories will first be plotted for several ionization times between 0 and T.
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
# In the MONOCHROMATIC CASE we take E_2 = E0*sin(wt+phi1)+E2*sin(2wt+phi2) and Laser_source_2
#####################################################################################################################

LASER_SOURCE=Laser_source_2
E=E_2



""" An array of ionization times """
Number_of_t_prime=11
T_prime=np.linspace(0,T,Number_of_t_prime)

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
    plt.plot(trajectory(T_prime[i],LASER_SOURCE,3*T)[0],trajectory(T_prime[i],LASER_SOURCE,3*T)[1][:,0],color=rgb_to_hex(Echelle_c[i]),label="$t_i$={0}T".format(T_prime[i]),linewidth=3)
plt.plot(np.linspace(0,3*T,10),np.linspace(0,0,10),color='k')  
plt.grid() 
plt.xlabel('Time $t$', fontsize=30)
plt.xticks( np.arange(5)*T, (r"$0$", r"$T$", r"$2T$", r"$3T$", r"$4T$"), fontsize=30)
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
T_Prime_converge=np.linspace(0.25*T,0.5*T,Number_of_t_prime) #The electron will go back to the atom
T_Prime_diverge=np.linspace(0*T,0.25*T,Number_of_t_prime)    #The electron will never go back to the atom
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_Prime_diverge[i],LASER_SOURCE,3*T)[0],trajectory(T_Prime_diverge[i],LASER_SOURCE,3*T)[1][:,0],color=rgb_to_hex(Echelle_d[i]),linewidth=2)
    plt.plot(trajectory(T_Prime_converge[i],LASER_SOURCE,3*T)[0],trajectory(T_Prime_converge[i],LASER_SOURCE,3*T)[1][:,0],color=rgb_to_hex(Echelle_c[i]),linewidth=2)
    plt.scatter(T_Prime_diverge[i],E(T_Prime_diverge[i]),color=rgb_to_hex(Echelle_d[i]),linewidth=3)
    plt.scatter(T_Prime_converge[i],E(T_Prime_converge[i]),color=rgb_to_hex(Echelle_c[i]),linewidth=3)
s=np.linspace(0,3*T,Number_of_points)
plt.plot(s,np.linspace(0,0,Number_of_points),color='k')
plt.plot(s,E(s),color='k',linewidth=2)
plt.xlabel('Ionization Time $t_i$', fontsize=30)
plt.ylabel('Trajectories and Electric Field', fontsize=20)
plt.xticks( np.arange(4)*T, (r"$0$", ur"$T$", r"$2T$", r"$3T$"), fontsize=30)
plt.yticks( np.arange(3)-E0, (r"$-E_{0}$" , r"$0$", r"$E_{0}$"), fontsize=30)
plt.grid()
plt.tight_layout() 
"""


########################################################################################################################
# Visualisation of the trajectories for a time of ionization in [T/2,T], the moment when the electron leaves the atom is highlighted with a colored bullet
########################################################################################################################

#We are interested here only in the case where t_prime is in [T/2,T].
"""
plt.figure()
T_Prime_converge=np.linspace(0.75*T,1*T,Number_of_t_prime)     #The electron will go back to the atom
T_Prime_diverge=np.linspace(0.5*T,0.75*T,Number_of_t_prime)    #The electron will never go back to the atom
for i in range(Number_of_t_prime):
    plt.plot(trajectory(T_Prime_diverge[i],LASER_SOURCE,3*T)[0],trajectory(T_Prime_diverge[i],LASER_SOURCE,3*T)[1][:,0],color=rgb_to_hex(Echelle_d[i]),linewidth=2)
    plt.plot(trajectory(T_Prime_converge[i],LASER_SOURCE,3*T)[0],trajectory(T_Prime_converge[i],LASER_SOURCE,3*T)[1][:,0],color=rgb_to_hex(Echelle_c[i]),linewidth=2)
    plt.scatter(T_Prime_diverge[i],E(T_Prime_diverge[i]),color=rgb_to_hex(Echelle_d[i]),linewidth=3)
    plt.scatter(T_Prime_converge[i],E(T_Prime_converge[i]),color=rgb_to_hex(Echelle_c[i]),linewidth=3)
s=np.linspace(0,3*T,Number_of_points)
plt.plot(s,np.linspace(0,0,Number_of_points),color='k')
plt.plot(s,E(s),color='k',linewidth=2)
plt.xlabel('Ionization Time $t_i$', fontsize=30)
plt.ylabel('Trajectories and Electric Field', fontsize=20)
plt.xticks( np.arange(4)*T, (r"$0$", r"$T$", r"$2T$", r"$3T$"), fontsize=30)
plt.yticks( np.arange(3)-E0, (r"$-E_{0}$" , r"$0$", r"$E_{0}$"), fontsize=30)
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
T_prime=np.linspace(0,T,100)
List_of_good_ionization_times=[]
for i in range(len(T_prime)):
    indice=idx_change_of_sign(trajectory(T_prime[i],LASER_SOURCE,3*T)[1][:,0])
    if indice>0 :
        List_of_good_ionization_times.append(0.01*i*T)
print(List_of_good_ionization_times)
"""
#the list we get is : [0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.5700000000000001, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.81, 0.8200000000000001, 0.8300000000000001, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.9400000000000001, 0.9500000000000001, 0.96, 0.97, 0.98]


"""Initialization of vectors of ionization time using the precedent result"""

#We can choose to work either on the reunion of [0.2*T,0.32*T], [0.42*T,0.65*T] and [0.82*T, 0.98*T]

TPRIME=np.linspace(0.2*T,0.32*T,100)
TPRIME_bis=np.linspace(0.42*T,0.65*T,100)
TPRIME_ter=np.linspace(0.82*T,0.99*T,100)


##############################################################################################################################
# RETURN TIME AS A FUNCTION OF THE IONIZATION TIME. 
# It should be a strictly decreasing function in the case of a monochromatic source which suggests that there's only one return time for one ionization time. (bijection)
# The theoretical solution is compared to the numerical solution
##############################################################################################################################

"""

for i in range(len(TPRIME)):
    plt.scatter(TPRIME[i],find_Treturn(TPRIME[i],LASER_SOURCE,3*T),marker='+',color='red',lw=3)
    plt.scatter(TPRIME_bis[i],find_Treturn(TPRIME_bis[i],LASER_SOURCE,3*T),marker='+',color='red',lw=3)
    plt.scatter(TPRIME_ter[i],find_Treturn(TPRIME_ter[i],LASER_SOURCE,3*T),marker='+',color='red',lw=3)

plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Ionization Time $t_i$ ', fontsize=30)
plt.xticks( np.arange(5)*T*0.25, (r"$0$",r"$\frac{T}{4}$", r"$\frac{T}{2}$",r"$\frac{3T}{4}$",r"$T$"), fontsize=30)
plt.yticks( 0.25*T+np.arange(7)*(T)*0.25, (r"$\frac{T}{4}$", r"$\frac{T}{2}$",r"$\frac{3T}{4}$", r"$T$", r"$\frac{5T}{4}$",r"$\frac{3T}{2}$",r"$\frac{7T}{4}$" ), fontsize=30)
plt.ylabel('Return time $t_r$', fontsize=30)
plt.tight_layout()

"""
#################################################################################
# TIME OF FLIGHT AS A FUNCTION OF THE RETURN TIME
#Here is plotted the time of flight (t_return-t_ionization) as a function of the return time.
#################################################################################
"""
plt.figure()
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i],LASER_SOURCE,3*T),find_Treturn(TPRIME[i],LASER_SOURCE,3*T)-TPRIME[i])
    plt.scatter(find_Treturn(TPRIME_bis[i],LASER_SOURCE,3*T),find_Treturn(TPRIME_bis[i],LASER_SOURCE,3*T)-TPRIME_bis[i])
    plt.scatter(find_Treturn(TPRIME_ter[i],LASER_SOURCE,3*T),find_Treturn(TPRIME_ter[i],LASER_SOURCE,3*T)-TPRIME_ter[i])
plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Return Time $t_r$ (T)', fontsize=30)
plt.ylabel('Time of flight $\\tau$', fontsize=30)
plt.yticks( np.arange(5)*T*0.25, (r"$0$",r"$\frac{T}{4}$", r"$\frac{T}{2}$", r"$\frac{3T}{4}$", r"$T$"), fontsize=30)
plt.xticks( 0.25*T+np.arange(7)*(T)*0.25, (r"$\frac{T}{4}$", r"$\frac{T}{2}$",r"$\frac{3T}{4}$", r"$T$", r"$\frac{5T}{4}$",r"$\frac{3T}{2}$",r"$\frac{7T}{4}$" ), fontsize=30)
plt.tight_layout()

"""

#################################################################################
# KINETIC RETURN ENERGY AS A FUNCTION OF THE RETURN TIME
# Here is plotted the Kinetic Return Energy as a function of the Return Time.
##############################################################################
"""

plt.figure()
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i],LASER_SOURCE,3*T),16*(np.pi)**2*Kinetic_energy(TPRIME[i],LASER_SOURCE,3*T))
    plt.scatter(find_Treturn(TPRIME_bis[i],LASER_SOURCE,3*T),16*(np.pi)**2*Kinetic_energy(TPRIME_bis[i],LASER_SOURCE,3*T))
    plt.scatter(find_Treturn(TPRIME_ter[i],LASER_SOURCE,3*T),16*(np.pi)**2*Kinetic_energy(TPRIME_ter[i],LASER_SOURCE,3*T))
plt.grid() 
plt.yticks(fontsize=30)
plt.xlabel('Return Time $t_r$ (T)', fontsize=30)
plt.ylabel('Kinetic Return Energy $\\frac{E_k}{U_p}$', fontsize=25)
plt.xticks( 0.25*T+np.arange(7)*(T)*0.25, (r"$\frac{T}{4}$", r"$\frac{T}{2}$",r"$\frac{3T}{4}$", r"$T$", r"$\frac{5T}{4}$",r"$\frac{3T}{2}$",r"$\frac{7T}{4}$" ), fontsize=30)
plt.tight_layout()
"""





plt.legend()
plt.show()
