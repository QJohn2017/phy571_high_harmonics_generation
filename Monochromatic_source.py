#####################################################################################################################
###################################DESCRIPTION OF THE FILE###########################################################
#####################################################################################################################
# This file is used to fully study the case of a monochromatic laser source <-> E(t)=E0*sin(wt+phi1)
# The trajectories will first be plotted for several ionization times between 0 and T.
# We'll focus on the [0,T/2] or [T/2,T] intervalls because the electric field is T/2 symmetric.
# These trajectories help us to understand that the electron doesn't always come back to the atom. It only does for some value of the ionization time.
# t_prime will be used as the time when the electron leaves the atomical orbit = ionization time.
# We then find the return time and study the time of flight and the return kinetic energy as functions of t_return
# Every function is under quotes, you need to unquote one by one and run the file litle by litle in order to get the good results
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################


from Useful_functions import * 


#####################################################################################################################
# The Electric field and the laser source needs to be chosen here.
# In the MONOCHROMATIC CASE we take E_1 = E0*sin(wt+phi1) and Laser_source_1
#####################################################################################################################

LASER_SOURCE=Laser_source_1
E=E_1



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


#In reality there are two cases which need to be distinguished.
#On the one hand, the case where the electron never goes back to the atom : It occurs for an ionization time in [0,T/4] or in [T/2,3T/4] (See the jupyter Notebook for a more precise explanation on how to find these intervalls).
#On the other hand, the case where the electron goes back to the atom : It occurs for an ionization time in [T/4,T/2] or in [3T/4,T]. We can then define a return time t_return
#This can clearly be seen on the precedent plot. When the trajectory crosses the horizontal axe, the electron comes back to the atom. 


########################################################################################################################
# Visualisation of the trajectories for a time of ionization in [0,T/2], the moment when the electron leaves the atom is highlighted with a colored bullet
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

"""Initialization of vectors of ionization time"""

#We can choose to work either on [T/4,T/2] or in [3T/4,T].
#The numerical solution can't be as precise as the theoretical one, the intervalls are thus smaller than [T/4,T/2] or in [3T/4,T].
TPRIME=np.linspace(0.2501*T,0.499*T,100)
TPRIME_bis=np.linspace(0.751*T,0.999*T,100)



##############################################################################################################################
# RETURN TIME AS A FUNCTION OF THE IONIZATION TIME. 
# It should be a strictly decreasing function in the case of a monochromatic source which suggests that there's only one return time for one ionization time. (bijection)
# The theoretical solution is compared to the numerical solution
##############################################################################################################################

#TEST for t' in [0.25,0.5]

"""
for i in range(len(TPRIME)):
    plt.scatter(TPRIME[i],find_Treturn(TPRIME[i],LASER_SOURCE,3*T),marker='+',color='red',lw=3)
    

Tseconde=find_Treturn_theoretical(TPRIME)
plt.plot(TPRIME,Tseconde,label='Theoretical solution',lw=2)
plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Ionization Time $t_i$ ', fontsize=30)
plt.xticks( 0.25*T+np.arange(2)*T*0.25, (r"$\frac{T}{4}$", r"$\frac{T}{2}$"), fontsize=30)
plt.yticks( 0.5*T+np.arange(4)*(T)*0.25, ( r"$\frac{T}{2}$",r"$\frac{3T}{4}$", r"$T$", r"$\frac{5T}{4}$"), fontsize=30)
plt.scatter(0.5,0.5,marker='+',color='red',lw=3, label ='Numerical solution')
plt.ylabel('Return time $t_r$', fontsize=30)
plt.tight_layout()
"""

#TEST for t' in [0.75,1]

"""
for i in range(len(TPRIME)):
    plt.scatter(TPRIME_bis[i],find_Treturn(TPRIME_bis[i],LASER_SOURCE,T),marker='+',color='red',lw=3)

plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Ionization Time $t_i$ ', fontsize=30)
plt.xticks( 0.75*T+np.arange(2)*T*0.25, (r"$\frac{3T}{4}$", r"$T$"), fontsize=30)
plt.yticks( 1*T+np.arange(4)*(T)*0.25, ( r"$T$",r"$\frac{5T}{4}$", r"$\frac{3T}{2}$", r"$\frac{7T}{4}$"), fontsize=30)
plt.scatter(1,1,marker='+',color='red',lw=3, label ='Numerical solution')
plt.ylabel('Return time $t_r$', fontsize=30)
plt.tight_layout()
"""

#################################################################################
# TIME OF FLIGHT AS A FUNCTION OF THE RETURN TIME
#Here is plotted the time of flight (t_return-t_ionization) as a function of the return time.
#TEST for t' in [0.25,0.5]
#################################################################################

"""
plt.figure()
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i],LASER_SOURCE,3*T),find_Treturn(TPRIME[i],LASER_SOURCE,3*T)-TPRIME[i])
    
plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Return Time $t_r$ (T)', fontsize=30)
plt.ylabel('Time of flight $\\tau$', fontsize=30)
plt.yticks( np.arange(5)*T*0.25, (r"$0$",r"$\frac{T}{4}$", r"$\frac{T}{2}$", r"$\frac{3T}{4}$", r"$T$"), fontsize=30)
plt.xticks( 0.5*T+np.arange(4)*(T)*0.25, ( r"$\frac{T}{2}$",r"$\frac{3T}{4}$", r"$T$", r"$\frac{5T}{4}$"), fontsize=30)
plt.tight_layout()
"""
#The same can be done for the interval [3T/4,T] using plt.scatter(find_Treturn(TPRIME_bis[i],LASER_SOURCE),find_Treturn(TPRIME_bis[i],LASER_SOURCE)-TPRIME_bis[i])



#################################################################################
# KINETIC RETURN ENERGY AS A FUNCTION OF THE RETURN TIME
# Here is plotted the Kinetic Return Energy as a function of the Return Time.
#TEST for t' in [0.25,0.5]
##############################################################################
"""
plt.figure()
for i in range(len(TPRIME)):
    plt.scatter(find_Treturn(TPRIME[i],LASER_SOURCE,3*T),16*(np.pi)**2*Kinetic_energy(TPRIME[i],LASER_SOURCE,3*T))
plt.grid() 
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Return Time $t_r$ (T)', fontsize=30)
plt.ylabel('Kinetic Return Energy $\\frac{E_k}{U_p}$', fontsize=25)
plt.plot(np.linspace(0.5,T+0.25,100),3.17*np.ones(100), color='red',label='$E_k=3.17U_p$',lw=2)
plt.xticks( 0.5*T+np.arange(4)*(T)*0.25, ( r"$\frac{T}{2}$",r"$\frac{3T}{4}$", r"$T$", r"$\frac{5T}{4}$"), fontsize=30)
plt.tight_layout()
"""
#The same can be done for the interval [3T/4,T] USING plt.scatter(find_Treturn(TPRIME_bis[i],LASER_SOURCE),find_Treturn(TPRIME_bis[i],LASER_SOURCE)-TPRIME_bis[i])


#################################################################################
# SPECTRUM FOR THE KINETIC RETURN ENERGY
# Here is plotted an histogram for the return kinetic energy which could be seen as a spectrum
# The cut-off is clearly seen
##############################################################################

T_Prime=np.linspace(0.251,0.499,100)
"""
E_kin=Kinetic_energy(T_Prime,LASER_SOURCE,3*T)
Max_Ekin=max(E_kin)
size_of_T=10
Hist=[np.floor(size_of_T*e/Max_Ekin) for e in E_kin]

plt.hist(Hist,normed=1)
plt.xticks( np.arange(5)*5,(r"0", r"$1.58U_p$", r"$3.17U_p$", r"$4.75U_p$", r"$6.34U_p$"), fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel('Kinetic Energy $E_k$' , fontsize=30)
plt.ylabel('Frequency', fontsize=30)
plt.tight_layout()

"""


plt.legend()
plt.show()
