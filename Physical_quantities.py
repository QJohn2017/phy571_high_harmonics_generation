# definition of the physical quantities in (u.a)
import numpy as np

"Charge of the electron"
q = -1  

"Mass of the electron"             
m = 1 

"Amplitudes of the incident function"                 
E0=1 
E2=1
E3=1
E0_env=0.0001 
  
"Period of the oscillation"               
T=1 
  
"Proper frequency"                  
omega=2*np.pi/T
omega_env=0.05

T_env=(2*np.pi)/omega_env

"Duration of the pulse"
tau=800

"Phase"	         
phi1=0													
phi2=0	
phi3=0	


#Numerical constants

Number_of_points=100000
