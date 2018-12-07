# phy571_high_harmonics_generation
**numerical physics project**

This repository is a project aiming at studying high harmonics generation by laser pulses.

**Groundwave.py** : Compute the groundstate of the soft-coulomb potential using wick's rotation, then write the resulting array in a separated .npy file.

**Groundstate** : An array representing the groundstate wavefunction of a soft-coulomb potential.

**Integration.py** : used to derive the solution from the movement equation. The **odeint** function has been chosen to perform the integration.The system needs to be written under the vectorial notation y=[x,x'] where one is looking for y'=[x',x'']

**Physical_quantities.py** : definition of the physical quantities in (u.a)

**Useful_functions.py** : contains all the functions which will be usefull for the project. e.g : Electric Field for the laser source, theoretical solution for the movement, etc...

**color.py** : used to build a gradient of N colors between two colors given in the RGB notation.

**trajectory_evolution.py :** used to plot the evolution of the trajectory with respect to the time. The indidence of the ionization time t' on the result is studied. t' is the time when the electron leaves the atomical orbit. It's the beginning of the classical solution. We define here a vector of t' of length Number_of_t_prime.

**find_Treturn.py** : calculation of the return time t'' when it exists.
