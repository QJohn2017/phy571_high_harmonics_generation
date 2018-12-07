# phy571_high_harmonics_generation
**numerical physics project**

This repository is a project aiming at studying high harmonics generation by laser pulses.

**Physical_quantities.py** : definition of the physical quantities in (u.a)

#######################################################################################################################
################################################## CLASSICAL PART #####################################################
#######################################################################################################################


**Useful_functions.py** : contains all the functions which will be usefull for the project. e.g : Electric Field for the laser source, theoretical solution for the movement, kinetic retrun energy, integration of the movement equation using **odeint**, determination of the time when the electron goes back to the atom etc..
All the other files import this one to run correctly.

**Monochromatic_source.py** : File to plot the trajectories, the return time, the time of flight and the kinetic return energy for a monochromatic source E(t)=E0*sin(wt).

**Bichromatic_source.py** : File to plot the trajectories, the return time, the time of flight and the kinetic return energy for a bichromatic source E(t)=E0*sin(wt)+E2*sin(2wt).

**Trichromatic_source.py** : File to plot the trajectories, the return time, the time of flight and the kinetic return energy for a trichromatic source EE(t)=E0*sin(wt)+E2*sin(2wt)+E3*sin(3wt).

**Enveloppe.py** : File to plot the trajectories, the return time, the time of flight and the kinetic return energy for a laser pulse source E(t)=E0*sin(pi*t/tau)^2*sin(wt) with tau the duration of the pulse.


**color.py** : used to build a gradient of N colors between two colors given in the RGB notation.

#######################################################################################################################
################################################## QUANTUM PART #####################################################
#######################################################################################################################

**Groundwave.py** : Compute the groundstate of the soft-coulomb potential using wick's rotation, then write the resulting array in a separated .npy file.

**Groundstate** : An array representing the groundstate wavefunction of a soft-coulomb potential.

