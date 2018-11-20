#This file is used to build a gradient of color

from __future__ import division
#!/usr/local/bin/python
# coding: latin-1
import os, sys
import matplotlib.pyplot as plt
import numpy as np

#This function wil be used to realise a rainbow of colours

def echelle(color_begin, color_end, n_vals):
    r1, g1, b1 = color_begin
    r2, g2, b2 = color_end
    degrade = []
    etendue = n_vals-1
    for i in range(0,n_vals):
       alpha = 1-i/etendue
       beta = i/etendue
       r = int(r1 * alpha + r2 * beta)
       g = int(g1 * alpha + g2 * beta)
       b = int(b1 * alpha + b2 * beta)
       degrade.append((r, g, b))
    return degrade

#Definition of the basis colour
jaune = (255, 255, 0)
rouge = (255, 0, 0)
bleu = (0, 0, 255)
vert = (0,255,0)


#function to convert RGB colors in Hexadecimal
def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb



#TEST des fonctions construites

"""t=np.linspace(0,1,100)
N=6
ECHELLE=echelle(jaune, bleu, N)
for i in range(1,N+1):
    plt.plot(t,t**i,color=rgb_to_hex(ECHELLE[i-1]),label='t** %s'%(i))
    plt.legend()
plt.show()"""



