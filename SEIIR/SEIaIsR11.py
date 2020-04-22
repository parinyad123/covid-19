#from collections import Counter
#from pprint import pprint

#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from EpiModel import *


# Asymptomatic individuals are often less infectious than those displaying symptoms by some fraction rbeta
rbeta = 2 #0.75
# A fraction pₐ of all of those Exposed
pa = .7  #0.4 (ไม่เกิน 1)
R0 =1.5  #2.0 (ห้ามน้อยกว่า 1)
epsilon =  0.5 # 0.4 (ห้ามน้อยกว่า 0)
mu = 4 # 0.1  #(ห้ามน้อยกว่า 0)

day = 365
# Number of population
N = 10 #100000 ต้องมากกว่า (S+Is+Ia+R+E)
# Initial number of Asymptomatic
I0 = 0.2 # (ห้ามน้อยกว่า 0)


SEIIR.integrate(day, S=N-I0, Ia=0, Is=I0, E=0, R=0)

beta = R0*mu/(pa*rbeta+(1-pa))

SEIIR = EpiModel()
SEIIR.add_interaction('S', 'E', 'Ia', rbeta*beta)
SEIIR.add_interaction('S', 'E', 'Is', beta)
SEIIR.add_spontaneous('E', 'Ia', epsilon*pa)
SEIIR.add_spontaneous('E', 'Is', epsilon*(1-pa))
SEIIR.add_spontaneous('Ia', 'R', mu)
SEIIR.add_spontaneous('Is', 'R', mu)


Sy, Ey, Iay, Isy, Ry = SEIIR.S/N, SEIIR.E/N, SEIIR.Ia/N, SEIIR.Is/N, SEIIR.R/N

Ef = Ey
Iaf = Iay+Ey
Isf = Isy+Iaf
Sf = Sy+Isf

t=np.arange(0,day-1)

fig= plt.figure(figsize=(15,6))
plt.style.use('dark_background')
plt.fill_between(t, 0, Sf,  label='Susceptible')
plt.fill_between(t, 0, Isf,  label='Symptomatic')
plt.fill_between(t, 0, Iaf,  label='Asymptomatic')
plt.fill_between(t, 0, Ef,  label='Exposed')
plt.fill_between(t, 1-SEIIR.R/N,1,  label='Recovered')
plt.legend(prop={"size":15})

plt.title("Symptomatic/Asymptomatic model", fontsize=18)
plt.xlabel("Day", fontsize=14)
plt.ylabel("Population", fontsize=14)
plt.tick_params(labelsize=15)