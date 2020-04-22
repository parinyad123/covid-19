
import numpy as np
import matplotlib.pyplot as plt
from EpiModel import *


# Asymptomatic individuals are often less infectious than those displaying symptoms by some fraction rbeta
rbeta = 0.5 #0.75
# A fraction of all of those Exposed
pa = .5  #0.4 (ไม่เกิน 1) (pa is between zero to one)
R0 = 1.5  #2.0
epsilon =  0.5 # 0.4 (ห้ามน้อยกว่า 0) (epsilon is more than zero)
mu = 0.5 # 0.1  #(ห้ามน้อยกว่า 0) (mu is more than zero)
# Number of population
N = 100 #100000 ต้องมากกว่า (N = S+Is+Ia+R+E)
# Initial number of symptomatic
Is_0 = 1 # (ห้ามน้อยกว่า 0) (I0 is more than zero)
E_0 = 4
R_0 = 2
Ia_0 = 0

day = 150

beta = R0*mu/(pa*rbeta+(1-pa))

SEIIR = EpiModel()
SEIIR.add_interaction('S', 'E', 'Ia', rbeta*beta)
SEIIR.add_interaction('S', 'E', 'Is', beta)
SEIIR.add_spontaneous('E', 'Ia', epsilon*pa)
SEIIR.add_spontaneous('E', 'Is', epsilon*(1-pa))
SEIIR.add_spontaneous('Ia', 'R', mu)
SEIIR.add_spontaneous('Is', 'R', mu)

SEIIR.integrate(day, S=N-Is_0-Ia_0-E_0-R_0, Ia=Ia_0, Is=Is_0, E=E_0, R=R_0)

Sy, Ey, Iay, Isy, Ry = SEIIR.S/N, SEIIR.E/N, SEIIR.Ia/N, SEIIR.Is/N, SEIIR.R/N

Ef = Ey
Iaf = Iay+Ey
Isf = Isy+Iaf
Sf = Sy+Isf

t=np.arange(0, day-1)
plt.style.use('dark_background')
fig= plt.figure(figsize=(15,6))
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
plt.show()
