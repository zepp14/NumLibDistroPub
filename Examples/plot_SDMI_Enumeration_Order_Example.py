#This example implements a 5 invader SDMI game with a Brute force Capture order
import sdmiFuncLib as sfl
import numpy as np
import matplotlib.pyplot as plt


#Invader Locations

#I0 => [-12.25709932   4.44973388]
#I1 => [ -6.6489077   21.23919568]
#I2 => [  0.41996017  24.16072049]
#I3 => [  8.51769598  11.40288417]

I_array = [[-12.25709932,   4.44973388],
           [ -6.6489077,   21.23919568],
           [  0.41996017,  24.16072049],
           [  8.51769598,  11.40288417],
           [  -6.6489077,  24.16072049]]

## Game Properties

r =0    #Capture range
alpha = 4  #Speed ratio
Defender = [0,0]  #Defender Position
G_circ = [0, -8,  8] #Target Region Properties

## Compute Capture order With Enumeration
bestOrder,_,meanVect, StdData,_ = sfl.computeBestOrderEnum(I_array, r, alpha, Defender, G_circ,eMethod=0)
print("Capture Order: ", bestOrder)
print("Number of Computations", sfl.NumeralCounter)
I_dict,_ = sfl.createInvDict(I_array)

## Plot Trial
fig, axs = plt.subplots()
TotEff, P_star, T = sfl.plotTrial(Defender, I_dict, bestOrder, G_circ, alpha, r, axisSz=[-10,10,-5,10], plotFlag=1, T = 0, effMethod = 0, plt=axs)

FL_score = sfl.flightPathScore(P_star, Defender)
NumInTarget_score = sfl.numInTarget(P_star, G_circ, bestOrder, I_dict)

print('Effeciency Score: ', TotEff)
print('Flight Path Score: ', FL_score)
print('Number of Invaders Captured: ', len(I_array ) - NumInTarget_score)

axs.set_aspect('equal', adjustable='datalim')
plt.show()
