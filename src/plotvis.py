#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def loadvis(filename):
    f=open(filename,"rb")
    nB=int.from_bytes(f.read(4),'little')
    nL=int.from_bytes(f.read(4),'little')
    return (nB,nL,np.reshape(np.frombuffer(f.read(nB*nL*16),"<f8",count=nB*nL*2),(nB,nL,2)))

v=loadvis("vis.dat")

figure,axis=plt.subplots(2,1)
for i in range(v[0]):
    axis[0].plot(v[2][i,:,0]*v[2][i,:,0]+v[2][i,:,1]*v[2][i,:,1])
    axis[1].plot(np.arctan2(v[2][i,:,0],v[2][i,:,1]))

axis[1].set(xlabel='Wavelength direction')

plt.show()


