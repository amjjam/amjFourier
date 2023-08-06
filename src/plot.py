#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def loadframe(filename):
    f=open(filename,"rb")
    nL=int.from_bytes(f.read(4),'little')
    nF=int.from_bytes(f.read(4),'little')
    return np.reshape(np.frombuffer(f.read(nL*nF*8),"<f8",count=nL*nF),(nL,nF))

frame=loadframe("frame.dat")

figure,axis=plt.subplots(2,1,gridspec_kw={'height_ratios':[3,1]})

image=axis[0].imshow(frame[:,100:220],aspect='auto',interpolation='none',extent=[100,220,256,0])
axis[0].set(ylabel="Wavelength direction")
axis[1].plot(np.arange(100,220),frame[150,100:220])
axis[1].set(xlabel="Fringe direction")
#plt.colorbar(image,fraction=0.028,cax=axis[0,1])
plt.show()

