#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import frame

f=open("frame.dat","rb")
d=frame.load(f)[0]
f.close()

print("nL=",d['nL']," nF=",d['nF']);

plt.imshow(d['image'])#,aspect='auto',interpolation='none')
plt.show()

plt.plot(d['image'][99,:])
plt.plot(d['image'][0,:])
plt.show()

