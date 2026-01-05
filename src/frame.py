import numpy as np

def load(f):
    b=f.read();
    nL=np.frombuffer(b,dtype=np.uint32,count=1,offset=11)[0]
    nF=np.frombuffer(b,dtype=np.uint32,count=1,offset=15)[0]
    return np.frombuffer(b,dtype=np.dtype([('yr',np.uint16),('mo',np.uint8),
                                           ('dy',np.uint8),('hr',np.uint8),
                                           ('mn',np.uint8),('se',np.uint8),
                                           ('ns',np.uint32),('nL',np.uint32),
                                           ('nF',np.uint32),
                                           ('image',np.float64,(nL,nF))]))
