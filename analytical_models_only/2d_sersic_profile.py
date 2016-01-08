import numpy as np

def sersic_2d(x1,x2,Ie,Re,bn,n):
    R = np.sqrt(x1*x1+x2*x2)
    res = Ie*np.exp(-bn*((R/Re)**(1.0/n)-1.0))
    return res
