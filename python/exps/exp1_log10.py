#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

dataMax=1
dataMinRef=dataMax/255
t1=np.linspace(0, dataMax, 1e6)
t2=t1+dataMinRef*0.33
t3=t1+dataMinRef*0.50
t4=t1+dataMinRef*1.00
t5=t1+dataMinRef*2.00
t6=t1+dataMinRef*4.00
f1=10*np.log10(t1/dataMinRef)
f2=10*np.log10(t2/dataMinRef)
f3=10*np.log10(t3/dataMinRef)
f4=10*np.log10(t4/dataMinRef)
f5=10*np.log10(t5/dataMinRef)
f6=10*np.log10(t6/dataMinRef)

end=1e5
#plt.plot(t1,f1,"+",t1,f2,"*")
plt.plot(t1[:end],f1[:end],"g", t1[:end],f2[:end],"y", t1[:end], f3[:end],"r")
plt.plot(t1[:end],f4[:end],"g--", t1[:end],f5[:end],"y--", t1[:end], f6[:end],"r--")
plt.show()

