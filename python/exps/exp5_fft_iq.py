#!/usr/bin/env python3

from numpy import *
from matplotlib.pyplot import *

#interactive(True)
#show()
#clf()

subplot(3,2,1)
sr = 20
t = linspace(0,1,sr)
s8 = sin(2*pi*8*t)
s1 = sin(2*pi*1*t)
tS = s1 + s8
tFft = abs(fft.fft(tS))
np.round(tFft)
plot(tFft)

tFreqs = fft.fftfreq(len(tS),1/sr)
#fft.fftfreq(200,1/200)
#fft.fftfreq(20,1/200)
subplot(3,2,2)
plot(tFreqs, tFft)



#### IQ go full hog wrt amp

# A+B, Sqrt(A^2+B^2), Amplitude from IQ
subplot(3,2,3)
a = linspace(0,10,11)
plot(sqrt(a**2+a**2))
plot(a+a)
legend(["hyp", "a+a"])


# IQ
t = linspace(0,pi,20)
tI=cos(2*pi*1*t)
tQ=sin(2*pi*1*t)
subplot(3,2,5)
plot(tI)
plot(tQ)
subplot(3,2,6)
plot(tI+tQ)
plot(sqrt(tI**2+tQ**2))
legend(["i+q", "hyp"])
show()

#input("Press any key to continue...")

