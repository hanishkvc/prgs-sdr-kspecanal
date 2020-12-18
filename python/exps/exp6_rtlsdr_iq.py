#!/usr/bin/env python3

from numpy import *
from matplotlib.pyplot import *
import rtlsdr


# FM channels available here around this freq
# 91.1 91.9 92.7 93.5 94.3 95.0
# Setting Fc to 91.7e6,
#   should see 91.9 and 92.7 to one side
#   and 91.1 to the other side.
#   As this will appear asymmetric around the mid point,
#     2 signals to one side, 1 freq/signal to other side
#   so should be able to make out and cross check visually.
#

centerFreq = 91.7e6
samplingRate = 2.4e6

rs = rtlsdr.RtlSdr()
rs.sample_rate = samplingRate
rs.center_freq = centerFreq
rs.gain = 7.1
samples = rs.read_samples(16*1024)
samples = rs.read_samples(4096)
samples = rs.read_samples(256*1024)
rs.close()

startFreq = centerFreq - samplingRate/2
endFreq = centerFreq + samplingRate/2

print("INFO: startFreq[{}] centerFreq[{}] endFreq[{}], samplingRate[{}]".format(startFreq, centerFreq, endFreq, samplingRate))

freqs = fft.fftfreq(len(samples),1/samplingRate) + centerFreq
print("min[{}] max[{}]".format(min(freqs), max(freqs)))

#plot(centerFreq+freqs, f)
#show()

win = kaiser(len(samples),6)

fR=abs(fft.fft(real(samples)))
fI=abs(fft.fft(imag(samples)))
fC=abs(fft.fft(samples)) # The thing I missed out and or didnt try for what ever reason 2+ years back or maybe wanted to avoid complex fft or ...
fCW=abs(fft.fft(samples*win))
fA=abs(fft.fft(abs(samples)))
fP=abs(fft.fft(real(samples)+imag(samples)))
fR[0] = 0
fI[0] = 0
fC[0] = 0
fCW[0] = 0
fA[0] = 0
fP[0] = 0
subplot(6,1,1)
plot(fR)
title("Real")
subplot(6,1,2)
plot(fI)
title("Imag")
subplot(6,1,3)
plot(fC)
title("Comp")
subplot(6,1,4)
plot(fCW)
title("CompWin")
subplot(6,1,5)
plot(fA)
title("Abs")
subplot(6,1,6)
plot(fP)
title("I+Q")
show()

#help(fft.fftshift)

