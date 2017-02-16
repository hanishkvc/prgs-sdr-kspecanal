#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

sr = 10e3
totalTime = 2
t = np.arange(0+0.137, totalTime+0.137, 1/sr)

s10k = np.sin(2*np.pi*1e3*t)
s20k = np.sin(2*np.pi*2e3*t)
s17k = np.sin(2*np.pi*1.7e3*t)

s=s10k+s20k+s17k

def do_fft(s, startFreq, endFreq, sr):
    f=np.fft.fft(s)
    f=abs(f)
    f=f/len(f)
    f=f[0:len(f)/2]
    fr=np.linspace(startFreq, endFreq,len(f))
    plt.plot(fr,f)
    plt.grid()
    plt.title("sampleLen[{}]".format(len(s)))

def exp_fft(s, startFreq, endFreq, sr):
    for i in range(1,9):
        plt.subplot(4,2,i)
        do_fft(s[0:len(s)/i], startFreq, endFreq, sr)
    plt.show()

exp_fft(s, 0, sr/2, sr)

