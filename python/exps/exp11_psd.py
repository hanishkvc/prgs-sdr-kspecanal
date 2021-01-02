#!/usr/bin/env python3
# Test PSD
# HanishKVC, 2021

import numpy as np
import matplotlib.pyplot as plt

times = np.linspace(0,1,512)

def gen_freq(freq, amp):
    s = amp*np.sin(2*np.pi*freq*times) + amp*np.cos(2*np.pi*freq*times)*1j
    return s

s1 = gen_freq(2, 1)
s2 = gen_freq(20, 0.25)
s3 = gen_freq(40, 1)
s = s1 + s2 + s3

normalise = 1
normalise = len(s)
win = np.hanning(len(s))
sW = s * win
winAdj = len(win)/np.sum(win)
ampSpec = np.fft.fftshift(np.abs(np.fft.fft(sW)/normalise))
ampSpec = ampSpec * winAdj
powSpec = ampSpec**2

psd = plt.psd(s, NFFT=len(s), window=None)
plt.show()
magSpecgram = plt.specgram(s, NFFT=512, window=None, Fs=2, mode='magnitude')
plt.show()

plt.subplot(4,2,1)
plt.plot(ampSpec)
plt.subplot(4,2,2)
plt.plot(10*np.log10(ampSpec))

plt.subplot(4,2,3)
plt.plot(powSpec)
plt.subplot(4,2,4)
plt.plot(10*np.log10(powSpec))

plt.subplot(4,2,5)
plt.plot(psd[0])
plt.subplot(4,2,6)
plt.plot(10*np.log10(psd[0]))

ampSpec = winAdj*np.fft.fftshift(np.abs(np.fft.fft(sW)))
powSpec = ampSpec**2
plt.subplot(4,2,7)
print(max(powSpec))
plt.plot(powSpec/(512*512))
plt.subplot(4,2,8)
plt.plot(magSpecgram[0])
plt.show()

