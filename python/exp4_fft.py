#!/usr/bin/env python3
# Test impact of partial signal cycle and windowing
#

import numpy as np
import matplotlib.pyplot as plt


# Setup parameters
sr = 80
f = 2
startTime = 0
endTime = 10


# Generate time units
t = np.arange(startTime, endTime, 1/sr)
t = np.linspace(startTime, endTime, sr*(endTime-startTime))


# set plot params
nRows = 11
nCols = 4


def plot_fft(fftN):
    '''
    Plot fft data, till Nyquist freq
    Also match it to its corresponding frequencies
    '''
    fftN = fftN[:int(len(fftN)/2)]
    freqs = np.linspace(0, sr/2, len(fftN))
    plt.plot(freqs, fftN)


# Generate the signal
s = np.sin(2*np.pi*f*t)
print("NumOfSamples:", len(s))
plt.subplot(nRows, nCols, 1)
plt.plot(s)


# Get the Normalised Fft
fftN = np.abs(np.fft.fft(s)/len(s))
plt.subplot(nRows, nCols, 2)
plot_fft(fftN)

# Plot the window functions
plt.subplot(nRows, nCols, 3)
plt.plot(np.kaiser(sr, 0))
plt.plot(np.hanning(sr))
plt.plot(np.kaiser(sr, 8))


# Check Fft with X.y cycles - ie partial cycle impact
for i in range(0, 10):
    endIndex = int(sr*(5+(i/20)))
    # Plot Raw signal and fft
    sT = s[0:endIndex]
    print(endIndex, len(sT))
    fftN = np.abs(np.fft.fft(sT)/len(sT))
    plt.subplot(nRows, nCols, 5+i*4)
    plt.plot(sT)
    plt.subplot(nRows, nCols, 6+i*4)
    plot_fft(fftN)
    print(max(fftN))
    # Plot Windowed signal and fft
    win = np.hanning(len(sT))
    sT = sT*win
    fftN = np.abs(np.fft.fft(sT)/len(sT))
    plt.subplot(nRows, nCols, 7+i*4)
    plt.plot(sT)
    plt.subplot(nRows, nCols, 8+i*4)
    plot_fft(fftN)
    print(max(fftN))


# Show the plots
plt.show()
