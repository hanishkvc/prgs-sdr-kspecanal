#!/usr/bin/env python3
# Test impact of partial signal cycle and windowing
#

import numpy as np
import matplotlib.pyplot as plt


# Setup parameters
sr = 80
f1 = 2
f2 = 8
f3 = 38
startTime = 0
endTime = 10


# Generate time units
t = np.arange(startTime, endTime, 1/sr)
t = np.linspace(startTime, endTime, sr*(endTime-startTime))


# set plot params
nRows = 11
nCols = 4
plt.figure(figsize=(4*nCols, 3*nRows))


def plot_fft(fftN):
    '''
    Plot fft data, till Nyquist freq
    Also match it to its corresponding frequencies
    '''
    fftN = fftN[:int(len(fftN)/2)]
    freqs = np.linspace(0, sr/2, len(fftN))
    plt.grid(True)
    plt.plot(freqs, fftN, ".-")


# Generate the signal
s1 = np.sin(2*np.pi*f1*t)
s2 = np.sin(2*np.pi*f2*t)
s3 = np.sin(2*np.pi*f3*t)
s = s1 + s2 + s3
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
plt.plot(np.kaiser(sr, 15))


rawFftMaxMin = 10
rawFftMaxMax = -10
winFftMaxMin = 10
winFftMaxMax = -10
# Check Fft with X.y cycles - ie partial cycle impact
# as the test signal starts at time 0, so also full seconds correspond to full cycles.
# while partial seconds correspond to partial cycle of the signal in the test sample.
fullSecs = 1
for i in range(0, 10):
    endIndex = int(sr*(fullSecs+(i/20)))
    # Plot Raw signal and fft
    sT = s[0:endIndex]
    print(endIndex, len(sT))
    fftN = np.abs(np.fft.fft(sT)/len(sT))
    plt.subplot(nRows, nCols, 5+i*4)
    plt.plot(sT)
    plt.subplot(nRows, nCols, 6+i*4)
    plot_fft(fftN)
    if rawFftMaxMin > max(fftN):
        rawFftMaxMin = max(fftN)
    if rawFftMaxMax < max(fftN):
        rawFftMaxMax = max(fftN)
    print(max(fftN))
    # Plot Windowed signal and fft
    #win = np.hanning(len(sT))
    win = np.kaiser(len(sT),15)
    sT = sT*win
    fftN = np.abs(np.fft.fft(sT)/len(sT))
    plt.subplot(nRows, nCols, 7+i*4)
    plt.plot(sT)
    plt.subplot(nRows, nCols, 8+i*4)
    #plot_fft(fftN*2.3)
    plot_fft(fftN*(len(sT)/win.sum()))
    if winFftMaxMin > max(fftN):
        winFftMaxMin = max(fftN)
    if winFftMaxMax < max(fftN):
        winFftMaxMax = max(fftN)
    print(max(fftN))
print("Raw", rawFftMaxMax/rawFftMaxMin)
print("Win", winFftMaxMax/winFftMaxMin)


# Show the plots
plt.savefig("/tmp/exp4_fft.png")
plt.show()

