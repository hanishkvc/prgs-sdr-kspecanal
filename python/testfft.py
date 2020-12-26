#!/usr/bin/env python3
# Generate samples for testing FFT logic of kspecanal
# HanishKVC, v20201226
#


import random
import numpy as np
import matplotlib.pyplot as plt


class RtlSdr():
    '''
    A simple simulated rtlsdr
    '''

    def __init__(self, fC = 92e6, fS = 2.4e6, gain = 0.5):
        self.center_freq = fC
        self.sample_rate = fS
        self.gain = gain


    def read_samples(self, size):
        '''
        Actual rtlsdr returns iq data in complex notation.
        This currently returns real data
        '''
        gainMult = 10**(self.gain/10)
        f1 = self.sample_rate/4
        f2 = 1
        f3 = -self.sample_rate/4
        dur = size/self.sample_rate
        tStart = random.random()
        tTimes = np.linspace(tStart, tStart+dur, self.sample_rate*dur)
        sT = []
        f = np.array([ f1, f2, f3 ])
        print("INFO:testfft_rtlsdr: freqs [{}], sampRate [{}], tStart [{}], dur [{}], len [{}]".format(self.center_freq - f, self.sample_rate, tStart, dur, len(tTimes)))
        s = np.zeros(len(tTimes), dtype=complex)
        for i in range(len(f)):
            sS = gainMult * np.sin(2*np.pi*f[i]*tTimes)
            sC = gainMult * np.cos(2*np.pi*f[i]*tTimes)
            sT.append(sS + sC*1j)
            s += sT[i]
        return s


    def close(self):
        pass


    def test(self):
        dataBuf = self.read_samples(random.randint(1,10)*2**16)
        freqs = np.fft.fftshift(np.fft.fftfreq(len(dataBuf),1/self.sample_rate)) + self.center_freq
        print("INFO:testfft_rtlsdr: freqs [{}] - [{}], fftSize [{}]".format(min(freqs), max(freqs), len(dataBuf)))
        fftCur = np.abs(np.fft.fft(dataBuf))/len(dataBuf)
        fftCur = np.fft.fftshift(fftCur)
        plt.subplot(2,1,1)
        plt.plot(freqs, fftCur)
        plt.subplot(2,1,2)
        plt.plot(freqs, 10*np.log10(fftCur))
        plt.show()


if __name__ == "__main__":
    sdr = RtlSdr()
    sdr.test()


