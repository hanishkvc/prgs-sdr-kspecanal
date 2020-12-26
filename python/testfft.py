# Generate samples for testing FFT logic of kspecanal
# HanishKVC, v20201226
#


import random
import numpy as np


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
        f1 = self.center_freq - self.sample_rate/4
        f2 = self.center_freq
        f3 = self.center_freq + self.sample_rate/4
        dur = size/self.sample_rate
        tStart = random.random()
        tTimes = np.linspace(tStart, tStart+dur, self.sample_rate*dur)
        sT = []
        f = [ f1, f2, f3 ]
        print("INFO:testfft_rtlsdr: freqs [{}], sampRate [{}], dur [{}], len [{}]".format(f, self.sample_rate, dur, len(tTimes)))
        s = np.zeros(len(tTimes), dtype=complex)
        for i in range(len(f)):
            sS = np.sin(2*np.pi*f[i]*tTimes)
            sC = np.cos(2*np.pi*f[i]*tTimes)
            sT.append(sS + sC*1j)
            s += sT[i]
        return s


    def close(self):
        pass


