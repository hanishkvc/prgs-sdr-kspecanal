#!/usr/bin/env python3
# Generate samples for testing FFT logic of kspecanal
# HanishKVC, v20201226
#


import sys
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


    def rel_freqs(self):
        '''
        Generate frequencies which are fixed relative to the configured
        center frequency.
        Avoids generating freq which falls into 0th fft bin.
        '''
        f1 = self.sample_rate/4
        f2 = 1
        f3 = -self.sample_rate/4
        return [ f1, f2, f3 ]


    def abs_freqs(self):
        '''
        Generate frequencies which are fixed in a absolute term, independent
        of the configured center frequency.

        It generates signals at each megahz.
        Generates freq, even if it falls in the 0th fft bin/position.
        '''
        startF = self.center_freq - self.sample_rate/2
        endF = self.center_freq + self.sample_rate/2
        sR = int(np.ceil(startF/1e6)*1e6)
        eR = int((endF//1e6)*1e6)+1
        f = []
        for cur in range(sR, eR, int(1e6)):
            fCur = self.center_freq - cur
            if fCur == 0:
                fCur = 0
            print("DBUG:absFreqs:{}={}".format(cur,fCur))
            f.append(fCur)
        return f


    def read_samples(self, size):
        '''
        Actual rtlsdr returns iq data in complex notation.
        This currently returns real data
        '''
        gainMult = 10**(self.gain/10)
        dur = size/self.sample_rate
        tStart = random.random()
        tTimes = np.linspace(tStart, tStart+dur, self.sample_rate*dur)
        sT = []
        f = np.array(self.rel_freqs())
        f = np.array(self.abs_freqs())
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
        win = np.kaiser(len(dataBuf), 15)
        win = np.hanning(len(dataBuf))
        freqs = np.fft.fftshift(np.fft.fftfreq(len(dataBuf),1/self.sample_rate)) + self.center_freq
        print("INFO:testfft_rtlsdr: freqs [{}] - [{}], fftSize [{}]".format(min(freqs), max(freqs), len(dataBuf)))
        fftCur = np.abs(np.fft.fft(dataBuf))/len(dataBuf)
        fftCur = np.fft.fftshift(fftCur)
        fftCurWin = np.abs(np.fft.fft(dataBuf*win))/len(dataBuf)
        fftCurWin = np.fft.fftshift(fftCurWin)
        plt.subplot(2,2,1)
        plt.plot(freqs, fftCur)
        plt.subplot(2,2,2)
        plt.plot(freqs, 10*np.log10(fftCur))
        plt.subplot(2,2,3)
        plt.plot(freqs, fftCurWin)
        plt.subplot(2,2,4)
        plt.plot(freqs, 10*np.log10(fftCurWin))
        plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        centerFreq = float(sys.argv[1])
    else:
        centerFreq = 92e6
    sdr = RtlSdr(fC=centerFreq)
    sdr.test()


