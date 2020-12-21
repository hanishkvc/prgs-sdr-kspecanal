#####################################
Playing with Looking at EM spectrum
#####################################
HanishKVC, 20201216

Gist: this logic, rtlsdr ...
################################

While design of a product, one is interested in understanding

* what and all frequency emissions are there from a given product
  and or at a given location in the product.

  Kind of probe/antenna used and at what distance, helps decide
  the locality or globality of the study results.

* how changing some aspect of the design or its setup inturn improves
  or worsens the situation. Relative study is good enough at one level.


So using a full-samplingrate-width window based overlapped sliding with
max amplitude captured over multisec data in each freq band being scanned,
potentially using rectangle window function, can give a picture of what
frequencies are at play. And potentially to what max amplitude levels.

While running the same with a reduced fft window size and sliding, will
allow one to scan through with less fidelity in a fast way.

And using water-fall/heat-map view of fft outputs of smoothing windowed
(hanning/kaiser) data, as one slides over the full data in a overlapped
manner, looking at a subset of the data at each given step, should help
guage how the frequencies involved are coming and going out of existance
with time.

ie getting a rought picture of all the frequencies involved and to what
relative amplitude levels and how they roughly change relatively on their
own and with changes to design or setup is what this tries to give at a
simple level. It doesnt go into nitty gritties beyond it.

It is a study in rough relatives and not absolutes,
and that is good enough / sufficiently useful, many a times ;-)


Use
#####

Supports two scan modes

Zero Span
===========

kspecanal.py zerospan centerFreq <theFreq>

This scans a frequency band centered at centerFreq, and spread over a
band width of 2.4MHz (decided based on the sampling rate limit of rtlsdr),
again and again.

It shows the normalised fft result magnitudes of the scan on a Log scale.

It also shows a heatmap | waterfall of the signal levels over a period of
time till then.


Scan
=======

kspecanal.py scan startFreq <theStartFreq> endFreq <theEndFreq>

One can specify a frequency range over which to scan. If the specified
range is larger than what is supported by the hardware in one go, then
it will step through the specified range, in steps.

The normalised fft result is clipped wrt low values (so that the noise
can be clipped to some extent) and then shown on a log scale.



NOTE
=======

Currently the logic is setup to apply fft on 2**14 samples at a time,
this gives a fft bin width / RBW of around 150Hz for 2.4e6 sampling rate.

Other Args
-----------

samplingRate <samplingRateFloat>
minAmp4Clip <float>
gain <gainFloat>
cumuMode <Avg|Max|Copy>
window <true|false>
fftSize <integer>
nonOverlap <float>




TODO
#######

Account -ve freqs of complex iq fft. [almost]

Put something similar to old dwelltime, but controlled using rbw
rather than dwell time. Along with windowing and some amount of limited
sliding. [5050]

Add Max based cumulation of fft result and provide option to switch
between average and Max.

Add the running heatmap/waterfall view.

Overlap across scan bands.

