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

The log plot is shown at the end.

TODO
#######

Add Max based cumulation of fft result and provide option to switch
between average and Max.

Add the running heatmap/waterfall view.

Overlap across scan bands.

