#####################################
Playing with Looking at EM spectrum
#####################################
HanishKVC, 20201226

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
max mode amplitude capture over multisecond data in each freq band being
scanned, potentially using rectangle window function, can give a picture
of what frequencies are at play. And potentially to what max amplitude
levels.

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

One can mix and match the options supported by the program to try and
explore the emissions/signals in a suitable way.


Use
#####

Supports two scan modes

Zero Span
===========

kspecanal.py zeroSpan centerFreq <theFreq>

This scans a frequency band centered at centerFreq, and spread over a
band width of 2.4MHz (decided based on the sampling rate limit of rtlsdr),
again and again.

It shows the normalised fft result magnitudes of the scan on a Log scale.
These are shown as 3 curves

        The Red curve which is the max values seen till then.

        The Green curve, which is the average of values seen till then.

        The Blue curve, which is the current scan's values.

It also shows a heatmap | waterfall of the signal levels over a period of
time till then.


Scan
=======

Scan over a large freq range, in steps and show the results both as a
normal signal level plot as well as a heat map with historic data.

The signal level plot contains the Max level seen till then and the
average level seen till then, by default.

Normal
--------

kspecanal.py scan startFreq <theStartFreq> endFreq <theEndFreq>

One can specify a frequency range over which to scan. If the specified
range is larger than what is supported by the hardware in one go, then
it will step through the specified range, in steps.

The normalised fft result is clipped wrt low values (so that the noise
can be clipped to some extent) and then shown on a log scale.

NOTE: If the freq range being scanned isn't a multiple of the sampling
rate, then endFreq will be adjusted to make it a multiple. User will
be alerted about the same and they have to press any key to continue,
in this case.

QuickFullScan
---------------

kspecanal.py quickFullScan

is a alias for

        kspecanal.py scan startFreq 30e6 endFreq 1.5e9

i.e this triggers a quick scan from 30e6 to 1.5e9

FMScan
--------

kspecanal.py fmScan

is a alias for

        kspecanaly.py scan startFreq 88e6 endFreq 108e6

If user doesnt specify any arguments, then the program defaults to this
mode.



UI Buttons
============

Quit - On pressing the Quit button, the btn label changes to QuitWait,
inturn the program finishes the current freq band scan and then exits
the scan loop and changes btn label to QuitPress. User can now either
explore the plots using the pan and zoom buttons in the gui, if they
so desire. Then on pressing any key in the console from where the prg
was started, the program will quit.

Pause - This toggles the pltHighsPause between enable and disable. If
enabled, then user requires to press any key in the console, to step
into next round of scan. Parallely the user can explore the plots
before pressing any key in the console.

Levels - This toggles the bPltLevels between enable and disable.

HeatMap - This toggles the bPltHeatMap between enable and disable.


NOTE
=======

The logic is setup to apply fft on fftSize samples at a time, which is
independent of the samplingRate specified. This in turn controls the fft
bin width | RBW to be around samplingRate/fftSize. Inturn what is shown
on the screen is also controlled by xRes, larger the xRes more finegrained
the amount of data shown on screen, provided the screen resolution is also
equally good.

There is processing and plotting delay between the repeating scans, so
any signal occuring at that time will be lost. Similarly when using scan
to scan through a large freq range (especially when doing beyond 2.4MHz
band) at any given time only a freq band equivalent to samplingRate is
what is being monitored, so any signals occuring in any other bands at
that time will not be captured.

If there is a error in setting up the sdr, then the value of that freq
band gets set to all 1s, this inturn leads to a level of around -25 or
so in the levels plot.

For real signal the curscan flow maintains the signal levels; while for
complex iq signal data, curscan flow adds 3dB to signal levels. Also
dont forget that the default pltCompress of Avg, eats into the siglevels
in general.


Other Args
-----------

samplingRate <samplingRateFloat>

        Default 2.4e6; this is a good value for rtlsdr. If you want,
        you can reduce it.

minAmp4Clip <float>

        Default (1/256)*0.33; Change it to control the forced noise floor.
        Any measured signal level below this in the freq domain will be
        set to this value.

gain <gainFloat>

        Default 19.1; Increase or reduce this depending on the strength
        of the signals being studied.

scanRangeCumuMode <Avg|Max|Raw>

        Default Avg; Change to Max, if one wants to know the max value
        noticed at any time during the scan.

window <true|false>

        Default: False; Controls whether a windowing function (hanning)
        is applied to the time domain samples, before fft is done. Helps
        get a better sense about the signals in a scan. Useful if only
        a limited scan is being done. However for small fft window size,
        overlapped sliding may be more useful.

fftSize <integer>

        Default: 2**14; The number of samples that is run through the fft
        in one go. This also decides the resolution bandwidth of the logic.
        Larger the fftSize, finer the freq resolution. Needs to be a power
        of two value, or else multiple of xRes.

curScanNonOverlap <float>

        Default: 0.1; As the small size fft window is slide over a larger
        signal sample dataset, this controls how much of the data is
        skipped during the overlapping. 0.1 means 90% overlapping 1.0
        means 0% overlapping. Overlapping normally helps get a better feel of
        the signal level, even thou only a fraction of a second worth of data
        is run through fft at a time.

curScanCumuMode <Avg|Max|Raw>

        Default Avg; Change to Max, if one wants to know the max value
        noticed at any time during the scan.

bPltLevels <true|false>

        Default: True; Control whether the current internal scan signal level
        is plotted or not. Disabling this will speed up the scan interval a bit.

bPltHeatMap <true|false>

        Default: True; Control whether the signal level history | heat map is
        plotted or not. Disabling this will speed up the scan interval a bit.

scanRangeNonOverlap <float>

        Default: 0.75; Change to control how much of the freq band is overlapped
        as the scan range logic scans/steps through a given range of frequencies.
        Set it to 1.0 to avoid overlapping, or set it to 0.5 to overlap 50% of the
        freq band, as the logic tunes to the next center freq to scan the next
        adjacent freq band. Could help overcome any non linearity in measuring
        within a freq band, to an extent.

prgLoopCnt <int>

        Default: A large value; Change to a smaller value, if you want to scan
        for a short amount of time like few minutes or so. As zooming or panning
        the plot, when the program is running and updating the plot is not easy
        and consistent, so one can scan for a short time, and then once the scan
        is finished look into the scan plot in detail, or else one will have to
        wait till the program stops after a long time.

pltCompress <Raw|Avg|Max>

        Default: Average; This allows one to control how finegrained or not is
        the signal levels across adjacent freqs that are shown. This along with
        fftSize and xRes, decides how finegrained is the freq resolution you see
        on the screen. NOTE: Using Avg will smooth the display, but will impact
        the signal levels seen.

xRes <int_poweroftwovalue>

        Default: 512; This controls the horizontal resolution (number of data
        points related to frequencies or groups of adjacent frequencies) of the
        data passed to the plotting logic. This needs to be a power of two value,
        or else a sub multiple of fftSize.

        To ensure that any signal data is not lost wrt the heatmap display, set
        this to match the actual display resolution of the heatmap on your screen.
        Or even a value which is smaller than the actual screen resolution of the
        heatmap will also do, but you will lose some amount of freq resolution wrt
        display. Based on your situation, you may be able to increase this value,
        and still not lose any value.

pltHighsNumMarkers <int>

        Default: 5; Control how many markers should be shown in the plot, wrt
        the high signal levels.

pltHighsDelta4Marking <float>

        Default: 0.025; Specify how much fraction of the plot's full freq range,
        is used as the delta needed between marked frequencies, when deciding
        whether to mark the high signal level freq on the plot or not.

pltHighsPause <boolean>

        Default: False; Specify whether the scan range plot should pause after
        each scan of the specified range of frequencies. THis allows the user
        to see the list of high signal level frequencies, on the plot.
        Independent of above, the list of high siglevel freqs is also printed
        on the console.


NOTE: Do look into the source to get the latest | current default setting for the
different options, and or to change as one sees fit.



Signal level display
------------------------

For more representative signal level display, use the following property values

ZeroSpan mode

        pltCompress raw <OR ELSE> pltCompress max

Scan mode

        # Start with avg to get a rough overview

        pltCompress avg

        # Switch to conv to get a more representative view

        pltCompress conv

        # Then use max or raw to get the more practical view

        pltCompress raw <OR ELSE> pltCompress max

        # U can also add scanRangeNonOverlap to the mix

        scanRangeNonOverlap 1.0

To ensure that heatmap doesnt eat up any signal data, set the xRes to match the
actual screen resolution of the heatmap and or lesser than it.

NOTE: HeatMap by default uses pltCompressHM mapped to Max logic for its data and
is Not user controllable from commandline.



TODO
#######

Account -ve freqs of complex iq fft. [Done]

Put something similar to old dwelltime, but controlled using rbw
rather than dwell time. Along with windowing and some amount of limited
sliding. [5050]

Add Max based cumulation of fft result and provide option to switch
between average and Max [Done].

Add the running heatmap/waterfall view [Done].

Overlap across scan bands [Done].

Use pygame or cairo or .. to do the plots. Heatmap with large freq bands and
default or large fftSize, could bring the program and the system to its knees.
And or parallely save into image with sufficient resolution. Also the imshow,
losses signal info, if the signal is surrounded by very weak or no signal in
the adjacent frequencies. Need to use implement my own logic, with max instead
of averaging when mapping multiple data points into individual pixels. [Done
Rather process the data by merging adjacent data points, before plotting them]

