## OBSclockdrift
This package estimates the clock drift rate of ocean bottom seismometers using ambient noise cross-correlation and applies the correction to raw miniseed files. The methods used here are based on Hable et al. 2018 (GJI) [doi:10.1093/gji/ggy236](https://academic.oup.com/gji/article/214/3/2014/5038378)

In order to successfully measure clock drift at a station of interest, a reference station with correct timing is needed. It is best to use a reference station that is as close as possible to the station of interest in order to maximize signal to noise (however, we have had luck using a reference station up to 200 km away). By default, only the pressure channel is used, filtered from 3--8 s period. In general, higher frequencies are favorable.

It should be noted that drift estimates produced by default here may be improved upon by (1) incorporating all available channels, (2) utilizing multiple reference station pairings per estimate, (3) averaging over multiple frequency bands.

----

**a1_mseed2sac_rmresp.py** : Preprocessing script that converts miniseed data to SAC format. The script also removes instrument response to displacement and provides the option to downsample.

**a2_ccf_ambnoise.m** : Performs ambient-noise cross-correlation functions (CCF) for specified station pairs. By default, only the pressure channel is considered. Based on the [MATnoise](https://github.com/jbrussell/MATnoise) package.

**a3_measure_ccf_drfit.m** : Cross-correlate daily CCFs with reference stack in order to determine systematic clock drift rate.

**a4_mseed_removedrfit.py** : Remove clock drift (i.e., correct the timing) for raw miniseed files.

---
###### Optional
**b1_writedrifts2txt.m** : Write summary text file containing clock drifts for each station.


