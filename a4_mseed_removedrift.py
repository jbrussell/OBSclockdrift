# %% markdown
#
# Load mseed data, correct clock drift, and save as mseed.
# To correct, must *subtract* the accumulated drift from the daily start time.
# 
# jbrussell
# %% codecell
from obspy import read
from obspy.io.sac import SACTrace
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import os
import glob
import pandas as pd
from scipy.io import loadmat
# %% codecell
# Setup paths

# 3 stations with invalid drift measurements
path2mseed = "/Volumes/Russell_2TB/irma6/gaherty/youngORCA/OBS/NoClock/"
path2mseed_out = "/Volumes/Russell_2TB/irma6/gaherty/youngORCA/OBS/NoClock_corrected/"
path2sta = "/Volumes/Russell_2TB/CHAMBO/RESEARCH/PROJ_YoungPacificORCA/DATA/stations_all.txt" # station file
path2drift = "./clock_drifts/"

stas1 = ['WC01','EC01'] # station without drift correction
stas2 = ['CC04','CC04'] # reference station (with correct timing)

coperiod = [3, 8]
strNAMEcomp = 'PP'

# Channels to correct drift for
# chs = ["CH3"]
chs = ["CH0", "CH1", "CH2", "CH3"]

# %% codecell
# Load station file
inventory = pd.read_csv(path2sta, delimiter= '\s+', index_col=False, names=['station','stla','stlo','stel'])
inventory = inventory.set_index('station')
inventory = inventory.loc[stas1].reset_index()
# %% codecell

# Loop over stations
for ista,sta in enumerate(inventory.station):
    if not os.path.exists(path2mseed_out+sta):
        os.makedirs(path2mseed_out+sta)
        
    path2mseedsta = path2mseed + sta + '/Mseed/'
    # Loop over channels
    for ch in chs:
        # Get component name (flip BH2 and BH1)
        if ch == "CH0":
            comp = "BH2"
        elif ch == "CH1":
            comp = "BH1"
        elif ch == "CH2":
            comp = "BHZ"
        elif ch == "CH3":
            comp = "BDH"
        pathlist = sorted(Path(path2mseedsta).glob('**/*'+ch+'*.msd'))
        
        # Loop over day files and extract days
        for ifil,path in enumerate(pathlist):   
            # Load data
            st = read(str(pathlist[ifil]))
            
            if ifil == 0:
                tstart_deployment = st[0].stats.starttime
            
            # Check if day already processed
            if len(glob.glob(path2mseed_out+sta+'/'+sta+'.'+str(st[0].stats.starttime.year)+'.'+'%03i'%(st[0].stats.starttime.julday)+'.*.'+comp+'.sac'))!=0:  
                print('Skipping... already processed  ' + str(path))
                continue
            print("working on: " + str(path))
            
            # Read matlab file for drift rate
            filename = path2drift+'clockdrift'+strNAMEcomp+'_'+sta+'_'+stas2[ista]+'_'+str(coperiod[0])+'_'+str(coperiod[1])+'s.mat'
            temp = loadmat(filename)
            secperday = temp['secperday'][0][0]
            
            # Correct start time for drift
            tstart_day = st[0].stats.starttime
            tdrift_sec = (tstart_day-tstart_deployment)/60/60/24 * secperday
            st[0].stats.starttime = tstart_day - tdrift_sec
            
            sta_str = sta+'.'+ch+st[0].stats.starttime.strftime(".%Y.%j.%H.%M.%S.msd")
            st.write(path2mseed_out+sta+'/'+sta_str, format="MSEED")  

