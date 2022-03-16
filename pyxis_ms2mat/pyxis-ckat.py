## Rodrigues compatibale simulation pipeline
## sphesihle makhathini [sphemakh@gmail.com]

import json
import math
import os
import time
import im
import im.argo
import im.lwimager
import imager
import lsm
import mqt
import ms
import numpy as np
import pyfits
import pyrap.tables
import Pyxis
import scipy.io as sio
import Tigger
from astropy.io import ascii, fits
from casacore.tables import *
from Pyxis.ModSupport import *

PI = math.pi
FWHM = math.sqrt(math.log(256))
C = 299792458


def getdata_ms_concat_bandwidth(
    msf_lowbandwidth="$MSLOW",
    msf_highbandwidth="$MSHIGH",
    srcname="$SRCNAME",
    srcid="$FIELDID",
    mstag="$MSTAG",
    freqFirst="$FIRSTCH",
    freqLast="$LASTCH",
):

    msname1 = II("%s" % msf_lowbandwidth)
    msname2 = II("%s" % msf_highbandwidth)
    srcname = II("%s" % srcname)
    mstag = II("%s" % mstag)
    try:
       freqFirst = int(II("%s" % freqFirst))
    except:
       freqFirst =0;
    try:
       freqLast = int(II("%s" % freqLast))
    except: 
       freqLast =[];
    try:
        srcid = int(II("%s" % srcid))
    except:
        srcid =0;
        
        
    fileid0 = 0 # channel ID counter
    fileidMid = getdata(
        msname=msname1,
        srcname=srcname,
        mstag=mstag,
        freqFirst=freqFirst,
        freqLast=freqLast,
        srcid=srcid,
        fileid0=fileid0,
    )
    fileidEnd = getdata(
        msname=msname2,
        srcname=srcname,
        mstag=mstag,
        freqFirst=freqFirst,
        freqLast=freqLast,
        srcid=srcid,
        fileid0=fileidMid,
    )
    info(" %s files extracted." % fileidEnd)


def getdata_ms( msf="$MS", mstag="$MSTAG", freqFirst="$FIRSTCH", freqLast="$LASTCH", srcid="$FIELDID",srcname="$SRCNAME"):
    msname = II("%s" % msf)
    mstag = II("%s" % mstag)
    srcname = II("%s" % srcname)
    try:
        freqFirst = int(II("%s" % freqFirst))
    except:
        freqFirst = []
    try:
        freqLast = int(II("%s" % freqLast))
    except:
        freqLast =[]
    try:
        srcid = int(II("%s" % srcid))
    except: 
        srcid =0
        info("Field ID is not provided, data of field ID 0 will be extracted")
        
    fileid0 = 0
    fileidEnd = getdata(
        msname=msname,
        srcname=srcname,
        mstag=mstag,
        freqFirst=freqFirst,
        freqLast=freqLast,
        srcid=srcid,
        fileid0=fileid0,
    )


def getdata(msname, mstag, freqFirst, freqLast, srcid, srcname, fileid0):
    
    data_parent_dir = "../data"
    x.sh("mkdir -p  %s" % data_parent_dir)
    data_dir = "%s/%s"  % (data_parent_dir, srcname)
    x.sh("mkdir -p  %s" % (data_dir))
    data_dir = "%s/%s"  % (data_dir, mstag)
    x.sh("mkdir -p  %s" % (data_dir))
    info("Data .mat files will be saved in %s" %data_dir)
    
    info(msname)
    tab = ms.ms(msname, write=False)
    info("MS table columns:", *(tab.colnames()))

    spwtab = ms.ms(msname, subtable="SPECTRAL_WINDOW")    
    bandwidth = spwtab.getcol("CHAN_WIDTH")[ms.SPWID, 0]
    freqStepVect = spwtab.getcol("CHAN_WIDTH")
    freqsVect = spwtab.getcol("CHAN_FREQ")
    spw_numch = spwtab.getcol("NUM_CHAN")
    
    nSpw = len(spw_numch)
    nChperSpw = spw_numch[0]
    
    if nSpw > 1:     
        info("%s spectral window detected, composed of %s channels each"%(nSpw,nChperSpw))
    else:
        info("%s spectral window, composed of %s channels" %(nSpw,nChperSpw))    
    spwtab.close()
    
    if not freqFirst :
           freqFirst =0;
           info("First channel ID  is not provided or set to 0. First channel extracted ID=%s"%freqFirst)
    if not freqLast : 
           freqLast = nChperSpw-1;
           info("Highest channel ID  is not provided. Last channel extracted ID=%s" %freqLast)
           
    # load remaining specs
    field = tab.getcol("FIELD_ID")
    srcrows = field == srcid
    field = field[srcrows]
    uvw = tab.getcol("UVW")
    nmeas = len(uvw[:, 0])
    info("Number of measurements per channel %s" % nmeas)
    uvw = uvw[srcrows, :]
    data_id = tab.getcol("DATA_DESC_ID")
    data_id = data_id[srcrows]
    #ant1 = tab.getcol("ANTENNA1")
    #ant1 = ant1[srcrows]
    #ant2 = tab.getcol("ANTENNA2")
    #ant2 = ant2[srcrows]
    #time = tab.getcol("TIME")
    #time = time[srcrows]
    #scantab = tab.getcol("SCAN_NUMBER")
    #scantab = scantab[srcrows]
    #integrationTimeVect = tab.getcol("EXPOSURE")
    #integrationTimeVect = integrationTimeVect[srcrows]
    #dt = tab.getcol("EXPOSURE", 0, 1)[0]
    #dtf = (tab.getcol("TIME", tab.nrows() - 1, 1) - tab.getcol("TIME", 0, 1))[0]
    flag_row = tab.getcol("FLAG_ROW")
    flag_row = flag_row == False
    flag_row = flag_row[srcrows]
    flag_row = flag_row.transpose()
    flag_row = flag_row.astype(float)  # flag_row==False
    npts = len(flag_row)
    info( "Number of measurements per channel associated with the src of interest: %s" %npts)

    # load data
    try:
        data_full = tab.getcol("CORRECTED_DATA")
    except:
        data_full = tab.getcol("DATA")
        info("CORRECTED_DATA not found, reading DATA instead ")
        
    data_full = data_full[srcrows, :, :]
    data_size = data_full.shape
    ncorr = data_size[2]
    if ncorr == 4:
        info("All four correlations are available")
    else:
        info("%s correlations available, Stokes I assumed" %ncorr)
        
    data_corr_1 = data_full[:, :, 0]
    data_corr_4 = data_full[:, :, ncorr - 1]
    data_full = []

    # load remaining flags
    flagAll = tab.getcol("FLAG")
    flagAll = flagAll == False
    flagAll = flagAll.astype(float)
    flagAll = flagAll[srcrows, :, :]
    flagAll = flagAll[:, :, 0] + flagAll[:, :, ncorr - 1]

    
    # load weights
    isavail_weight_spectrum = 0
    try:
        weight_ch = tab.getcol("WEIGHT_SPECTRUM")
        weight_ch = weight_ch[srcrows, :, :]
        weight_corr_1 = weight_ch[:, :, 0]
        weight_corr_4 = weight_ch[:, :, ncorr - 1]
        weight_ch = []
        isavail_weight_spectrum = 1

    except:
        info("WEIGHT_SPECTRUM not available --> using WEIGHT")
        weight = tab.getcol("WEIGHT")
        weight = weight[srcrows]
        weight_corr_1 = weight[:, 0]
        weight_corr_4 = weight[:, ncorr - 1]
        weight = []

    # load imaging weights
    isavail_imaging_weight = 0
    try:
        weightImaging = tab.getcol("IMAGING_WEIGHT")
        weightImaging = weightImaging[srcrows]
        isavail_imaging_weight = 1
        weightImaging_corr_1 = weightImaging[:, 0]
        # typically same weights for all corr
        weightImaging = []
        info("IMAGING_WEIGHTS is  available")
    except:
        try:
            weightImaging = tab.getcol("IMAGING_WEIGHT_SPECTRUM")
            weightImaging = weightImaging[srcrows, :, :]
            isavail_imaging_weight = 1
            weightImaging_corr_1 = weightImaging[:, :, 0]
            # typically same weights for all corr
            weightImaging = []
            info("IMAGING_WEIGHTS_SPECTRUM is  available")
        except:
            weightImaging = []

    for ifreq in range(freqFirst, freqLast + 1):
        # imaging weights
        wimagI = ""
        if isavail_imaging_weight > 0:
            try:
                wimagI = weightImaging_corr_1[:, ifreq]
            except:
                wimagI = weightImaging_corr_1
        # natural weights
        if isavail_weight_spectrum > 0:
            w1 = weight_corr_1[:, ifreq]
            w4 = weight_corr_4[:, ifreq]
        else:
            w1 = weight_corr_1
            w4 = weight_corr_4
        # Stokes I
        data = (w1 * data_corr_1[:, ifreq] + w4 * data_corr_4[:, ifreq]) / (w1 + w4)
        # data = (data_corr_1[:,ifreq] +data_corr_4[:,ifreq])*0.5
        weightI = w1 + w4
        # flag
        flag = flag_row + flagAll[:, ifreq] + np.absolute(data)
        flag = flag != 0  # np.array(flag0,dtype=bool);
        # save data
        info("Reading data and writing file..Freq %s" % (ifreq))

        if nSpw > 1:
            for iSpw in range(0, nSpw):
                frequency = freqsVect[iSpw, ifreq]
                #info("Spectral window %s, current freq id %s: %s MHz"%( iSpw, ifreq,frequency/1e6))
                rows_slice = data_id == iSpw
                flag_slice = flag[rows_slice]
                #data
                y = data[rows_slice]
                y = y[flag_slice]
                #1/sigma
                nW = weightI[rows_slice]
                nW = np.sqrt(nW[flag_slice])
                #imaging weights if available               
                try:
                    nWimag = wimagI[rows_slice]
                    nWimag =np.sqrt(nWimag[flag_slice])
                except:
                    nWimag= []
                #uvw in units of the wavelength
                uvw_slice = uvw[rows_slice] / (C / frequency)
                u = uvw_slice[flag_slice, 0]
                v = uvw_slice[flag_slice, 1]
                w = uvw_slice[flag_slice, 2]
                uvw_slice = []
                #max projected baseline (needed for pixelsize)
                maxProjBaseline = np.sqrt(max(u ** 2 + v ** 2))

                fileid = 1 + fileid0 + (iSpw * nChperSpw) + ifreq
                dataFileName = "%s/data_ch_%s.mat" % (data_dir, fileid)
                #info("Channel id: %d. File saved at %s" % (fileid, dataFileName))
                sio.savemat(
                    dataFileName,
                    {
                        "frequency": frequency,                      
                        "y": y,# data (Stokes I)
                        "u": u,# u coordinate (in units of the wavelength)
                        "v": v,# u coordinate (in units of the wavelength)
                        "w": w,# u coordinate (in units of the wavelength)                       
                        "nW": nW,# 1/sigma
                        "nWimag": nWimag, #imaging weights if available (Briggs or uniform) 
                        "maxProjBaseline": maxProjBaseline,#max projected baseline                         
                    },
                )
              
                info("Spectral window %s, current channel id %s: %s MHz. File saved at %s"%(iSpw,ifreq, frequency/1e6,dataFileName))
                

        else:
            # applying flags 
            frequency = freqsVect[0,ifreq]
            y = data[flag] ;  
            u = uvw[flag, 0]; 
            v = uvw[flag, 1]; 
            w = uvw[flag, 2];       
            nW = np.sqrt(weightI[flag]) 
            try:
                nWimag = np.sqrt(wimagI[flag])
            except:
                nWimag =[]
            maxProjBaseline = np.sqrt(max(u ** 2 + v ** 2)) 
            
            dataFileName = "%s/data_ch_%s.mat" % (data_dir, ifreq + 1)
            sio.savemat(
                dataFileName,
                {
                    "frequency": frequency,
                    "y": y,# data (Stokes I)
                    "u": u,# u coordinate (in units of the wavelength)
                    "v": v, # v coordinate (in units of the wavelength)
                    "w": w, # w coordinate  (in units of the wavelength)              
                    "nW": nW, # 1/sigma
                    "nWimag": nWimag,  #imaging weights if available (Briggs or uniform)         
                    "maxProjBaseline": maxProjBaseline,#max projected baseline 
                },
            )
            
            info("Spw 0, current channel id %s: %s MHz. File saved at %s"%(ifreq, frequency/1e6,dataFileName))

            fileid = ifreq + 1
            

    tab.close()
    fileid0 = fileid
    info("Last file id %s" % fileid0)
    return fileid0
