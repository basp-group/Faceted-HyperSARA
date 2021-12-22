## Rodrigues compatibale simulation pipeline 
## sphesihle makhathini [sphemakh@gmail.com]

import Pyxis
import ms
import mqt 
import im
import imager
import lsm
import im.argo
import im.lwimager
from Pyxis.ModSupport import *
from astropy.io import fits
from astropy.io import ascii
from casacore.tables import *
import scipy.io as sio
import pyrap.tables
import pyfits
import Tigger
import numpy as np
import os
import math
import json
import time


PI = math.pi
FWHM = math.sqrt( math.log(256) )
C = 299792458

def getdata_ms_concat_bandwidth(msf_bandwidth1='$MSBW1',msf_bandwidth2='$MSBW2',freqFirst='$CHF',freqLast='$CHL',srcid='$FIELDID',mstag='$OUTPUT'):
    
    msname1 = II('%s'%msf_bandwidth1)
    msname2 = II('%s'%msf_bandwidth2)
    mstag=II('%s'%mstag)
    freqFirst = int(II('%s'%freqFirst))
    freqLast  = int(II('%s'%freqLast))
    srcid = int(II('%s'%srcid))
    fileid0=0; 
    fileidMid=getdata(msname=msname1,mstag=mstag,freqFirst=freqFirst,freqLast=freqLast,srcid=srcid,fileid0=fileid0);
    fileidEnd=getdata(msname=msname2,mstag=mstag,freqFirst=freqFirst,freqLast=freqLast,srcid=srcid,fileid0=fileidMid);
    info(" %s files extracted."%fileidEnd) 

def get_data_ms(msf='$MS',mstag='$OUTPUT',freqFirst='$CHF',freqLast='$CHL',srcid='$FIELDID'):
    msname = II('%s'%msf)
    mstag=II('%s'%mstag)
    freqFirst = int(II('%s'%freqFirst))
    freqLast  = int(II('%s'%freqLast))
    srcid = int(II('%s'%srcid))
    fileid0=0;
    fileidEnd=getdata(msname=msname,mstag=mstag,freqFirst=freqFirst,freqLast=freqLast,srcid=srcid,fileid0=fileid0)


   
def getdata(msname,mstag,freqFirst,freqLast,srcid,fileid0):

    msname = II('%s'%msname)
    data_parent="./data"
    info("Data .mat files will be saved in %s"%data_parent)
    x.sh("mkdir -p  %s"%data_parent)
    x.sh("mkdir -p  %s/%s"%(data_parent,mstag))    
    
    
    path_data= '%s/%s/'%(data_parent,mstag)
    x.sh("mkdir -p  %s"%path_data)  
    
    info(msname)
    tab = ms.ms(msname,write=False);
    info("columns are",*(tab.colnames()));
    
    spwtab = ms.ms(msname,subtable="SPECTRAL_WINDOW");
    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,freqFirst] ;#get the first freq
    freqStepVect = spwtab.getcol("CHAN_WIDTH")
    freqsVect=spwtab.getcol("CHAN_FREQ")
    nspw =  len(freqsVect[0,:]);
    nSpw = len(freqsVect[:,1]);
    if nSpw >1:
        info("more than one spectral window detected.")
    if nSpw==1:
        info("First physical channel %s, last physical channel %s"%(freqFirst,freqLast)) 
    else:
        info("First SPW %s, last  SPW %s"%(freqFirst,freqLast))    

    bandwidth = spwtab.getcol("CHAN_WIDTH")[ms.SPWID,0]
    spwtab.close()
    # load remaining specs
    uvw = tab.getcol("UVW");
    nmeas =len(uvw[:,0]);
    info('Number of measurements per channel %s'%nmeas);

    field=tab.getcol("FIELD_ID");
    srcrows=field==srcid; # get src data
    field = field[srcrows];
    data_id=tab.getcol("DATA_DESC_ID");
    data_id=data_id[srcrows];
    ant1=tab.getcol("ANTENNA1");
    ant1= ant1[srcrows];
    ant2=  tab.getcol("ANTENNA2");
    ant2= ant2[srcrows] ;
    time=tab.getcol("TIME");
    time = time[srcrows];
    scantab= tab.getcol("SCAN_NUMBER")
    scantab= scantab[srcrows] ;
    integrationTimeVect=tab.getcol("EXPOSURE");
    integrationTimeVect= integrationTimeVect[srcrows]; 
    dt = tab.getcol("EXPOSURE",0,1)[0]
    dtf = (tab.getcol("TIME",tab.nrows()-1,1)-tab.getcol("TIME",0,1))[0]
    npts=len(ant1);
    info("Number of measurements per channel associated with the src of interest: %s"%npts)
    uvw=uvw[srcrows,:];
    #load flags
    flag_row = tab.getcol("FLAG_ROW")
    flag_row=flag_row==False
    flag_row = flag_row[srcrows];
    flag_row = flag_row.transpose()
    flag_row = flag_row.astype(float) #flag_row==False
    flagAll= tab.getcol("FLAG")
    flagAll=flagAll==False
    flagAll= flagAll.astype(float);
    flagAll= flagAll[srcrows,:,:]
    flag_doneAll = flagAll[:,:,0]+flagAll[:,:,3]

    flagAll=[];
    # load data
    try:
        data_full = tab.getcol("CORRECTED_DATA");
    except:
        data_full = tab.getcol("DATA");
        info("CORRECTED_DATA not found --> will use DATA ")
    data_full=data_full[srcrows,:,:];
    data_size = data_full.shape
    ncorr = data_size[2]
    info("%s available correlations"%ncorr)
    data_corr_1 =  data_full[:,:,0];
    data_corr_4 =  data_full[:,:,ncorr-1];
    data_full =[];
    
    if ncorr== 4:
       info("All correlations are available")
    else:
       info("Only Stokes I corr are available")

    #load weights
    isavail_weight_spectrum=0;
    try:
        weight_ch = tab.getcol("WEIGHT_SPECTRUM");
        weight_ch=weight_ch[srcrows,:,:];
        weight_corr_1 =  weight_ch[:,:,0];
        weight_corr_4 =  weight_ch[:,:,ncorr-1];
        weight_ch =[];
        isavail_weight_spectrum=1;

    except: 
        info('WEIGHT_SPECTRUM not available --> using WEIGHT')
        weight = tab.getcol("WEIGHT");
        weight =weight[srcrows];
        weight_corr_1 =  weight[:,0];
        weight_corr_4 =  weight[:,ncorr-1];
        weight =[];
    
    #load imaging weights
    isavail_imaging_weight =0
    try:
        weightImaging= tab.getcol("IMAGING_WEIGHT");
        weightImaging = weightImaging[srcrows];
        isavail_imaging_weight=1;
        weightImaging_corr_1 = weightImaging[:,0]; # typically same weights for all corr
        weightImaging=[];
        info('IMAGING_WEIGHTS is  available')
    except:
        try:
            weightImaging= tab.getcol("IMAGING_WEIGHT_SPECTRUM");
            weightImaging=weightImaging[srcrows,:,:];
            isavail_imaging_weight=1;
            weightImaging_corr_1 = weightImaging[:,:,0]; # typically same weights for all corr
            weightImaging=[];
            info('IMAGING_WEIGHTS_SPECTRUM is  available')
        except:  
            weightImaging=[];
   
   
    for ifreq in range(freqFirst, freqLast+1):
        #imaging weights
        wimagI = ''
        if isavail_imaging_weight>0:
            try:
                wimagI = weightImaging_corr_1[:,ifreq] ;
            except:
                wimagI = weightImaging_corr_1;
        #natural weights
        if isavail_weight_spectrum>0:
            w1 = weight_corr_1[:,ifreq];
            w4 = weight_corr_4[:,ifreq];
        else:
            w1 = weight_corr_1;
            w4 = weight_corr_4;
        #Stokes I 
        data = (w1 * data_corr_1[:,ifreq] + w4 * data_corr_4[:,ifreq]) / (w1+w4);
        # data = (data_corr_1[:,ifreq] +data_corr_4[:,ifreq])*0.5
        weightI = w1 +w4
        # flag
        flag=flag_row+flag_doneAll[:,ifreq]+ numpy.absolute(data); 
        flag=flag!=0 #numpy.array(flag0,dtype=bool);
        # save data
        info("Reading data and writing file..Freq %s"%(ifreq))
        
        if nSpw>1:    
            for ich in range(0,nSpw):
               rows_slice=data_id==ich ;
               data_slice = data[rows_slice];
               weightI_slice = weightI[rows_slice];
               flag_slice =flag[rows_slice];
               try: 
                   wimagI_slice=wimagI[rows_slice];
                   info("Imaging weights available")
               except:
                   wimagI_slice=[];
               frequency=freqsVect[ich,ifreq]
               info("Current freq %s of spectral window %s: %s MHz"%(ifreq,ich,frequency))
               uvw_slice=uvw[rows_slice]/(C/frequency);
               u = uvw_slice[:,0];
               v = uvw_slice[:,1];
               w = uvw_slice[:,2];
               uvw_slice = [];
               flag_slice =flag[rows_slice];
               y=data_slice[flag_slice];
               data_slice=[];
               u=u[flag_slice];
               info(u.shape)
               v=v[flag_slice];
               w=w[flag_slice];
               maxProjBaseline=numpy.sqrt(max(u**2+v**2));
               nW=numpy.sqrt(weightI_slice[flag_slice]);
               try:
                  wimag=numpy.sqrt(wimagI_slice[flag_slice]);
               except:
                  wimag=[];
               fileid =1+fileid0+ (ich*16)+ifreq;# fileid0+1;
               dataFileName="%s/data_ch_%s.mat"%(path_data,fileid)
               info('Channel id: %d. File saved at %s'%(fileid,dataFileName))
               sio.savemat(dataFileName,{'maxProjBaseline':maxProjBaseline,'u':u,'v':v,'w':w,'y':y,'nW':nW,'nWimag':wimag,'frequency':frequency}); 

        else:
            frequency=freqsVect[ifreq];
            y=data[flag];
            data=[];
            u=uvw[flag,0];
            v=uvw[flag,1];
            w=uvw[flag,2];
            uvw=[];
            nW=numpy.sqrt(weightI[flag]);
            wimagI=numpy.sqrt(wimagI[flag]);
            maxProjBaseline=numpy.sqrt(max(u**2+v**2));
            dataFileName="%s/data_ch_%s.mat"%(path_data,ifreq+1)
            info('Channel id: %d. File saved at %s'%(ifreq+1,dataFileName))
            sio.savemat(dataFileName,{'maxProjBaseline':maxProjBaseline,'u':u,'v':v,'w':w,'y':y,'nW':nW,'nWimag':wimag,'frequency':frequency});
            fileid =ifreq+1;
            info("Current freq %s"%ifreq)
        
    tab.close() 
    fileid0 = fileid;
    info('Last file id %s'%fileid0)
    return fileid0;
