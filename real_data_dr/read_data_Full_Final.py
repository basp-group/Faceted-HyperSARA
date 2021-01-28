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
from pyrap.tables import table as tbl
import pyfits
import Tigger
import numpy
import os
import math
import json
import time

PI = math.pi
FWHM = math.sqrt( math.log(256) )
C = 2.998e8

msfreq = ms.ms("CYG-A-7-8M10S.MS::SPECTRAL_WINDOW",write=True)
print msfreq.colnames()

Freqs = msfreq.getcol("CHAN_FREQ")
print Freqs.shape
RefFreq = msfreq.getcol("REF_FREQUENCY")

info("Freq loaded !!!!!!")
    
msfreq.close()

    
tab = ms.ms("CYG-A-7-8M10S.MS",write=True)
print tab.colnames()

uvw = tab.getcol("UVW")
info("UVW !!!!!!")
info(uvw.shape)

uvw_new = uvw[:,[0,1]];

info(uvw_new.shape)
    
vis = tab.getcol("DATA")
info("DATA !!!!!!")
info(vis.shape)

vis_new = 0.5*(vis[:,:,0]+vis[:,:,1]);

info(vis_new.shape)

Flags = tab.getcol("FLAG")
info("FLAGS !!!!!!")
info(Flags.shape)

Flags_new  = numpy.logical_or(Flags[:,:,0],Flags[:,:,1]);

info(Flags_new.shape)

flag_row = tab.getcol("FLAG_ROW")
info("FLAG_ROW !!!!!!")
info(flag_row.shape)

weights_ch = tab.getcol("WEIGHT_SPECTRUM")
info("WEIGHT_SPECTRUM !!!!!!")
info(weights_ch.shape)

weights_ch_new = numpy.sqrt(4./((1./weights_ch[:,:,0])+(1./weights_ch[:,:,1])));

info(weights_ch_new.shape)
    
times = tab.getcol("TIME")
info("TIME !!!!!!")
info(times.shape)
    
data_id = tab.getcol("DATA_DESC_ID")
info("DATA_DESC_ID !!!!!!")
info(data_id.shape)
    
field_id = tab.getcol("FIELD_ID")
info("FIELD_ID !!!!!!")
info(field_id.shape)
    
tab.close()

sio.savemat('CYG-A-7-FULL.mat', {'vis':vis_new, 'flaging':Flags_new, 'weights_ch':weights_ch_new, 'Freqs':Freqs, 'RefFreq':RefFreq, 'uvw':uvw, 'flag_row':flag_row, 'times':times, 'data_id':data_id, 'field_id':field_id})
