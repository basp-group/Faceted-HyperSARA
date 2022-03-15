*Requirements:
Casacore https://github.com/casacore/casacore
Meqtrees https://github.com/ratt-ru/meqtrees/wiki/Installation

*General notes:
1. Output .mat files contain the following fields:
"frequency": frequency, # channel frequency in MHz                      
"y": y,# data (Stokes I)
"u": u,# u coordinate (in units of the wavelength)
"v": v,# u coordinate (in units of the wavelength)
"w": w,# u coordinate (in units of the wavelength)                       
"nW": nW,# 1/sigma
"nWimag": nWimag, #imaging weights if available (Briggs or uniform), empty otherwise
"maxProjBaseline": maxProjBaseline,#max projected baseline (in units of the wavelength)

2. Functions in pyxis_ckat.py are examples. Users are encouraged to tailor the functions if their measurement set is organised differently.

*Examples:
1. Single measurement set (MS)
Extracting data from the channels indexed from  0 to 12 of the source with field ID 0.
Data will be saved in the directory "./data/mysource_name/"

"pyxis MS=myms.ms  OUT=mysource_name CHF=0  CHL=12 FIELDID=0 get_data_ms"

2. Combine two measurement sets at two consecutive frequency bandwidths, with same spectral window specs
Extracting data from the channels indexed from  0 to 15 of the source with field ID 2, from ALL spectral windows.
Data will be saved in the directory "./data/mysource_name/"

"pyxis MSLOW=mylowband_ms.ms MSHIGH=myhighband_ms.ms CHF=0 CHL=15 FIELDID=2 OUT=mysource_name getdata_ms_concat_bandwidth"



