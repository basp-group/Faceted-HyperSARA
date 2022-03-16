# Requirements
1. Casacore https://github.com/casacore/casacore
2. Meqtrees https://github.com/ratt-ru/meqtrees/wiki/Installation

# General notes
1. The user must provide the name/path to the measurement set. The remaining optional inputs are:
```bash
"OUT" # nametag of the source of interest, it will define the directory of the extracted data, default OUT=""
"FIELDID" # field ID of the source of interest, default FIELDID=0
"CHF" # ID of the first channnel to be extracted, default CHF=0
"CHL" # ID of the last channnel to be extracted, default  CHL=NUM_CHAN.
```
2. Extracted .mat files are saved in `./data/`. The files contain the following fields:
```bash
"frequency" # channel frequency                       
"y"  # data (Stokes I)
"u"  # u coordinate (in units of the wavelength)
"v"  # v coordinate (in units of the wavelength)
"w"  # w coordinate (in units of the wavelength)                       
"nW"  # 1/sigma
"nWimag" # imaging weights if available (Briggs or uniform), empty otherwise
"maxProjBaseline" # max projected baseline (in units of the wavelength)
```
3. Channels indices and hence their associated filenames are incremented by 1 for MATLAB indexing purposes.
4. Functions in `pyxis_ckat.py` are examples. The user is encouraged to tailor these functions if the measurement set is organised differently.
5. In case of multiple spectral windows, the user is encouraged to use the  `split` fonctionality in CASA to image specific spectral windows and extract the data from the resulting subset measurement set.

# Examples
#### 1. Single measurement set (MS)

Extracting data from the channels indexed from  `0` to `12` of the source with field ID `0`.
```bash
pyxis MS=myms.ms  OUT=mysource_nametag CHF=0  CHL=12 FIELDID=0 getdata_ms
```
Data will be saved in the directory `./data/mysource_nametag/`.

#### 2. Combine two measurement sets at two consecutive frequency bandwidths, with same spectral window specs
Extracting data from the channels indexed from  `0` to `15` of the source with field ID `2`, from **all** spectral windows.
```bash
pyxis MSLOW=mylowband_ms.ms MSHIGH=myhighband_ms.ms CHF=0 CHL=15 FIELDID=2 OUT=mysource_nametag getdata_ms_concat_bandwidth
```
 Data will be saved in the directory `./data/mysource_nametag/`
#### 3. Combine different measurment sets spanning the same frequency bandwidth, with same spectral window specs
Extracting data from the channels indexed from  `0` to `15` of the source with field ID `0`, from **all** spectral windows.

```bash
pyxis MS=myms1.ms  OUT=mysource_nametag1 CHF=0  CHL=12 FIELDID=0 getdata_ms
pyxis MS=myms2.ms  OUT=mysource_nametag2 CHF=0  CHL=12 FIELDID=0 getdata_ms
pyxis MS=mymsn.ms  OUT=mysource_nametagn CHF=0  CHL=12 FIELDID=0 getdata_ms
```
Data sets will be saved in seperate directories. 
```bash
 /data/mysource_nametag1/
 /data/mysource_nametag2/
 .
 .
 /data/mysource_nametagn/
```
To combine the data sets during imaging, the user should provide the nametags of the different sets in a cell as input.
```bash
datasetNames = {'mysource_nametag1', 'mysource_nametag2',..,{mysource_nametagn}};
```
