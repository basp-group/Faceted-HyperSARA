# Requirements
1. Casacore https://github.com/casacore/casacore
2. Meqtrees https://github.com/ratt-ru/meqtrees/wiki/Installation

# General notes
1. Data extraction is performed using `pyxis`. The task should be performed in `./Faceted-Hyper-SARA/pyxis_ms2mat` which contains the python script `pyxis_ckat.py`.
2. Extracted .mat files are saved in `../data/`. The files contain the following fields
```bash
"frequency" # channel frequency                       
"y"  # data (Stokes I)
"u"  # u coordinate (in units of the wavelength)
"v"  # v coordinate (in units of the wavelength)
"w"  # w coordinate (in units of the wavelength)                       
"nW"  # sqrt(weights)
"nWimag" # imaging weights if available (Briggs or uniform), empty otherwise
"maxProjBaseline" # max projected baseline (in units of the wavelength)
```
3. Channels' indices and their associated filenames are incremented by 1 for MATLAB indexing purposes.
4. Functions in `pyxis_ckat.py` are examples. The user is encouraged to tailor these functions if the measurement set is organised differently.
5. In case of multiple spectral windows, the user is advised to use the  `split` task in CASA to image specific spectral windows and extract the data from the resulting measurement subset.

# Examples
#### 1. Single measurement set (MS)

Extracting data from the channels indexed from  `0` to `3` of the source with field ID `0`, from **all** spectral windows.

The user must provide the name/path to the measurement set `$MS`. The following inputs are optional.
```bash
"SRCNAME" # default SRCNAME="" : source nametag which defines the main directory of the extracted data.
"MSTAG" # default MSTAG="" : nametag of the data set which defines the subdirectory of the extracted data. It needs to be set when multiple data sets of the source of interest are available. 
"FIELDID" # default FIELDID=0 : field ID of the source of interest, 
"FIRSTCH" # default FIRSTCH=0 : ID of the first channnel to be extracted, 
"LASTCH" # default  LASTCH=NUM_CHAN : ID of the last channnel to be extracted
```
From the terminal launch:
```bash
pyxis MS=myms.ms  SRCNAME=cyga FIRSTCH=0  LASTCH=3 FIELDID=0 getdata_ms
```

Data will be saved as .mat files in the directory `../data/mysource_nametag/`. In the case of a single spectral window, the outcome is as follows
```bash
../data/cyga/data_ch_1.mat
.
.
../data/cyga/data_ch_4.mat
```
In `imaging/main_input_imaging.m`, if the data set nametag (`$MSTAG`) was not provided during data extraction, as in the example above, the user should set `datasetNames` as
```bash
datasetNames = {''};
```
#### 2. Combine two measurement sets at two consecutive frequency bandwidths, with same spectral window specs
Extracting data from the channels indexed from  `0` to `15` of the source with field ID `2`, from **all** spectral windows.

The user must provide the name/path to the two measurement sets `$MSLOW` and `$MSHIGH`. The following inputs are optional.
```bash
"SRCNAME" # default SRCNAME="" : source nametag which defines the main directory of the extracted data.
"MSTAG" # default MSTAG="" : nametag of the data set which defines the subdirectory of the extracted data. It needs to be set when multiple data sets of the source of interest are available. 
"FIELDID" # default FIELDID=0 : field ID of the source of interest, 
"FIRSTCH" # default FIRSTCH=0 : ID of the first channnel to be extracted, 
"LASTCH" # default  LASTCH=NUM_CHAN : ID of the last channnel to be extracted
```
From the terminal launch:
```bash
pyxis MSLOW=mylowbandms.ms MSHIGH=myhighbandms.ms SRCNAME=cyga FIRSTCH=0 LASTCH=15 FIELDID=2  getdata_ms_concat_bandwidth
```
Data will be saved as .mat files in the directory `../data/cyga/`.

In the case of a single spectral window, the outcome is as follows
```bash
../data/cyga/data_ch_1.mat  % channel 1 from mylowband_ms.ms
.
../data/cyga/data_ch_16.mat % channel 16 from mylowband_ms.ms
../data/cyga/data_ch_17.mat % channel 1 from myhighband_ms.ms
.
../data/cyga/data_ch_32.mat % channel 16 from myhighband_ms.ms
```
#### 3. Combine different measurment sets spanning the same frequency bandwidth, with same spectral window specs
Extracting data from the channels indexed from  `0` to `15` of the source with field ID `0`, from **all** spectral windows.
From the terminal launch:
```bash
pyxis MS=myms1.ms  SRCNAME=cyga MSTAG=dataset1 FIRSTCH=0  LASTCH=12 FIELDID=0 getdata_ms \
pyxis MS=myms2.ms  SRCNAME=cyga MSTAG=dataset2 FIRSTCH=0  LASTCH=12 FIELDID=0 getdata_ms \
pyxis MS=mymsn.ms  SRCNAME=cyga MSTAG=datasetn FIRSTCH=0  LASTCH=12 FIELDID=0 getdata_ms \
```
Data sets will be saved in the sub-directories. 
```bash
 ../data/cyga/dataset1/
 ../data/cyga/dataset2/
 ../data/cyga/datasetn/
```
In `imaging/main_input_imaging.m`, the user should provide the nametags of the different sets in a cell as input to combine the data sets during imaging.
```bash
datasetNames = {'dataset1', 'dataset2','datasetn'};
```

