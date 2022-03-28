# Note on the structure of the .mat files to be passed to the solver

- [Note on the structure of the .mat files to be passed to the solver](#note-on-the-structure-of-the-mat-files-to-be-passed-to-the-solver)
  - [List of input .mat files](#list-of-input-mat-files)
  - [Data (visibilities) files](#data-visibilities-files)
  - [Pre-processing](#pre-processing)
    - [$\ell_2$ bounds](#ell_2-bounds)
    - [Calibration](#calibration)
      - [DIE calibration](#die-calibration)
      - [DDE calibration](#dde-calibration)
    - [(Optional) Gridding matrices (NUFFT)](#optional-gridding-matrices-nufft)
    - [(Optional) Initial image file](#optional-initial-image-file)
  - [List of output .mat files](#list-of-output-mat-files)
    - [Operator norm](#operator-norm)
    - [Temporary results (checkpoints)](#temporary-results-checkpoints)

## List of input .mat files

One file per channel / effective channel. One dataset for each different VLA configuration (if relevant).

## Data (visibilities) files

`maxProjBaseline`
`u`
`v`
`w`
`nW`
`y`
`frequency`

## Pre-processing

### $\ell_2$ bounds

`sigmac`
`l2bounds`
`RESIDUAL_DATA`

### Calibration

#### DIE calibration

`DIES`

#### DDE calibration

`DDES`

### (Optional) Gridding matrices (NUFFT)

### (Optional) Initial image file

Image initialized to 0 by default if no input initial image is found.

## List of output .mat files

### Operator norm

### Temporary results (checkpoints)

Name structure: ...
Variables: ...
