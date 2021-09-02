# Dev notes

## Installation

**Installation**: To properly clone the project with the submodules, you may need to do follow one of set of instructions:

- updating from an existing `Faceted-Hyper-SARA` repository:

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
git submodule update --remote --merge # fetch and merge latest state of the submodule
```

If you get an access error, fix the `.gitmodules` file as follows

```bash
[submodule "lib/faceted-wavelet-transform"]
        path = lib/faceted-wavelet-transform
        url = https://github.com/basp-group-private/faceted-wavelet-transform.git
[submodule "lib/measurement-operator"]
        path = lib/measurement-operator
        url = https://github.com/basp-group-private/measurement-operator.git
```

- cloning the repository from scratch

```bash
git clone --recurse-submodules https://github.com/basp-group-private/Faceted-Hyper-SARA.git
```

In the long term, the `Faceted-Hyper-SARA` project might be updated to use `subtrees` instead of `submodules` (though it seems to require a recent version of git (not installed by default on `CIRRUS`)).

## Note: reading selected channels from a large `.fits` file in Matlab

```matlab
% extract first and last channels from the fits file
info        = fitsinfo('filename.fits');
rowend      = info.PrimaryData.Size(1);
colend      = info.PrimaryData.Size(2);
sliceend    = info.PrimaryData.Size(3);
data = fitsread('filename.fits','primary',...
          'Info', info,...
          'PixelRegion',{[1 1 rowend], [1 1 colend], [1 sliceend-1 sliceend]});
```

## Code formatting

Try the [`miss_hit`](https://github.com/florianschanda/miss_hit) package (to be checked first on a few examples before moving to the full project).

```bash
pip install miss_hit
mh_style my_file.m # only analyse the file
mh_style --fix my_file.m # fix the file based on the rules given in miss_hit.cfg
mh_style folder/
mh_style --fix .
```

## New interface (generic measurement  operator, lambda function)

Functions to be modified:

- `sara`: l. 172-181, 311-333, 442-452, 527-536
- solvers `hs`, `fhs`: update call to the 2 functions below
- `src/`: `compute_residual_images`, `update_dual_fidelity`
