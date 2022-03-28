# Parallel faceted imaging in radio interferometry via proximal splitting

[![license](https://img.shields.io/badge/license-GPL--3.0-brightgreen.svg)](LICENSE)
[![docs-page](https://img.shields.io/badge/docs-latest-blue)](https://basp-group.github.io/Faceted-HyperSARA/_imaging/imaging.html)

## Description

This repository contains Matlab codes associated with the wideband radio-interferometry imaging approach approach described in

>P.-A. Thouvenin, A. Abdulaziz, M. Jiang, A. Dabbech, A. Repetti, A. Jackson, J.-P. Thiran, Y. Wiaux, Parallel faceted imaging in radio interferometry via proximal splitting (Faceted HyperSARA), submitted, [preprint available online](https://arxiv.org/abs/2003.07358), Jan. 2020.  

**Authors:** P.-A. Thouvenin, A. Dabbech, M. Jiang, A. Abdulaziz.

**Dependencies:** the present codes depend on the content of the [`RI-measurement-operator`](https://github.com/basp-group/RI-measurement-operator) and [`SARA-dictionary`](https://github.com/basp-group/SARA-dictionary) github repositories, loaded as github `submodule`s (see [installation](install) section). These modules contains codes associated with the following publications

> A. Dabbech, L. Wolz, L. Pratley, J. D. McEwen and Y. Wiaux, [The w-effect in interferometric imaging: from a fast sparse measurement operator to superresolution](http://dx.doi.org/10.1093/mnras/stx1775), *Mon. Not. Roy. Astron. Soc.*, 471(4):4300-4313, 2017.
>
> J. A. Fessler and B. P. Sutton, Nonuniform Fast Fourier Transforms Using Min-Max Interpolation, *IEEE Trans. Image Process.*, vol. 51, n. 2, pp. 560--574, Feb. 2003.
>
> A. Onose, A. Dabbech and Y. Wiaux, [An accelerated splitting algorithm for radio-interferometric imaging: when natural and uniform weighting meet](http://dx.doi.org/10.1093/mnras/stx755), *Mon. Not. Roy. Astron. Soc.*, 469(1):938-949, 2017.
> 
> Z. Prusa, Segmentwise discrete wavelet transform, *Brno university of technology*, 2012.

## Installation

To properly clone the project with the submodules, you may need to choose one of following set of instructions:

- cloning the repository from scratch

```bash
git clone --recurse-submodules https://github.com/basp-group/Faceted-HyperSARA.git
```
If you are using an SSH key for github rather than a personal token, then you will need to clone the repository as follows instead:

```bash
git clone git@github.com:basp-group/Faceted-HyperSARA.git
```

You will then also need to update the local repository configuration to use this approach for the sub-modules and update the submodules separately as detailed below. If you have problems updating the submodule access configuration open a issue on this repository.

- submodules update: updating from an existing `Faceted-HyperSARA` repository:

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
git submodule update --remote --merge # fetch and merge latest state of the submodule
```

- to generate the synthetic data cubes used in the spatial and spectral experiments, you will need to download the `S_DDE_MODEL.fits` image file associated with the following paper.

> A. Dabbech, A. Repetti, R. A. Perley, O. M. Smirnov, Y. Wiaux, [Cygnus A jointly calibrated and imaged via non-convex optimization from VLA data](https://doi.org/10.1093/mnras/stab1903), *Mon. Not. Roy. Astron. Soc.*, 506(4):4855-4876, Oct. 2021.

To do so, issue the following commands in a terminal.

```bash
# if on MAC: 
# brew install wget
cd path/to/Faceted-HyperSARA
cd data
wget -P . https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits
```

Please refer to the [online documentation](https://basp-group.github.io/Faceted-HyperSARA/_imaging/imaging.html) to configure and run the imaging script (see [documentation](doc) section if you prefer building the documentation locally).

<!-- ## Configuration

To reproduce the experiments on the Cirrus HPC system ([https://www.cirrus.ac.uk](https://www.cirrus.ac.uk)), configure the `.csv` file contained in `imaging/spatial` or `imaging/spectral` as required and run the following

```bash
cd path/to/Faceted-HyperSARA
module load anaconda/python3
# configure / update the python job script
vi job_spatial.py
python job_spatial.py
```
Cirrus is configured with Matlab and python installed using the module system and uses the Slurm batch system. If your system is configured differently then the batch files (`run_simulation.slurm` in the `imaging/spatial` or `imaging/spectral` directories) will need to be altered to fit your system setup. -->

## Documentation <a name="doc"></a>

- A full [documentation](https://basp-group.github.io/Faceted-HyperSARA/_imaging/imaging.html) is directly accessible online.

- To build the documentation, make sure the following Python packages have been installed, and issue the appropriate buid command.

    ```bash
    # setup conda environment to build the documentation
    conda env create --name fhs-doc --file environment.yml

    ## or using conda/pip
    # conda create -n fhs-doc
    # conda activate fhs-doc
    # conda install pip
    # pip install -r requirement.txt

    # building the documentation in html format
    cd docs
    make html
    ```

All the generated ``.html`` files are contained in the ``docs/build`` folder.

If needed, you can delete the `conda` environment as follows

```bash
conda env remove -n fhs-doc
```

## Note to contributors

### Pushing the documentation online

Add a `worktree` from the `master` branch to be able to push the documentation, once built locally.

```bash
# make sure the folder html does not exist before running the command
git worktree add docs/build/html gh-pages
cd docs/build/html
git add .
git commit -m "Build documentation as of $(git log '--format=format:%H' master -1)"
git push origin gh-pages
# delete the worktree
cd ../
git worktree remove html
```

### Pull request: Matlab code formatting

Make sure codes in any pull request have been properly formatted with the [`miss_hit`](https://pypi.org/project/miss-hit/) package using the `miss_hit.cfg` file provided

```bash
# activate fhs-doc environment (see previous paragraph)
conda activate fhs-doc
# run the following command from the root of the package (where the miss_hit.cfg file is)
mh_style --fix .
```
