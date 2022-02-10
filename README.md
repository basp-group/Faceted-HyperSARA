# Parallel faceted imaging in radio interferometry via proximal splitting

[![license](https://img.shields.io/badge/license-GPL--3.0-brightgreen.svg)](LICENSE)
[![docs-page](https://img.shields.io/badge/docs-latest-blue)](#)
<!-- [![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit) -->

## Description

This repository contains Matlab codes associated with the approach described in

>P.-A. Thouvenin, A. Abdulaziz, M. Jiang, A. Dabbech, A. Repetti, A. Jackson, J.-P. Thiran, Y. Wiaux, Parallel faceted imaging in radio interferometry via proximal splitting (Faceted HyperSARA), submitted, [preprint available online](https://arxiv.org/abs/2003.07358), Jan. 2020.  

**Authors:** P.-A. Thouvenin, A. Abdulaziz, M. Jiang, A. Dabbech.

**Documentation:** a full documentation is currently under development, and will be made available in the present repository.

**Dependencies:** the present codes depend on the content of the `measurement-operator` github repository, loaded as a github `submodule`. This module contains codes associated with the following publications

> J. A. Fessler and B. P. Sutton, Nonuniform Fast Fourier Transforms Using Min-Max Interpolation, *IEEE Trans. Image Process.*, vol. 51, n. 2, pp. 560--574, Feb. 2003.
>
> A. Dabbech, L. Wolz, L. Pratley, J. D. McEwen and Y. Wiaux, [The w-effect in interferometric imaging: from a fast sparse measurement operator to superresolution](http://dx.doi.org/10.1093/mnras/stx1775), *Mon. Not. Roy. Astron. Soc.*, 471(4):4300-4313, 2017.
>
> A. Onose, A. Dabbech and Y. Wiaux, [An accelerated splitting algorithm for radio-interferometric imaging: when natural and uniform weighting meet](http://dx.doi.org/10.1093/mnras/stx755), *Mon. Not. Roy. Astron. Soc.*, 469(1):938-949, 2017.

## Installation

To properly clone the project with the submodules, you may need to choose one of following set of instructions:

- cloning the repository from scratch

```bash
git clone --recurse-submodules https://github.com/basp-group-private/Faceted-Hyper-SARA.git
```
If you are using an SSH key for github rather than a personal token, then you will need to clone the repository as follows instead:

```bash
git clone git@github.com:basp-group-private/Faceted-Hyper-SARA.git
```

You will then also need to update the local repository configuration to use this approach for the sub-modules and update the submodules separately as detailed below. If you have problems updating the submodule access configuration open a issue on this repository.

- submodules update: updating from an existing `Faceted-Hyper-SARA` repository:

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
cd path/to/Faceted-Hyper-SARA
cd data
wget -P . https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits
```

You can then execute the `experiments/sim_script_cyga_cubes.m` script to generate all the different synthetic data cubes considered in the paper.

## Configuration

To reproduce the experiments on the Cirrus HPC system ([https://www.cirrus.ac.uk](https://www.cirrus.ac.uk)), configure the `.csv` file contained in `imaging/spatial` or `imaging/spectral` as required and run the following

```bash
cd path/to/Faceted-Hyper-SARA
module load anaconda/python3
# configure / update the python job script
vi job_spatial.py
python job_spatial.py
```
Cirrus is configured with Matlab and python installed using the module system and uses the Slurm batch system. If your system is configured differently then the batch files (`run_simulation.slurm` in the `imaging/spatial` or `imaging/spectral` directories) will need to be altered to fit your system setup.


To use the scripts reported in `figures`, make sure the [`Cubehelix`](...) package is installed.

## Contributions

### Building the documentation

To build the documentation, make sure the following Python packages have been installed, and issue the appropriate buid command.

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

### Pushing the documentation online

Add a `worktree` from the `master` branch

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

### Matlab code formatting

Make sure any pull request has been properly formatted with the [`miss_hit`](https://pypi.org/project/miss-hit/) package using the `miss_hit.cfg` file provided

```bash
# activate fhs-doc environment (see previous paragraph)
conda activate fhs-doc
# run the following command from the root of the package (where the miss_hit.cfg file is)
mh_style --fix .
```
