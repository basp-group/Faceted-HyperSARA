# Parallel faceted imaging in radio interferometry via proximal splitting

![language](https://img.shields.io/badge/language-MATLAB-orange.svg)
[![license](https://img.shields.io/badge/license-GPL--3.0-brightgreen.svg)](LICENSE)
[![docs-page](https://img.shields.io/badge/docs-latest-blue)](https://basp-group.github.io/Faceted-HyperSARA/index.html)

- [Parallel faceted imaging in radio interferometry via proximal splitting](#parallel-faceted-imaging-in-radio-interferometry-via-proximal-splitting)
  - [Description](#description)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
    - [Cloning the project](#cloning-the-project)
    - [Generating a synthetic dataset to test the library](#generating-a-synthetic-dataset-to-test-the-library)
  - [Documentation](#documentation)
  - [Note to contributors](#note-to-contributors)
    - [Pushing the documentation online](#pushing-the-documentation-online)
    - [Pull request: Matlab code formatting](#pull-request-matlab-code-formatting)

## Description

``Faceted-HyperSARA`` is a [fully documented]((https://basp-group.github.io/Faceted-HyperSARA/index.html)) MATLAB library for radio-interferometric wideband intensity imaging. The library offers a collection of utility functions and scripts from data extraction from an RI measurement set ``MS Table`` to the reconstruction of a wideband intensity image over the field of view and frequency range of interest. The library contains an implementation of algorithms from the "SARA" family (SARA,
HyperSARA and Faceted HyperSARA), and more specifically described in the following publications. 

>P.-A. Thouvenin, A. Abdulaziz, A. Dabbech, A. Repetti, Y. Wiaux, Parallel faceted imaging in radio interferometry via proximal splitting (Faceted HyperSARA): I. Algorithm and simulations, submitted, [preprint available online](https://arxiv.org/abs/2003.07358), March 2022.  
>
>P.-A. Thouvenin, A. Dabbech, M. Jiang, A. Abdulaziz, J.-P. Thiran, A. Jackson, Y. Wiaux, Parallel faceted imaging in radio interferometry via proximal splitting (Faceted HyperSARA): II. Code and real data proof of concept, submitted, April  2022.

The following features, crucial to achieve high precision imaging from large data volumes, are supported:

- data dimensionality reduction via visibility gridding ([Kartik2017](https://academic.oup.com/mnras/article/468/2/2382/3061880));
- estimation of the effective noise level when reliable noise estimates are not available, due to imperfect calibration.
- correction of the `w`-term via `w`-projection ([Dabbech2017](https://academic.oup.com/mnras/article/471/4/4300/3965853));
- incorporation of available compact Fourier models of the direction dependent effects (DDEs) in the measurement operator ([Dabbech2021](https://academic.oup.com/mnras/article-abstract/506/4/4855/6315336?redirectedFrom=fulltext)).

## Dependencies 

``Faceted-HyperSARA`` relies on two auxiliary submodules (see [installation](#installation) section):

1. [`RI-measurement-operator`](https://github.com/basp-group/RI-measurement-operator) for the formation of the radio-interferometric measurement operator;
2. [`SARA-dictionary`](https://github.com/basp-group/SARA-dictionary) for the implementation of the sparsity priors.

These modules contains codes associated with the following publications

> A. Dabbech, L. Wolz, L. Pratley, J. D. McEwen and Y. Wiaux, [The w-effect in interferometric imaging: from a fast sparse measurement operator to superresolution](http://dx.doi.org/10.1093/mnras/stx1775), *Mon. Not. Roy. Astron. Soc.*, 471(4):4300-4313, 2017.
>
> J. A. Fessler and B. P. Sutton, Nonuniform Fast Fourier Transforms Using Min-Max Interpolation, *IEEE Trans. Image Process.*, vol. 51, n. 2, pp. 560--574, Feb. 2003.
>
> A. Onose, A. Dabbech and Y. Wiaux, [An accelerated splitting algorithm for radio-interferometric imaging: when natural and uniform weighting meet](http://dx.doi.org/10.1093/mnras/stx755), *Mon. Not. Roy. Astron. Soc.*, 469(1):938-949, 2017.
> 
> Z. Prusa, Segmentwise discrete wavelet transform, *Brno university of technology*, 2012.

## Installation

To properly clone the project with the submodules, you may need to choose one of following set of instructions.

### Cloning the project

- If you plan to clone the project from `https`, you will first need to edit the `.gitmodules` file accordingly, replacing the `git@github.com` addresses with the `https` counterpart. That is

```bash
[submodule "lib/SARA-dictionary"]
	path = lib/SARA-dictionary
	url = https://github.com/basp-group/SARA-dictionary.git
[submodule "lib/RI-measurement-operator"]
	path = lib/RI-measurement-operator
	url = https://github.com/basp-group/RI-measurement-operator.git
```

- Cloning the repository from scratch. If you used `https`, issue the following command

```bash
git clone --recurse-submodules https://github.com/basp-group/Faceted-HyperSARA.git
```

If you are using an SSH key for github rather than a personal token, then you will need to clone the repository as follows instead:

```bash
git clone --recurse-submodules git@github.com:basp-group/Faceted-HyperSARA.git
```

You will then also need to update the local repository configuration to use this approach for the sub-modules and update the submodules separately as detailed below. If you have problems updating the submodule access configuration open a issue on this repository.

- Submodules update: updating from an existing `Faceted-HyperSARA` repository.

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
git submodule update --remote --merge # fetch and merge latest state of the submodule
```

<!-- git submodule init -->
<!-- git submodule update --remote -->

### Generating a synthetic dataset to test the library

To generate the synthetic data cubes used in the spatial and spectral experiments, you will need to download the `S_DDE_MODEL.fits` image file associated with the following paper.

> A. Dabbech, A. Repetti, R. A. Perley, O. M. Smirnov, Y. Wiaux, [Cygnus A jointly calibrated and imaged via non-convex optimization from VLA data](https://doi.org/10.1093/mnras/stab1903), *Mon. Not. Roy. Astron. Soc.*, 506(4):4855-4876, Oct. 2021.

To do so, issue the following commands in a terminal.

```bash
# if on MAC: 
# brew install wget
cd path/to/Faceted-HyperSARA
cd data
wget -P . https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits
```

Please refer to the [online documentation](https://basp-group.github.io/Faceted-HyperSARA/index.html) to configure and run the imaging script (see [documentation](#documentation) section if you prefer building the documentation locally).

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

## Documentation

- Full library [documentation is accessible online](https://basp-group.github.io/Faceted-HyperSARA/index.html).

- To build the documentation on your machine, make sure the following Python packages have been installed, and issue the appropriate build command.

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
