# Parallel faceted imaging in radio interferometry via proximal splitting

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

- submodules update: updating from an existing `Faceted-Hyper-SARA` repository:

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
git submodule update --remote --merge # fetch and merge latest state of the submodule
```

## Configuration

To reproduce the experiments (on CIRRUS), configure the `.csv` file contained in ..., and run the following

```bash
cd path/to/Faceted-Hyper-SARA
module load anaconda/python3
# configure / update the python job script
vi job_spatial.py
python job_spatial.py
```

To use the scripts reported in `figures`, make sure the [`Cubehelix`](...) package is installed.

## Real data development

Most of the instructions reported below refer to the content of the `experiment/real` folder. To help you in configuring the final scripts, you can use the scripts used for the first submission of the paper (`experiments/archive/mnras_real`) or those proposed by Arwa at an earlier stgae of the development (`experiments/archive/real_data`).

- All the source files for the solvers are defined in `src/`. Any modification to the measurement operator needs to be defined separately in the `measurement-operator` repository, later synchronized fron this project (see the submodules update instructions in the [installation](#installation) section).
- All user-defined and default parameters are defined in two separate files:
    1. `parameters_problem.m`: problem-specific parameters (blocking, nufft, ...)
    2. `parameters_solver.m`: solver-specific parameters (pdfb, reweighting, ...)
Note that real-data parameters, (weighting scheme, ...) if any, still need to be added in `parameters_problem.m`.
- Before running: update and configure the sections delimited by `# !` in following files
    1. `main.m`
    2. `job_cyga.py`
    3. `run_cyga.slurm`

## TODO

- [x] Prepare draft Python scripts for real data examples (only keep useful values, see if anything needs to be added)
- [ ] Replace FB by FISTA for the projection onto the ellipsoids
- [ ] Update function interface + name of variables (if necessary)
- [ ] Adding H matrices to `measurement-operator`
- [ ] Adapt synth data scripts to real data
- [ ] Documenting all functions (ongoing)
