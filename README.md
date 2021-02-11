# Parallel faceted imaging in radio interferometry via proximal splitting

**Description:** Matlab codes associated with the method described in

>P.-A. Thouvenin, A. Abdulaziz, M. Jiang, A. Dabbech, A. Repetti, A. Jackson, J.-P. Thiran, Y. Wiaux, Parallel faceted imaging in radio interferometry via proximal splitting (Faceted HyperSARA), submitted, [preprint available online](https://arxiv.org/abs/2003.07358), Jan. 2020.  

**Authors:** P.-A. Thouvenin, A. Abdulaziz, M. Jiang, A. Dabbech.

**Documentation:** a full documentation is currently under development, and will be made available in the present repository.

**Dependencies:** the present codes depend on the content of the `measurement-operator` github repository, loaded as a github `submodule`. This module contains codes associated with the following publications

> J. A. Fessler and B. P. Sutton, Nonuniform Fast Fourier Transforms Using Min-Max Interpolation, *IEEE Trans. Image Process.*, vol. 51, n. 2, pp. 560--574, Feb. 2003.
>
> A. Dabbech, L. Wolz, L. Pratley, J. D. McEwen and Y. Wiaux, [The w-effect in interferometric imaging: from a fast sparse measurement operator to superresolution](http://dx.doi.org/10.1093/mnras/stx1775), *Mon. Not. Roy. Astron. Soc.*, 471(4):4300-4313, 2017.
> 
> A. Onose, A. Dabbech and Y. Wiaux, [An accelerated splitting algorithm for radio-interferometric imaging: when natural and uniform weighting meet](http://dx.doi.org/10.1093/mnras/stx755), *Mon. Not. Roy. Astron. Soc.*, 469(1):938-949, 2017.

**Installation**: To properly clone the project with the submodules, you may need to do follow one of set of instructions:

- updating from an existing `Faceted-Hyper-SARA` repository:

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
# git submodule update --remote # if only a submodule has changed
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
