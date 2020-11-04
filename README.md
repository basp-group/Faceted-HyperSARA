# Parallel faceted imaging in radio interferometry via proximal splitting

**Description:** Matlab codes associated with the method described in 

>P.-A. Thouvenin, A. Abdulaziz, M. Jiang, A. Dabbech, A. Repetti, A. Jackson, J.-P. Thiran, Y. Wiaux -
<strong>Parallel faceted imaging in radio interferometry via proximal splitting (Faceted HyperSARA)</strong>, submitted, <a href="https://researchportal.hw.ac.uk/en/publications/parallel-faceted-imaging-in-radio-interferometry-via-proximal-spl">preprint available online</a>, Jan. 2020.  

**Authors:** P.-A. Thouvenin (pierre-antoine.thouvenin@centralelille.fr), A. Abdulaziz (aa61@hw.ac.uk), M. Jiang (ming.jiang@epfl.ch).

**Documentation:** a full documentation is currently under development, and will be made available in the present repository.

**Dependencies:** the present codes includes a slightly modified version of the MATLAB NUFFT algorithm available at http://web.eecs.umich.edu/~fessler/irt/fessler.tgz, described in

> J. A. Fessler and B. P. Sutton - 
<strong>Nonuniform Fast Fourier Transforms Using Min-Max Interpolation</strong>, <em>IEEE Trans. Image Process.</em>, vol. 51, n. 2, pp. 560--574, Feb. 2003.

**Installation** To properly clone the project with the submodules, you may need to do follow one of set of instructions:

- updating from an existing `Faceted-Hyper-SARA` repository:

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
```

If you get an access error, fix the `.gitmodules` file as follows

```bash
[submodule "lib/faceted-wavelet-transform"]
        path = lib/faceted-wavelet-transform
        url = https://github.com/basp-group-private/faceted-wavelet-transform.git
```

- cloning the repository from scratch

```bash
git clone --recurse-submodules https://github.com/basp-group-private/Faceted-Hyper-SARA.git
```

In the long term, the `Faceted-Hyper-SARA` project might be updated to use `subtrees` instead of `submodules` (though it seems to require a recent version of git (not installed by default on `CIRRUS`)).
