Faceted-HyperSARA
==================

``Faceted-HyperSARA`` is a MATLAB library for radio-interferometric wideband intensity imaging. The library offers a collection of utility functions and scripts from data extraction from an RI measurement set ``MS Table`` to the reconstruction of a wideband intensity image over the field of view and frequency range of interest. The library contains an implementation of algorithms from the "SARA" family (SARA :cite:p:`Carrillo2012,Onose2017`,
HyperSARA :cite:p:`Abdulaziz2019` and Faceted HyperSARA
:cite:p:`Thouvenin2021`). The follwing features, crucial to achieve high precision imaging from large data volumes, are supported.

- data dimensionality reduction via visibility gridding :cite:p:`Kartik2017`;
- correction of the :math:`w`-term via :math:`w`-projection :cite:p:`Dabbech2017`;
- incorporation of available compact Fourier models of the direction dependent effects (DDEs) in the measurement operator :cite:p:`Dabbech2021`.

``Faceted-HyperSARA`` relies on two auxiliary submodules:

1. `RI-measurement-operator <https://github.com/basp-group/RI-measurement-operator>`_ for the formation of the radio-interferometric measurement operator;
2. `SARA-dictionary <https://github.com/basp-group/SARA-dictionary>`_ for the implementation of the sparsity priors.


.. warning::

   This project is still under active development.


.. toctree::
   :maxdepth: 2
   :caption: Installation & references

   setup
   contributors
   biblio

.. toctree::
   :maxdepth: 2
   :caption: API

   default
   _src/src
   _lib/lib
   _imaging/imaging
   _pyxisMs2mat/pyxisMs2mat


Indices and tables
==================

* :ref:`genindex`
* :ref:`mat-modindex`
* :ref:`search`
