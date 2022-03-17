Faceted-Hyper-SARA
==================

``Faceted-Hyper-SARA`` is a MATLAB wideband imaging library for radio-interferometry. The library contains an implementation of
algorithms from the "SARA" family (SARA :cite:p:`Carrillo2012,Onose2017`,
HyperSARA :cite:p:`Abdulaziz2019` and Faceted HyperSARA 
:cite:p:`Thouvenin2021`). The proposed implementation can accommodate

- visibility gridding :cite:p:`Kartik2017`;
- :math:`w`-correction :cite:p:`Dabbech2018`;
- a compact Fourier model for the direction dependent effects (DDEs) :cite:p:`Dabbech2021`.

``Faceted-Hyper-SARA`` relies on two auxiliary submodules:

1. `measurement-operator <https://github.com/basp-group-private/measurement-operator>`_ for the implementation of the radio-interferometry measurement operator;
2. `faceted-wavelet-transform <https://github.com/basp-group-private/faceted-wavelet-transform>`_ for the implementation of the priors.


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

.. note::

   - This project is under active development.



Indices and tables
==================

* :ref:`genindex`
* :ref:`mat-modindex`
* :ref:`search`
