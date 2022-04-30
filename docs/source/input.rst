Optional input datasets
=======================

.. note::

    The full path to the ``Faceted-HyperSARA`` repository is referred to as
    ``$FHS`` in the following lines.


Directory description
---------------------

The ``$FHS/data`` directory contains the ``.mat`` files generated during the data extraction step using ``pyxis`` packages.
It can also contain additional user input from a pre-processing step such as monochromatic joint (DDE/DIE) calibration and imaging step. 

The pre-processing input can include corrected noise estimates, initial image estimates, and DDE/DIE estimates. Considering a target source to be imaged tagged using ``$SRCNAME``, the pre-processing input is stored in the sub-directory ``$FHS/data/$SRCNAME/pre_processing/``.


General remarks
---------------

1. The user can opt for a different organisation of the auxilliary files, provided that this organisation is reflected in the filenames indicated in ``$FHS/imaging/main_input_imaging.m``.

2. All ``.mat`` files are associated with  a frequency channel and a dataset (if many), similarly to the data files 
generated during the data extraction step. 
 
3. If not available, the user should keep the filenames associated with the auxiliary input empty in ``$FHS/imaging/main_input_imaging.m``.


Input datasets
--------------

Corrected noise estimates
^^^^^^^^^^^^^^^^^^^^^^^^^

The different parameters involved in the definition of the optimisation task are derived from the noise statistics. More precisely, the :math:`\ell_2`-bounds of the data fidelity terms, and the regularisation parameters associated with the different model priors are linked to the standard deviation of the zero-mean noise, affecting the data. In practice, calibration errors can dominate the theoretical noise.  In such case, the user is advised, when possible, to provide "corrected" noise levels, which are to be applied on the **naturally** weighted visibilities. For instance, these estimates can be obtained as the standard deviation of the **naturally** weighted **residual** visibilities resulting from  a pre-processing joint calibration and imaging step.
When available, the noise estimates are stored ``.mat`` files.

.. code-block:: bash

    # Example: calibration: physical channel frequency id: 1, dataset nametag: $MSTAG 
    $FHS/data/$SRCNAME/pre_processing/sigma/$MSTAG/chs1_sigma.mat 

The ``.mat`` files contain the following field

.. code-block:: matlab

    "sigma" % (double) the corrected standard deviation of the noise.    
    "RESIDUAL" % (double complex vector) !! optional !! residual visibilities to be considered as a noise vector if the visibility gridding-based scheme is enabled.               

   
.. note::

   In the case of imperfect calibration, and non availability of reliable noise level estimates, the user is advised to activate the functionality for adaptive estimation of the noise level by setting the following parameter in  ``$FHS/imaging/main_input_imaging.m``.

   .. code-block:: bash

      param_global.adjust_flag_noise = true;
  
Image estimates
^^^^^^^^^^^^^^^

Estimated model images can be exploited to initialise the iterative algorithm, for acceleration purposes. These are saved in ``.fits`` files.  

.. code-block:: bash

    # Example: effective channel combining physical channels indexed from 1 to 16
    $FHS/data/$SRCNAME/pre_processing/model_images/chs1-16_model_image.fits 


Direction-independent effect (DIE) estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DIE estimates from a pre-processing step to be incorporated in the measurement operator. These are stored in ``.mat`` files.

.. code-block:: bash

    # Example: calibration: DIE, physical channel frequency id: 1, dataset nametag: $MSTAG 
    $FHS/data/$SRCNAME/pre_processing/cal_solutions/$MSTAG/chs1_cal_solutions.mat 

The ``.mat`` files contain the following field

.. code-block:: matlab

    "DIEs" % (double complex vector) multiplicative gains, of the same size as the data.                     


Direction-dependent effect (DDE) estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DDE estimates from a pre-processing can be exploited to obtain a more accurate measurement operator. These are 1D kernels in the spatial Fourier domain, provided for each visibility as the result of the convolution of the calibration kernels of the associated antenna pair. Since DDE estimates are provided in the spatial Fourier domain, the  DDE calibration pre-processing step must have been performed over the exact same field of view of interest. On a further note, the Fourier shifting convention must be exactly the same as the one adopted for the NUFFT interpolation kernel. Finally, DDE solutions have the same support size.  

DDE estimates are stored in ``.mat`` files.

.. code-block:: bash

    # Example: calibration: DDE, physical channel frequency id: 1, dataset nametag: $MSTAG 
    $FHS/data/$SRCNAME/pre_processing/cal_solutions/$MSTAG/chs1_cal_solutions.mat 

The ``.mat`` files contain the following field

.. code-block:: matlab

    "DDEs" % (complex array), each row corresponds to a DDE kernel, associated with a visibility.                   
