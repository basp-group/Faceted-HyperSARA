Default parameters
==================

A number of parameters can be accessed by the user to apply one of the three implemented imaging approaches.

.. https://stackoverflow.com/questions/40748886/how-can-i-document-a-constant-module-level-variable-with-sphinx-docstring-wit
.. https://stackoverflow.com/questions/9162891/define-mark-up-for-generic-sphinx-admonitions-with-a-specific-title
.. https://docutils.sourceforge.io/docs/ref/rst/directives.html#generic-admonition

User-defined parameters
-----------------------

All parameters which explicitly require to be set by the user are defined in the :mat:scpt:`imaging.main_input_imaging` script. Additional instructions on how to set these values are provided in the documentation below and in the source code of this script, to be modified by the user.

.. mat:autoscript:: imaging.main_input_imaging

Default values of the parameters listed below are included in ``imaging/default_parameters.json``. These will be utilised unless the user updates them directly or alternatively specifies them in :mat:scpt:`imaging.main_input_imaging`, in such case, they will be overwritten.

.. admonition:: Global parameters

    reg_flag_reweighting: true
        Flag to activate the re-weighting procedure. 
    measop_flag_wproj: false
        Flag to activate `w`-projection.
    measop_flag_visibility_gridding: false
        Flag to activate visibility gridding for data dimensionality reduction.
    adjust_flag_noise: false
        Flag to estimate the effective noise level on the fly, activate only when reliable noise estimates are not available.
    data_flag_apply_imaging_weights: false
        Flag to apply imaging weights to the data (such as uniform or Briggs, if available) and to incorporate them in the measurement operator. Natural weighting is adopted by default.
    algo_flag_computeOperatorNorm: true 
        Flag to compute the spectral norm of the measurement operator.
    algo_flag_saveOperatorNorm: true
        Flag to save the spectral norm of the measurement operator.
    algo_flag_solveMinimization: true
        Flag to run imaging.
    parcluster: "local"
        Name of the parcluster profile.
    preproc_filename_noise_std: []
        Filename of the estimated noise level.
    preproc_filename_cal_solutions: []
        Filename of the DDE calibration kernels. 
    preproc_filename_model: []
        Filename of the initial image estimate.
            
  


Inner algorithmic parameters
----------------------------

The parameters listed below, defined in the ``imaging/default_parameters.json`` file, specify the configuration in which the reweighting scheme and the primal-dual forward-backward (PDFB) algorithm are applied in the selected imaging approach. They also include the parameters involved in the `w`-projection to update the measurement operator, and the parameters involved in the estimation of the noise level on the fly in the PDFB algorithm, when available noise statistics are unreliable. 


Default values can be overwritten by providing the name of a new ``.json`` file to the :mat:scpt:`imaging.main_input_imaging` configuration script (see the variable ``json_filename``). Note that all the fields specified in ``imaging/default_parameters.json`` are required to be able to run the imaging pipeline.

.. warning::
    The user is encouraged to use the default values provided in ``imaging/default_parameters.json`` (particularly the list provided below) unless precisely aware of their impact on the reconstruction process.



.. admonition:: Reweighting scheme :cite:p:`Candes2009`

    min_iter : 1
        Minimum number of reweighting iterations.
    max_iter : 10
        Maximum number of reweighting iterations.
    rel_var : 1e-4
        Stopping criterion (relative variation of the objective function over two consecutive reweighting steps).
    alpha : 1
        Additional 
    backup_frequency : 1
        Number of iterations after which all the variables within the algorithm are saved (for a warm-restart).
    verbose : 2
        Run the algorithm in verbose mode (0: no information displayed, 1: partial information displayed, 2: most verbose state).

    See :mat:func:`src.sara.sara`, :mat:func:`src.hs.hyperSARA` and :mat:func:`src.fhs.facetHyperSARA` for further details.


.. admonition:: Primal-dual forward-backward algorithm (PDFB) :cite:p:`Condat2013,Vu2013,Pesquet2014`
    
    min_iter : 10
        Minimum number of iterations.
    max_iter : 2000
        Maximum number of iterations.
    rel_var : 1e-5
        Stopping criterion (relative variation of the objective function). This criterion needs to be satisfied jointly with the data-fidelity criterion encoded in ``fidelity_tolerance``.
    fidelity_tolerance : 1.01
        Data fidelity cirterion (norm of the residual in each data block is within ``fidelity_tolerance`` of the radius of the associated :math:`\ell_2`-ball constraint). This criterion needs to be satisfied jointly with the relative variation criterion encoded in ``rel_var``.
    rel_var_low : 5e-6
        Additional relative variation criterion, taking precedence over ``rel_var`` and ``fidelity_tolerance``. This criterion avoids configurations in which the iterates are trapped in a local minimum, and very slowly evolve to satisfy the data-fidelity constraints.

    See :mat:func:`src.sara.sara`, :mat:func:`src.hs.hyperSARA` and :mat:func:`src.fhs.facetHyperSARA` for further details.


.. admonition:: Non-uniform fast Fourier transform (NUFFT) :cite:p:`Fessler2003`
    
    ox : 2
        Fourier oversampling factor along the axis x.
    oy : 2
        Fourier oversampling factor along the axis y.
    Kx : 7
        Size of the NUFFT interpolation kernel along axis x.
    Ky : 7
        Size of the NUFFT interpolation kernel along axis y.
    kernel : "minmax:tuned"
        Name of the selected NUFFT interpolation kernel. Possible options include ``kaiser``, ``minmax:kb`` and ``minmax:tuned``. See associated documentation of the `lib.operators.op_nufft function <https://basp-group.github.io/RI-measurement-operator/_lib/lib.operators.html#lib.operators.op_nufft>`_ from the RI-measurement-operator module.


.. admonition:: Preconditioning (instrumental in PDFB) :cite:p:`Onose2017`
    
    gen_uniform_weight_matrix : true
        Flag to activate the generation of uniform weights (to be kept active)
    uniform_weight_sub_pixels : 1
        Parameter to consider sub-pixel weights in the uniform weighting scheme.

    See associated documentation of the `lib.utils.util_gen_preconditioning_matrix function <https://basp-group.github.io/RI-measurement-operator/_lib/lib.utils.html#lib.utils.util_gen_preconditioning_matrix>`_ from the RI-measurement-operator module.


.. admonition:: Ellipsoid projection (instrumental in PDFB) :cite:p:`Onose2017`

    min_iter : 1
        Minimum number of iterations.
    max_iter : 20
        Maximum number of iterations.
    eps : 1e-8
        Stopping criterion based on the relative variation of the objective function (associated with the projection problem).

    See associated documentation of the `lib.utils.solver_proj_elipse_fb <https://basp-group.github.io/RI-measurement-operator/_lib/lib.utils.html#lib.utils.solver_proj_elipse_fb>`_ from the RI-measurement-operator module.


.. admonition:: Wavelet dictionary (SARA dictionary by default) :cite:p:`Carrillo2012`

    basis : ["db1", "db2", "db3", "db4", "db5", "db6", "db7", "db8", "self"]
        Name of the wavelet dictionaries considered ("self" corresponding to the Dirac basis). By default, contains the list of wavelets defining the SARA dictionary :cite:p:`Carrillo2012`. Whenever used, the Dirac basis needs to be specified in last.
    nlevel : 4
        Number of decomposition scales considered.
    filter_length : [2, 4, 6, 8, 10, 12, 14, 16, 0]
        Length of the filters corresponding to the selected wavelet dictionaries, where by convention 0 corresponds to the Dirac basis.

    See associated documentation in the `SARA-dictionary <https://basp-group.github.io/SARA-dictionary/index.html>`_ module.

.. admonition:: Noise estimation on the fly

    min_iter : 100
        Minimum number of iterations.
    rel_var : 1e-3
        Tolerance to adjust the noise estimate.
    start_iter: 500
        Iteration numbers to force triggering the noise estimation.
    change_percentage: 0.5
        The weight of the update w.r.t the l2 norm of the residual data.
    start_change_percentage: 0.1
        The weight of the update w.r.t the l2 norm of the residual data, if noise level update triggered at `start_iter` .
    




.. admonition:: :math:`w`-projection parameters :cite:p:`Dabbech2017`

    measop_wprojCEnergyL2 : 0.9999
        Sparsification levels for the :math:`w` kernel, to be selected in the interval :math:`[0.99, 1]`. 
    measop_wprojGEnergyL2 : 0.9999
        Sparsification levels for convolution kernels involved in the degridding matrix :math:`G` after :math:`w`-correction. To be selected in the interval :math:`[0.99, 1]`.



