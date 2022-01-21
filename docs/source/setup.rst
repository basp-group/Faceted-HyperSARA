Setup
=====

Installation
------------

Environment setup
^^^^^^^^^^^^^^^^^

Clone the current repository with all the submodules as follows

.. code-block:: bash

   # Cloning the repo. with all the submodules:
   git clone --recurse-submodules https://github.com/basp-group-private/Faceted-Hyper-SARA.git
   cd Faceted-Hyper-SARA


Updating submodules (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From an existing `faceted-hyper-sara` repository, issue the following commands

.. code-block:: bash

    cd path/to/repo
    git pull
    git submodule sync --recursive # update submodule address, in case the url has changed
    git submodule update --init --recursive # update the content of the submodules
    git submodule update --remote --merge # fetch and merge latest state of the submodule


Getting started
---------------


Reconstructing an image cube from an MS-Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. To reconstruct an image from an ``MS Table``, first use the Python script 
   ``pyxis_ms2mat/pyxis-ckat.py`` to extract the data in the form of a 
   collection of ``.mat`` (an example is provided in
   ``pyxis_ms2mat/job_example_slurm``). This requires the follownig libraries: 
   ``owlcat``, ``pyxis``, ``py-casacore``, ``tigger``, ``meqtrees-cattery``,
   ``meqtrees-timba``.

2. (Optional) Make a copy of ``experiments/default_parameters.json``, and 
   update the main algorithm parameters specified in this file (all values are 
   set to default, reliable values).

3. Configure ``experiments/main_input_exp.m`` following the instructions 
   provided in the file. Additional documentation is available in the 
   documentation of the :mod:`experiments` module.

4. Run the ``experiments/main_input_exp.m`` with MATLAB.

.. .. code-block:: bash

..    conda activate async_sampling
..    cd path/to/repo
..    # to check the list of input parameters
..    # ...

Testing the library with a synthetic dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test the imaging pipeline with a synthetic data test.

1. Generate a synthetic dataset using 
   :func:`experiments/main_generate_data`, and run the (properly
   configured) `experiments/main_input_exp.m` script (see 
   :ref:`previous section<Reconstructing an image cube from an MS-Table>`).

2. Configure ``experiments/main_input_exp.m`` following the instructions 
   provided in the file. Additional documentation is available in the 
   documentation of the :mod:`experiments` module.


Contributing
------------

- Issue Tracker: `https://github.com/basp-group-private/Faceted-Hyper-SARA/issues <https://github.com/basp-group-private/Faceted-Hyper-SARA/issues>`_
- Source Code: `https://github.com/basp-group-private/Faceted-Hyper-SARA <https://github.com/basp-group-private/Faceted-Hyper-SARA>`_

To contribute to the project, make sure the following elements are properly
configured before submitting any pull request (PR).


Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Make sure any new functionality is properly documented using the ``numpy``
  docstring style.
- To build the documentation, issue the folowing commands.

.. code-block:: bash

   # setup conda environment to build the documentation
   conda create -n fhs-doc
   conda activate fhs-doc
   pip install sphinx sphinx_rtd_theme sphinxcontrib-bibtex sphinxcontrib-matlabdomain
   # building the documentation in html format
   cd docs
   make html

- All the generated ``.html`` files are contained in the ``docs/build`` folder.
- If needed, you can delete the ``conda`` environment as follows

.. code-block:: bash
   
   conda env remove -n fhs-doc
