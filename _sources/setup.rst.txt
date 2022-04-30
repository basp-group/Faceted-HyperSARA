Setup
=====

Installation
------------

Environment setup
^^^^^^^^^^^^^^^^^

- If you use ``https``, first edit the ``.gitmodules`` file as follows

.. code-block:: rst

   [submodule "lib/SARA-dictionary"]
      path = lib/SARA-dictionary
      url = https://github.com/basp-group/SARA-dictionary.git
   [submodule "lib/RI-measurement-operator"]
      path = lib/RI-measurement-operator
      url = https://github.com/basp-group/RI-measurement-operator.git


You can then clone the repository with all the submodules as follows

.. code-block:: bash

   # Cloning the repo. with all the submodules:
   git clone --recurse-submodules  https://github.com/basp-group/Faceted-HyperSARA.git
   cd Faceted-HyperSARA


- If you are using a properly configured ``ssh`` key for Github, you can clone the repository with all the submodules as follows

.. code-block:: bash

   # Cloning the repo. with all the submodules:
   git clone --recurse-submodules git@github.com:basp-group/Faceted-HyperSARA.git
   cd Faceted-HyperSARA


.. important::
   
   The full path to the ``Faceted-HyperSARA`` repository is referred to as ``$FHS`` in the rest of the documentation.


Updating submodules (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From ``$FHS``, issue the following commands

.. code-block:: bash

   cd $FHS
   git pull
   git submodule sync --recursive # update submodule address, in case the url has changed
   git submodule update --init --recursive # update the content of the submodules
   git submodule update --remote --merge # fetch and merge latest state of the submodule

.. warning::

   A simple ``pull`` instruction from the ``Faceted-HyperSARA`` repository does not update the submodules if the version on which it relies has been modified. If this is the case, update the submodules with the instructions provided above.



Getting started
---------------


Reconstructing an image cube from an MS-Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Start by extracting the data from an ``MS Table`` using the Python script 
   ``$FHS/pyxisMs2mat/pyxis4DataExtraction.py``, that is to be launched from  the sub-directory ``$FHS/pyxisMs2mat/``. The data will be extracted as a
   collection of ``.mat``, saved in the sub-directory ``$FHS/data``. Instructions and examples are provided in the :doc:`_pyxisMs2mat/pyxisMs2mat` page (also available in the standalone ``$FHS/pyxisMs2mat/ReadMe.md`` file). Data extraction requires the `casacore <https://github.com/casacore/casacore>`_ and `meqtrees <https://github.com/ratt-ru/meqtrees/wiki/Installation>`_ libraries.

.. code-block:: matlab

   %%  the data files contain the following fields. 
   % Note that flagged data are already discarded.
   "frequency"  % channel frequency                       
   "y"  % data (Stokes I)
   "u"  % u coordinate (in units of the wavelength)
   "v"  % v coordinate (in units of the wavelength)
   "w"  % w coordinate (in units of the wavelength)                       
   "nW"  % sqrt(weights)
   "nWimag" % (optional field) imaging weights if available (Briggs or uniform), empty otherwise
   "maxProjBaseline"  % max projected baseline (in units of the wavelength)


2. (Optional) Make a copy of ``$FHS/imaging/default_parameters.json`` and 
   update the main algorithm parameters specified in this file (all values are 
   set to default, reliable values). Documentation for all the parameters involved is given in the :doc:`default` page.

3. (Optional) In order to use pre-estimated calibration kernels, wideband cube
   and noise level estimates, read carefully the instructions reported in the
   :doc:`input` page.

4. Configure :mat:scpt:`imaging.main_input_imaging` following the instructions
   provided in the documentation and the file itself. Blocks of the variables to be configured are preceded with the instruction ``%/ TODO: to be adjusted by the user``.

5. Run ``main_input_imaging.m`` from ``$FHS/imaging/`` sub-directory with MATLAB. Results will be saved in the sub-directory ``$FHS/results/``.


Testing the library with a synthetic dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test the imaging pipeline with a synthetic data test similar to the one considered in :cite:p:`Thouvenin2021`, follow the instructions below.

1. Retrieve the 
   `S_DDE_MODEL.fits <https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits>`_ image file associated with :cite:p:`Dabbech2021`.
   Place the file in the ``$FHS/data/`` folder at the root of the library.

   .. code-block:: bash
 
      # if on MAC:
      # brew install wget
      cd $FHS/
      mkdir data && cd data
      wget -P . https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits

2. Generate one (or all) synthetic wideband image cube used in
   :cite:p:`Thouvenin2021` using the
   :mat:scpt:`imaging.main_generate_cyga_cubes` script. Data cubes and auxiliary matlab files will be saved in ``$FHS/data/``.

   .. code-block:: matlab

      main_generate_cyga_cubes

3. Generate a synthetic dataset using 
   :mat:func:`imaging.main_generate_data`. The two datasets considered in :cite:p:`Thouvenin2021` can be generated by running the following MATLAB instructions

   .. code-block:: matlab

      % generate data for the spatial faceting experiment
      main_generate_data('default_parameters.json', 'cygA', 8, ...
      '../data/msSpecs.mat', 'spatial', 2, 40, false, ...
      "local", false)

      % generate data for the spectral faceting experiment
      main_generate_data('default_parameters.json', 'cygA', 8, ...
      '../data/msSpecs.mat', 'spectral', 2, 40, false, ...
      "local", false)

4. (Optional) Make a copy of ``$FHS/imaging/default_parameters.json``, and 
   update the main algorithm parameters specified in this file (all values are 
   set to default, reliable values). Documentation for all the parameters involved is given in the :doc:`default` page.

5. Configure :mat:scpt:`imaging.main_input_imaging` following the instructions
   provided in the documentation and the file itself. Blocks of variables to be configured are preceded with the instruction ``% TODO: to be adjusted by the user``. Example configuration used for the experiments reported in :cite:p:`Thouvenin2021` is provided in the ``$FHS/imaging/main_input_imaging_synth.m`` script.

6. Run ``main_input_imaging.m`` from ``$FHS/imaging/`` sub-directory with MATLAB. Results will be saved in the sub-directory ``$FHS/results/``.


Contributing
------------

- Issue Tracker: `https://github.com/basp-group/Faceted-HyperSARA/issues <https://github.com/basp-group/Faceted-HyperSARA/issues>`_
- Source Code: `https://github.com/basp-group/Faceted-HyperSARA <https://github.com/basp-group/Faceted-HyperSARA>`_

To contribute to the project, make sure the following elements are properly
configured before submitting any pull request (PR).


Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Make sure any new functionality is properly documented using the ``numpy``
  docstring style.
- To build the documentation, issue the folowing commands.

.. code-block:: bash

   # setup conda environment to build the documentation
   conda env create --name fhs-doc --file environment.yml 

   # alternative using conda/pip
   # conda create -n fhs-doc
   # conda activate fhs-doc
   # conda install pip
   # pip install miss_hit
   # pip install -r requirement.txt

   # building the documentation in html format
   cd docs
   make html

- All the generated ``.html`` files are contained in the ``$FHS/docs/build`` folder.
- If needed, you can delete the ``conda`` environment as follows

.. code-block:: bash
   
   conda env remove -n fhs-doc


Pushing the documentation online
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add a ``worktree`` from the ``master`` branch

.. code-block:: bash

   # make sure the folder html does not exist before running the command
   git worktree add $FHS/docs/build/html gh-pages
   cd $FHS/docs/build/html
   git add .
   git commit -m "Build documentation as of $(git log '--format=format:%H' master -1)"
   git push origin gh-pages
   # delete the worktree
   cd ../
   git worktree remove html


Code layout
^^^^^^^^^^^

If you contribute code to the library (through a `pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests>`_), make sure any submitted code is properly formatted with the `miss_hit <https://pypi.org/project/miss-hit/>`_ package using the provided ``miss_hit.cfg`` configuration file

.. code-block:: bash

   # activate fhs-doc environment (see previous paragraph)
   conda activate fhs-doc
   # run the following command from the root of the package (where the miss_hit.cfg file is)
   mh_style --fix .
