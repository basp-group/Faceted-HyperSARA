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

From an existing ``Faceted-Hyper-SARA`` repository, issue the following commands

.. code-block:: bash

   cd path/to/repo
   git pull
   git submodule sync --recursive # update submodule address, in case the url has changed
   git submodule update --init --recursive # update the content of the submodules
   git submodule update --remote --merge # fetch and merge latest state of the submodule

.. warning::

   Note that a simple ``pull`` instruction from the ``Faceted-Hyper-SARA`` repository does not update the submodules if the version on which it relies has been modified. If this is the case, update the submodules with the ``git submodule``instructions provided above.


Getting started
---------------


Reconstructing an image cube from an MS-Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. To reconstruct an image from an ``MS Table``, first use the Python script 
   ``pyxis_ms2mat/pyxis-ckat.py`` to extract the data in the form of a 
   collection of ``.mat`` (an example is provided in
   ``pyxis_ms2mat/job_example_slurm``). This requires the follownig libraries: 
   ``casacore`` and ``meqtrees``.

2. (Optional) Make a copy of ``imaging/default_parameters.json``, and 
   update the main algorithm parameters specified in this file (all values are 
   set to default, reliable values). Documentation for all the parameters involved is given in the :doc:`default` page.

3. Configure :mat:scpt:`imaging.main_input_imaging` following the instructions
   provided in the documentation and the file itself. Blocks of variables to be configured are preceded with the instruction ``% TODO: to be adjusted by the user``.

4. Run ``./imaging/main_input_imaging.m``` with MATLAB.


Testing the library with a synthetic dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test the imaging pipeline with a synthetic data test.

1. Retrieve the 
   `S_DDE_MODEL.fits <https://researchportal.hw.ac.uk/files/43645966/S_DDE_MODEL.fits>`_ image file associated with :cite:p:`Dabbech2021`.
   Place the file in the ``data/`` folder at the root of the library.

2. Generate one (or all) synthetic wideband image cube used in
   :cite:p:`Thouvenin2021` using the
   :mat:scpt:`imaging.main_generate_cyga_cubes` script, following the
   instructions.

3. Generate a synthetic dataset using 
   :mat:func:`imaging.main_generate_data`.

4. (Optional) Make a copy of ``imaging/default_parameters.json``, and 
   update the main algorithm parameters specified in this file (all values are 
   set to default, reliable values). Documentation for all the parameters involved is given in the :doc:`default` page.

5. Configure :mat:scpt:`imaging.main_input_imaging` following the instructions
   provided in the documentation and the file itself. Blocks of variables to be configured are preceded with the instruction ``% TODO: to be adjusted by the user``.

6. Run :mat:scpt:`imaging.main_input_imaging` with MATLAB.


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

- All the generated ``.html`` files are contained in the ``docs/build`` folder.
- If needed, you can delete the ``conda`` environment as follows

.. code-block:: bash
   
   conda env remove -n fhs-doc


Pushing the documentation online
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add a ``worktree`` from the ``master`` branch

.. code-block:: bash

   # make sure the folder html does not exist before running the command
   git worktree add docs/build/html gh-pages
   cd docs/build/html
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

   # activate sdwt-doc environment (see previous paragraph)
   conda activate sdwt-doc
   # run the following command from the root of the package (where the miss_hit.cfg file is)
   mh_style --fix .
