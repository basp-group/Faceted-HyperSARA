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


Experiments
-----------

To run a representative example of the experiments, configure and run the `main.py` script.

.. code-block:: bash

   conda activate async_sampling
   cd path/to/repo
   # to check the list of input parameters
   # ...


Contributing
------------

- Issue Tracker: `https://github.com/basp-group-private/Faceted-Hyper-SARA/issues <https://github.com/basp-group-private/Faceted-Hyper-SARA/issues>`_
- Source Code: `https://github.com/basp-group-private/Faceted-Hyper-SARA <https://github.com/basp-group-private/Faceted-Hyper-SARA>`_

To contribute to the project, make sure the following elements are properly configured before submitting any pull request (PR).


Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Make sure any new functionality is properly documented using the ``numpy`` docstring style.
- As soon as ``sphinx`` is installed (``conda install -c anaconda sphinx``), issue the following commands.

.. code-block:: bash

   conda install -c anaconda sphinx
   conda install -c conda-forge sphinxcontrib-bibtex
   pip install -U sphinxcontrib-matlabdomain
   cd build docs/build/html
   make html
