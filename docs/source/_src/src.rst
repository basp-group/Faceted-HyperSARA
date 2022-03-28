src package
===========

Core implementation of the SARA, HyperSARA and Faceted HyperSARA imaging algorithms.

.. contents:: :local:

.. toctree::
   :maxdepth: 2
   :caption: List of modules
   
   src.fhs
   src.heuristics
   src.hs
   src.sara

API src
-------

.. mat:automodule:: src
    :members: apply_adjoint_operator, apply_direct_operator, apply_scaled_fourier_transform, compute_residual_images, compute_sara_prior, nuclear_norm, proj_l2ball, sanity_check, update_dual_data_fidelity
