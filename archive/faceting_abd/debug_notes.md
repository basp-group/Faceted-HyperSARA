# [21/10/2020] List of functions / folders used in Abdullah's older codes (mnras branch, test7_20_triangular_d256)

**Caveat**: apparently Abdullah re-generated the image cube and uv-coverage with the older functions (which may explain why he doesn't see the same issue as I currently do).

## List of functions used (per script) (to be checked and compared to the current ones)

  1. Generate_cube
  2. Generate_undersampled_cube_new
  3. interleaved_facets
  4. domain_decomposition
  5. Generate_data_new ->
      - generate_uv_coverage2
      - util_gen_preconditioning_matrix
      - util_gen_block_structure
      - op_p_nufft
  6. Generate_Measurements ->
      - util_gen_measurements
      - util_gen_data_fidelity_bounds
  7. Solver_simulated_data_FINAL_clean ->
      - HS_forward_operator_precond_G
      - HS_adjoint_operator_precond_G
      - fb_nnls_blocks
      - domain_decomposition
      - pdfb_LRJS_precond_NL21_sdwt2_spmd4_weighted

## List of functions used (per folder)

- `addpath ../hypersara-clean/lib/`
- `addpath ../hypersara-clean/lib/generate_data/`
- `addpath ../hypersara-clean/lib/operators/`
- `addpath(irt_library)`
- `addpath ../sdwt2/`
- `addpath ../src/`
- `addpath ../data/`
- `addpath ../hypersara-clean/lib/Proximity_operators/code/matlab/indicator/`
- `addpath ../hypersara-clean/lib/Proximity_operators/code/matlab/multi/`

## Comparison between the different files and folder (manual)

- `sdwt2/`: domain_decomposition to be checked (all the rest is identical, up to comments)
- `nufft/`: all files are identical, but my version has all the files required (including `.mat`)
- `operators/`: HS_forwad_operator_precond_G, op_p_nufft, so_fft2_adj
- `generate_data/`: Generate_cube, Generate_cube_new, generate_uv_coverage -> ok, no significant difference (most difference in the way the uv-coverage is obtained)
- `src/`: only functions to be compared for the moment are thos used in the weigthed version of the algo
  1. domain_decomposition: ok (sdwt2)
  2. generate_segdwt_indices: ok (sdwt2)
  3. domain_decomposition_overlap2: ok (sdwt2)
  4. initialize_dual_variables_prior_cst_overlap: ok
  5. update_primal: ok
  6. comm2d_update_ghost_cells: ok (sdwt2)
  7. update_nuclear_spmd_weighted -> update_dual_nuclear: ok
  8. update_l21_spmd -> update_dual_l21: ok
  9. comm2d_reduce: ok (sdwt2)
  10. update_data_fidelity: ok -> update_dual_fidelity
  11. prior_overlap_spmd_cst3_weighte -> spmd/weighted -> d: ok
  12. update_weights_overlap2_weighted: ok
  13. compute_residual_images: ok
  14. pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted: same inner functions used (those mentioned above), only differences in warm restart, use of some auxiliary variables

## Remaining differences (out of the solver itself)

[21/10/2020]

- uv-coverage? (Abdullah currently modifying, using meqtrees instead of the current MATLAB codes)
- image cube: possible issue with the spectra? (Abdullah currently running a new test to check this point)
- type of prior? (need to add "facets" on the border of the field of view) -> would need to rerun the real data expriments

[28/10/2020]

- results from Abdullah are really good: make sure l2 constraint is "sufficiently" satisfied before starting the first reweight (increase max number of PPD iterations + consider l2 constraint in the stopping criterion for PPD)
  - compute duality gap to properly check convergence for each PPD algo?

- function within the solver to be debugged: (co_w version, /weighted/) [done]
  - spmd/update_primal
  - spmd/update_nuclear_spmd_weighted -> update_dual_nuclear
  - spmd/update_l21_spmd -> update_dual_l21
  - spmd/prior_overlap_spmd_cst3_weighted -> spmd/weighted ->  -> compute_facet_prior_overlap
  - spmd/update_weights_overlap2_weighted -> spmd/weighted -> update_weights_overlap
  - spmd/initialize_dual_variables_prior_cst_overlap -> initialize_dual_overlap
  - update_data_fidelity -> update_dual_fidelity
  - update_epsilon
  - compute_residual_images

- co version (/weighted/) [done]
  - spmd/initialize_dual_variables_prior_cst_overlap -> initialize_dual_overlap
  - spmd/update_primal
  - spmd/update_nuclear_spmd -> update_dual_nuclear
  - spmd/update_l21_spmd -> update_dual_l21
  - spmd/weighted/prior_overlap_spmd_cst -> compute_facet_prior_overlap
  - spmd/weighted/update_weights_cst_overlap -> update_weights_overlap
  - update_epsilon
  - update_data_fidelity -> update_dual_fidelity
  - compute_residual_images

- w version (/weighted/) [trashed] [moved to trash, less general than co_weighted]
  - spmd/weighted/initialize_dual_variables_prior_overlap -> initialize_dual_overlap
  - spmd/weighted/prior_overlap_spmd_weighted -> compute_facet_prior_overlap
  - spmd/weighted/update_weights_overlap_weigthed
  - spmd/update_primal
  - spmd/update_nuclear_spmd_weighted -> update_dual_nuclear -> create a single function, with auxiliary parameter
  - spmd/update_l21_spmd -> update_dual_l21
  - update_data_fidelity -> update_dual_fidelity
  - update_epsilon
  - compute_residual_images

---

- no overlap [done]
  - spmd/no/initialize_dual_variables_prior -> initialize_dual_variables
  - spmd/no/prior_value_spmd -> compute_facet_prior
  - spmd/no/update_weights
  - spmd/update_primal
  - spmd/update_nuclear_spmd -> update_dual_nuclear
  - spmd/update_l21_spmd -> update_dual_l21
  - update_data_fidelity -> update_dual_fidelity
  - update_epsilon
  - compute_residual_images

---

- standard (same overlap for nuclear and l21) [done]
  - spmd/standard/prior_overlap_spmd -> compute_facet_prior_so
  - spmd/standard/update_weights_overlap -> update_weights_so
  - spmd/initialize_dual_variables_prior_overlap -> initialize_dual_overlap
  - spmd/update_primal
  - spmd/update_nuclear_spmd -> update_dual_nuclear
  - spmd/update_l21_spmd -> update_dual_l21
  - update_data_fidelity -> update_dual_fidelity
  - update_epsilon
  - compute_residual_images

---

- serial (can be kept as is) [done]
  - serial/initialize_l21_serial
  - serial/run_par_nuclear -> update_dual_nuclear_serial
  - serial/run_par_l21 -> update_dual_l21_serial
  - serial/update_weights_nuclear_serial
  - serial/update_weights_l21_serial
  - serial/nuclear_norm
  - serial/l21_norm_sara -> compute_sara_prior
  - update_data_fidelity -> update_dual_fidelity
  - update_epsilon
  - compute_residual_images

---

- from the sdwt2 library
  - domain_decomposition -> split_range
  - generate_segdwt_indices
  - domain_decomposition_overlap2 -> split_range
  - comm2d_update_ghost_cells -> comm2d_update_borders
  - comm2d_reduce

---

[29/10/2020] Start debugging src_new