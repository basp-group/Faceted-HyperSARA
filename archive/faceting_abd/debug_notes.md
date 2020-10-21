# [21/10/2020] List of functions / folders used in Abdullah's older codes (mnras branch, test7_20_triangular_d256)

**Caveat**: apparently Abdullah re-generated the image cube with the older functions (which may explain why he doesn't see the same issue as I do)

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
   1. ...
   2. ...
   3. ...

- `addpath ../hypersara-clean/lib/generate_data/`
- `addpath ../hypersara-clean/lib/operators/`
- `addpath ../hypersara-clean/lib/CubeHelix/`
- `addpath(irt_library)`
- `addpath ../sdwt2/`
- `addpath ../src/`
- `addpath ../data/`
- `addpath ../hypersara-clean/lib/Proximity_operators/code/matlab/indicator/`
- `addpath ../hypersara-clean/lib/Proximity_operators/code/matlab/multi/`
