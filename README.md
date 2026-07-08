# Randomized Sketching in Quadrature-Based Restarting for Matrix Functions

## About

This MATLAB code implements a sketch-and-restart framework for computing  $f(A)b$, where $A$ is a large sparse matrix, $f$ is a suitable matrix function, and $b$ is a vector.

This implementation is based on the original [funm_quad](https://github.com/guettel/funm_quad) code for quadrature-based restarts for matrix functions.

## Reproducing Section 6

The numerical tests in Section 6 can be reproduced by running the following
simple driver files from the repository root:

```matlab
test_conv_diff
test_qcd
test_frac_laplacian
```

Each `test_*` driver loads the corresponding settings file, runs all methods, prints the LaTeX table to the MATLAB command window, and regenerates the paper figures. The randomized sketching is seeded in `run_single_method.m` with
`rng(1)`, so the random choices are reproducible across runs on the same MATLAB
version and platform.

If the `.mat` result files already exist and only the table or figures need to
be regenerated, use the lighter helper scripts:

```matlab
print_conv_diff
plot_conv_diff

print_qcd
plot_qcd

print_frac_laplacian
plot_frac_laplacian
```

## Main Code Organization

- `funm_quad/`: core restarted Krylov and quadrature routines.
- `run_single_method.m`: dispatches one method and fixes the random seed.
- `run_methods.m`: runs all non-thick-restart methods for one example.
- `run_methods_thick_restart.m`: runs the thick-restart variants.
- `plot_figures.m` and `plot_figures_thick_restart.m`: generate `.eps` and
  `.png` figures.
- `print_table.m`: prints the LaTeX table rows used in Section 6.
