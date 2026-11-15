"""Level 2 - Lambert solver on the GPU.

Plan: compile ``itokawa_opt.astro.lambert.lambert_izzo`` (plain arithmetic, bounded
iterations) as a ``numba.cuda.jit`` device function, one thread per Lambert problem;
then a CuPy RawKernel version in CUDA C++. Compare FP64 vs FP32, measure warp divergence.
"""
