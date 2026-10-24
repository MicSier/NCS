# CPU emulation build

Runs NCS and its examples on a machine without an NVIDIA GPU (or without nvcc),
executing the CUDA kernels in `NCS.h` unchanged. Useful for development on a
laptop, CI, and regression tests of the kernels themselves.

```sh
make cpu          # generates build-cpu/ and builds *_cpu binaries with g++ + OpenMP
make cpu-test     # kernel / integrator tests (reductions, barrier semantics, quadrature)
./singlet_cpu     # same output as ./singlet on a GPU
```

## How it works

- `cpu/cuda_runtime.h` stands in for the CUDA runtime: `__global__`/`__device__`/
  `__constant__` become plain C++, `cudaMalloc`/`cudaMemcpy` become `malloc`/`memcpy`,
  and the device reports 1024 threads per block.
- `cpu/cudaify.py` copies `NCS.h`, `gl.h`, `examples/` and `cpu/tests/` into
  `build-cpu/`, rewriting `kernel<<<grid, block[, shmem]>>>(args)` into
  `cpu_launch(cpu_cfg(grid, block[, shmem]), kernel, args)` and
  `extern __shared__ T x[];` into a per-block buffer. Sources in the repo are not touched.
- Each block runs on one OS thread (blocks are spread over OpenMP threads). Each CUDA
  thread runs on a small fiber, and `__syncthreads()` is a real barrier. If thread 0
  finishes without reaching a barrier, the rest of the block runs as a plain loop.
- Linux x86-64 uses a hand-written context switch; other platforms use `<ucontext.h>`.

Speed: ~7 ms per `sc_integrate1D_gl_gpu` call on 2 cores (`IntegrationBenchmark_cpu`).

Limitations: no warp-level primitives (`__shfl*`, `__ballot*`), atomics, streams or
events; a `__syncthreads()` that only some threads of a block reach is not detected.

## Regression against the thesis tables

`regress_rashba_cpu` recomputes rows of `data/rashba/tc/coupled+{s,t}+...dat`
(Vs=0.5, Vt=2, nl=1.5) with the current coupled solver:

```sh
./regress_rashba_cpu "data/rashba/tc/coupled+t+V_s=0.5+V_t=2+nl=1.5+t'=0.dat" t 5
NCS_BRACKET=1.005 ./regress_rashba_cpu "data/rashba/tc/coupled+s+V_s=0.5+V_t=2+nl=1.5+t'=0.dat" s 1
```

By default it uses `coupled_get_res`'s `[x/2, 2x]` brackets; `NCS_BRACKET=f` solves in
`[x/f, x*f]` instead, which isolates one branch where singlet and triplet Tc are close.

## Independent check

`cpu/tests/check_singlet_tryplet.py` re-solves the gap + occupation equations with
scipy (adaptive quadrature with the Fermi points as breakpoints), sharing no code
with `NCS.h`:

```sh
python3 cpu/tests/check_singlet_tryplet.py singlet singlet.txt Vs=2 nl=1.1 tp=0
```
