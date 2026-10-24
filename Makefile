# Makefile

all: singlet tryplet coupled

singlet:
	nvcc examples/singlet/singlet.cu -o singlet -O3

tryplet:
	nvcc examples/tryplet/tryplet.cu -o tryplet -O3 

coupled:
	nvcc examples/coupled/coupled.cu -o coupled -O3 


# ---- CPU emulation (no GPU / no nvcc needed) --------------------------------
# Builds the same examples with g++ through cpu/cuda_runtime.h; see cpu/README.md.
CPU_CXX   ?= g++
CPU_FLAGS ?= -std=c++17 -O3 -fopenmp -x c++ -Icpu
CPU_EXAMPLES = singlet tryplet coupled IntegrationBenchmark

cpu: cpu-gen $(addsuffix _cpu,$(CPU_EXAMPLES)) test_kernels_cpu regress_rashba_cpu

cpu-gen:
	python3 cpu/cudaify.py build-cpu

singlet_cpu tryplet_cpu coupled_cpu: %_cpu: cpu-gen
	$(CPU_CXX) $(CPU_FLAGS) build-cpu/examples/$*/$*.cu -o $@

IntegrationBenchmark_cpu: cpu-gen
	$(CPU_CXX) $(CPU_FLAGS) build-cpu/examples/IntegrationBenchmark/integration_benchmark.cu -o $@

test_kernels_cpu: cpu-gen
	$(CPU_CXX) $(CPU_FLAGS) build-cpu/cpu/tests/test_kernels.cu -o $@

regress_rashba_cpu: cpu-gen
	$(CPU_CXX) $(CPU_FLAGS) build-cpu/cpu/tests/regress_rashba.cu -o $@

cpu-test: test_kernels_cpu
	./test_kernels_cpu

.PHONY: all singlet tryplet coupled cpu cpu-gen cpu-test
