// Checks the CUDA kernels in NCS.h through the CPU emulation:
//  1. barrier semantics of the emulator itself (shared-memory reversal),
//  2. sc_reduce_sum_gpu (reduce_custom + its host loop) on many sizes,
//  3. reduceShared (currently unused by the library),
//  4. sc_integrate1D_gl_gpu for every integrand vs. a serial long-double
//     composite Gauss-Legendre sum over the same panels.
// Build: python3 cpu/cudaify.py && g++ -std=c++17 -O2 -fopenmp -x c++ -Icpu build-cpu/cpu/tests/test_kernels.cu
#define NCS_IMPLEMENTATION
#include "../../NCS.h"

#include <cstdio>
#include <random>

static int failures = 0;
static void check(bool ok, const char* what, double got, double want) {
    double rel = std::fabs(got - want) / std::max(1.0, std::fabs(want));
    std::printf("%-58s %s  got % .17g  want % .17g  rel %.2e\n", what, ok ? "ok  " : "FAIL", got, want, rel);
    if (!ok) ++failures;
}

__global__ void reverse_shared(const int* in, int* out, int n) {
    __shared__ int buf[1024];
    int t = threadIdx.x;
    if (t < n) buf[t] = in[t];
    __syncthreads();
    if (t < n) out[t] = buf[n - 1 - t] * 10 + (t > 0 ? buf[t - 1] % 10 : 0);
}

static double reduce_custom_host_loop(std::vector<double> v) {
    // The library's own host loop around reduce_custom
    int n = (int)v.size();
    double *d_data, *d_scratch;
    cudaMalloc((void**)&d_data, n * sizeof(double));
    cudaMalloc((void**)&d_scratch, n * sizeof(double));
    cudaMemcpy(d_data, v.data(), n * sizeof(double), cudaMemcpyHostToDevice);
    double r = sc_reduce_sum_gpu(d_data, d_scratch, n);
    cudaFree(d_data); cudaFree(d_scratch);
    return r;
}

static double reduce_shared_full(std::vector<double> v, int bs) {
    int n = (int)v.size();
    double *a, *b;
    cudaMalloc((void**)&a, n * sizeof(double));
    cudaMalloc((void**)&b, n * sizeof(double));
    cudaMemcpy(a, v.data(), n * sizeof(double), cudaMemcpyHostToDevice);
    while (n > 1) {
        int blocks = (n + 2 * bs - 1) / (2 * bs);
        reduceShared<<<blocks, bs, bs * sizeof(double)>>>(a, b, n);
        std::swap(a, b);
        n = blocks;
    }
    double r;
    cudaMemcpy(&r, a, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(a); cudaFree(b);
    return r;
}

static long double serial_gl(Under_Integral_Func f, double a, double b, double Tc, double mi, double g0, double tp) {
    double h = (b - a) / n;
    long double s = 0.0L;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < k; ++j) {
            double z = ((2.0 * a + h * (2.0 * i + 1.0)) - h * xk[j]) * 0.5;
            s += (long double)ak[j] * f(z, Tc, mi, g0, tp);
        }
    return s * h * 0.5L;
}

int main() {
    // 1. emulator barrier sanity
    {
        const int N = 1000;
        std::vector<int> in(N), out(N);
        for (int i = 0; i < N; ++i) in[i] = i;
        int *d_in, *d_out;
        cudaMalloc((void**)&d_in, N * sizeof(int)); cudaMalloc((void**)&d_out, N * sizeof(int));
        cudaMemcpy(d_in, in.data(), N * sizeof(int), cudaMemcpyHostToDevice);
        reverse_shared<<<1, 1024>>>(d_in, d_out, N);
        cudaMemcpy(out.data(), d_out, N * sizeof(int), cudaMemcpyDeviceToHost);
        bool ok = true;
        for (int t = 0; t < N; ++t) ok &= out[t] == (N - 1 - t) * 10 + (t > 0 ? (t - 1) % 10 : 0);
        check(ok, "emulator: __syncthreads orders shared-memory phases", ok, 1);
        cudaFree(d_in); cudaFree(d_out);
    }

    // 2./3. reductions
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> U(-1.0, 2.0);
    for (int sz : {1, 2, 3, 5, 1023, 1024, 2047, 2048, 2049, 4097, 10000, 100000}) {
        std::vector<double> v(sz);
        long double exact = 0;
        for (auto& x : v) { x = U(rng); exact += x; }
        char label[96];
        double r1 = reduce_custom_host_loop(v);
        std::snprintf(label, sizeof label, "sc_reduce_sum_gpu n=%d", sz);
        check(std::fabs(r1 - (double)exact) <= 1e-12 * std::max(1.0, (double)std::fabs(exact)), label, r1, (double)exact);
        double r2 = reduce_shared_full(v, 256);
        std::snprintf(label, sizeof label, "reduceShared    n=%d (256 threads)", sz);
        check(std::fabs(r2 - (double)exact) <= 1e-12 * std::max(1.0, (double)std::fabs(exact)), label, r2, (double)exact);
    }

    // 4. full GPU integration pipeline vs serial quadrature
    struct { UI_Func_Index i; Under_Integral_Func f; const char* name; } fs_[] = {
        { i_sc_uifunc, sc_uifunc, "sc_uifunc" }, { i_singlet_uifunt, singlet_uifunt, "singlet_uifunt" },
        { i_tryplet_uifunt, tryplet_uifunt, "tryplet_uifunt" }, { i_coupled_fp, coupled_fp, "coupled_fp" },
        { i_coupled_fmdk, coupled_fmdk, "coupled_fmdk" }, { i_coupled_fpdk2, coupled_fpdk2, "coupled_fpdk2" } };
    double params[][4] = { { 0.1, 0.1, 0.1, 0.0 }, { 0.05, -0.7, 0.5, 0.0 }, { 0.02, 0.3, 1.5, 0.2 } };
    for (auto& F : fs_)
        for (auto& p : params) {
            double got = sc_integrate1D_gl_gpu(F.i, 0.0, pi, p[0], p[1], p[2], p[3]);
            long double want = serial_gl(F.f, 0.0, pi, p[0], p[1], p[2], p[3]);
            char label[96];
            std::snprintf(label, sizeof label, "integrate %-15s Tc=%.2f mu=%.2f g0=%.1f", F.name, p[0], p[1], p[2]);
            check(std::fabs(got - (double)want) <= 1e-12 * std::max(1.0, (double)std::fabs(want)), label, got, (double)want);
        }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures, failures == 1 ? "" : "s");
    return failures != 0;
}
