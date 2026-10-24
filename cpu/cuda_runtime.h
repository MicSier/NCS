// cpu/cuda_runtime.h -- minimal CUDA-on-CPU emulation for NCS
//
// Lets NCS.h and the examples build with g++/clang++ and run on a machine
// without an NVIDIA GPU. Kernels are executed unchanged:
//   * each block runs on one OS thread (blocks are spread over OpenMP threads),
//   * each CUDA thread runs on a lightweight fiber, and __syncthreads() is a real
//     barrier (all threads of the block finish phase k before any starts k+1),
//   * blocks whose first thread finishes without hitting a barrier are run as a
//     plain loop (fast path, e.g. integrateKernel).
// `kernel<<<g, b[, shmem]>>>(args)` launches and `extern __shared__` arrays are
// rewritten by cpu/cudaify.py, which generates the CPU build tree.
//
// Linux x86-64 uses a tiny hand-written context switch; other platforms fall
// back to <ucontext.h> (correct, but much slower for kernels with barriers).
#pragma once
#ifndef NCS_CPU_CUDA_RUNTIME_H
#define NCS_CPU_CUDA_RUNTIME_H

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <memory>

#if !(defined(__x86_64__) && defined(__linux__))
#include <ucontext.h>
#define NCS_CPU_UCONTEXT 1
#endif

#ifndef NCS_CPU_MAX_THREADS_PER_BLOCK
#define NCS_CPU_MAX_THREADS_PER_BLOCK 1024
#endif
#ifndef NCS_CPU_FIBER_STACK
#define NCS_CPU_FIBER_STACK (64 * 1024)
#endif

// ---- qualifiers -------------------------------------------------------------
#define __global__
#define __device__
#define __host__
#define __constant__
#define __forceinline__ inline
#define __restrict__
#define __shared__ static thread_local

// ---- vector types -----------------------------------------------------------
struct uint3 { unsigned int x, y, z; };
struct dim3 {
    unsigned int x, y, z;
    dim3(unsigned int x_ = 1, unsigned int y_ = 1, unsigned int z_ = 1) : x(x_), y(y_), z(z_) {}
};

inline thread_local uint3 threadIdx;
inline thread_local uint3 blockIdx;
inline thread_local dim3  blockDim;
inline thread_local dim3  gridDim;

// ---- runtime API subset -----------------------------------------------------
typedef int cudaError_t;
enum { cudaSuccess = 0, cudaErrorLaunchFailure = 4 };
enum cudaMemcpyKind { cudaMemcpyHostToHost, cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost,
                      cudaMemcpyDeviceToDevice, cudaMemcpyDefault };
enum cudaDeviceAttr { cudaDevAttrMaxThreadsPerBlock = 1 };
struct cudaDeviceProp { char name[256]; int maxThreadsPerBlock; int multiProcessorCount; int warpSize; };

inline thread_local cudaError_t ncs_cpu_last_error = cudaSuccess;

inline cudaError_t cudaMalloc(void** p, size_t sz) { *p = std::malloc(sz ? sz : 1); return *p ? cudaSuccess : 2; }
template <class T> inline cudaError_t cudaMalloc(T** p, size_t sz) { return cudaMalloc((void**)p, sz); }
inline cudaError_t cudaFree(void* p) { std::free(p); return cudaSuccess; }
inline cudaError_t cudaMemcpy(void* d, const void* s, size_t sz, cudaMemcpyKind) { std::memcpy(d, s, sz); return cudaSuccess; }
inline cudaError_t cudaMemset(void* d, int v, size_t sz) { std::memset(d, v, sz); return cudaSuccess; }
inline cudaError_t cudaGetDeviceCount(int* n) { *n = 1; return cudaSuccess; }
inline cudaError_t cudaDeviceSynchronize() { return cudaSuccess; }
inline cudaError_t cudaGetLastError() { cudaError_t e = ncs_cpu_last_error; ncs_cpu_last_error = cudaSuccess; return e; }
inline const char* cudaGetErrorString(cudaError_t e) { return e == cudaSuccess ? "no error" : "cpu emulation error"; }
inline cudaError_t cudaGetDeviceProperties(cudaDeviceProp* p, int) {
    std::memset(p, 0, sizeof(*p));
    std::snprintf(p->name, sizeof(p->name), "NCS CPU emulation");
    p->maxThreadsPerBlock = NCS_CPU_MAX_THREADS_PER_BLOCK; p->multiProcessorCount = 1; p->warpSize = 32;
    return cudaSuccess;
}
inline cudaError_t cudaDeviceGetAttribute(int* v, cudaDeviceAttr a, int) {
    *v = (a == cudaDevAttrMaxThreadsPerBlock) ? NCS_CPU_MAX_THREADS_PER_BLOCK : 0;
    return cudaSuccess;
}

// ---- fibers -----------------------------------------------------------------
namespace ncs_cpu {

#ifndef NCS_CPU_UCONTEXT
extern "C" void ncs_cpu_ctx_switch(void** save_sp, void* load_sp);
asm(R"(
    .text
    .weak ncs_cpu_ctx_switch
    .type ncs_cpu_ctx_switch,@function
ncs_cpu_ctx_switch:
    pushq %rbp
    pushq %rbx
    pushq %r12
    pushq %r13
    pushq %r14
    pushq %r15
    movq  %rsp, (%rdi)
    movq  %rsi, %rsp
    popq  %r15
    popq  %r14
    popq  %r13
    popq  %r12
    popq  %rbx
    popq  %rbp
    ret
    .size ncs_cpu_ctx_switch, .-ncs_cpu_ctx_switch
)");
#endif

struct Fiber {
    std::unique_ptr<char[]> stack;
#ifdef NCS_CPU_UCONTEXT
    ucontext_t ctx;
#else
    void* sp = nullptr;
#endif
    uint3 tid{};
    bool done = false;
    bool hit_barrier = false;
};

struct Scheduler {
#ifdef NCS_CPU_UCONTEXT
    ucontext_t main_ctx;
#else
    void* main_sp = nullptr;
#endif
    std::vector<Fiber> pool;
    Fiber* current = nullptr;
    void (*body)(void*) = nullptr;
    void* body_arg = nullptr;
    std::vector<char> dyn_shared;
};

inline thread_local Scheduler sched;

[[noreturn]] inline void fiber_entry() {
    Fiber* f = sched.current;
    threadIdx = f->tid;
    sched.body(sched.body_arg);
    f->done = true;
#ifdef NCS_CPU_UCONTEXT
    swapcontext(&f->ctx, &sched.main_ctx);
#else
    void* dummy;
    ncs_cpu_ctx_switch(&dummy, sched.main_sp);
#endif
    std::abort();  // a finished fiber is never resumed
}
#ifdef NCS_CPU_UCONTEXT
inline void fiber_entry_uc() { fiber_entry(); }
#endif

inline void fiber_init(Fiber& f, uint3 tid) {
    if (!f.stack) f.stack.reset(new char[NCS_CPU_FIBER_STACK]);
    f.tid = tid; f.done = false; f.hit_barrier = false;
#ifdef NCS_CPU_UCONTEXT
    getcontext(&f.ctx);
    f.ctx.uc_stack.ss_sp = f.stack.get();
    f.ctx.uc_stack.ss_size = NCS_CPU_FIBER_STACK;
    f.ctx.uc_link = nullptr;
    makecontext(&f.ctx, fiber_entry_uc, 0);
#else
    // Layout so that the pops + ret in ncs_cpu_ctx_switch "call" fiber_entry
    // with a SysV-conforming stack (rsp % 16 == 8 at function entry).
    uintptr_t top = (reinterpret_cast<uintptr_t>(f.stack.get()) + NCS_CPU_FIBER_STACK) & ~uintptr_t(15);
    void** s = reinterpret_cast<void**>(top - 16);
    s[0] = reinterpret_cast<void*>(&fiber_entry);        // return address
    s[1] = nullptr;                                       // fake caller frame
    f.sp = reinterpret_cast<void*>(top - 16 - 6 * sizeof(void*));
    std::memset(f.sp, 0, 6 * sizeof(void*));              // r15 r14 r13 r12 rbx rbp
#endif
}

inline void fiber_resume(Fiber& f) {
    sched.current = &f;
    threadIdx = f.tid;
#ifdef NCS_CPU_UCONTEXT
    swapcontext(&sched.main_ctx, &f.ctx);
#else
    ncs_cpu_ctx_switch(&sched.main_sp, f.sp);
#endif
    sched.current = nullptr;
}

inline void barrier() {
    Fiber* f = sched.current;
    if (!f) return;  // plain-loop mode: kernel had no barrier in thread 0
    f->hit_barrier = true;
#ifdef NCS_CPU_UCONTEXT
    swapcontext(&f->ctx, &sched.main_ctx);
#else
    ncs_cpu_ctx_switch(&f->sp, sched.main_sp);
#endif
    threadIdx = f->tid;
}

inline uint3 linear_to_3d(unsigned i, const dim3& d) {
    return uint3{ i % d.x, (i / d.x) % d.y, i / (d.x * d.y) };
}

// Runs one block. Thread 0 is probed on a fiber; if it never reaches a barrier,
// the remaining threads run as a plain loop, otherwise every thread gets a fiber
// and they are resumed round-robin, one barrier phase at a time.
inline void run_block(void (*body)(void*), void* arg, unsigned nthreads) {
    sched.body = body; sched.body_arg = arg;
    if (sched.pool.size() < nthreads) sched.pool.resize(nthreads);

    fiber_init(sched.pool[0], linear_to_3d(0, blockDim));
    fiber_resume(sched.pool[0]);
    if (sched.pool[0].done && !sched.pool[0].hit_barrier) {
        for (unsigned t = 1; t < nthreads; ++t) { threadIdx = linear_to_3d(t, blockDim); body(arg); }
        return;
    }
    // Thread 0 is parked at its first barrier: bring every other thread there too
    // before advancing anyone to the next phase.
    for (unsigned t = 1; t < nthreads; ++t) {
        fiber_init(sched.pool[t], linear_to_3d(t, blockDim));
        fiber_resume(sched.pool[t]);
    }
    bool pending = true;
    while (pending) {
        pending = false;
        for (unsigned t = 0; t < nthreads; ++t) {
            Fiber& f = sched.pool[t];
            if (f.done) continue;
            fiber_resume(f);
            if (!f.done) pending = true;
        }
    }
}

struct LaunchCfg { dim3 grid, block; size_t shmem; };

template <class F>
void launch(const LaunchCfg& cfg, F&& f) {
    const unsigned nthreads = cfg.block.x * cfg.block.y * cfg.block.z;
    const long nblocks = long(cfg.grid.x) * cfg.grid.y * cfg.grid.z;
    if (nthreads == 0 || nthreads > NCS_CPU_MAX_THREADS_PER_BLOCK || nblocks <= 0) {
        std::fprintf(stderr, "[ncs-cpu] invalid launch config: %ld blocks x %u threads\n", nblocks, nthreads);
        ncs_cpu_last_error = cudaErrorLaunchFailure;
        return;
    }
    using Fn = std::remove_reference_t<F>;
    auto tramp = [](void* p) { (*static_cast<Fn*>(p))(); };
#pragma omp parallel for schedule(dynamic)
    for (long b = 0; b < nblocks; ++b) {
        gridDim = cfg.grid; blockDim = cfg.block;
        dim3 g = cfg.grid;
        blockIdx = uint3{ unsigned(b % g.x), unsigned((b / g.x) % g.y), unsigned(b / (long(g.x) * g.y)) };
        if (sched.dyn_shared.size() < cfg.shmem) sched.dyn_shared.resize(cfg.shmem);
        run_block(+tramp, const_cast<void*>(static_cast<const void*>(&f)), nthreads);
    }
}

template <class T> inline T* dyn_shared() { return reinterpret_cast<T*>(sched.dyn_shared.data()); }

} // namespace ncs_cpu

inline void __syncthreads() { ncs_cpu::barrier(); }

inline ncs_cpu::LaunchCfg cpu_cfg(dim3 g, dim3 b, size_t shmem = 0, int /*stream*/ = 0) { return { g, b, shmem }; }

template <class... P, class... A>
inline void cpu_launch(const ncs_cpu::LaunchCfg& cfg, void (*kernel)(P...), A&&... args) {
    auto call = [&]() { kernel(args...); };
    ncs_cpu::launch(cfg, call);
}

#endif // NCS_CPU_CUDA_RUNTIME_H
