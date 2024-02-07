#include <cmath>
#include <cuda_runtime.h>
#include <string>
#include <fstream>

using namespace std;

const double xk[]={0.9061798459386641,-0.9061798459386641,0.538469310105683,-0.5384693101056829,0.0};
const double ak[]={0.2369268850561876,0.2369268850561876,0.47862867049936647,0.47862867049936586,0.5688888888888889};

typedef double (*Under_Integral_Func)(double,double,double,double,double);

template <typename T>
struct cudaCallableFunctionPointer
{
public:
  cudaCallableFunctionPointer(T* f_)
  {
    T* host_ptr = (T*)malloc(sizeof(T));
    cudaMalloc((void**)&ptr, sizeof(T));

    cudaMemcpyFromSymbol(host_ptr, *f_, sizeof(T));
    cudaMemcpy(ptr, host_ptr, sizeof(T), cudaMemcpyHostToDevice);
    
    cudaFree(host_ptr);
  }

  ~cudaCallableFunctionPointer()
  {
    cudaFree(ptr);
  }

  T* ptr;
};




double occ_solve1D_zbr(double (*f)(double,double,double,int,int,double,double),double a, double b,double tol,double mi,double gamma0,int n,int k,double Vs,double tp);

double* sc_solve1D_zbr(double (*f)(double,double,double,int,int,double,double),double a, double b,double tol,double gamma0,int n,int k,double Vs,double sxc,double syc,double nl,double tp);

double coupled_solve1D_zbr(double (*f)(double,double,double,int,int,double,double,double),double a, double b,double tol,double gamma0,int n,int k,double Vs,double Vt,double sxc,double syc,double nl,double tp);

double sc_solve1D_new(double (*f)(double,double,double,int,int,double),double a,double tol,double mi,double gamma0,int n,int k,double Vs,double tp);

double coupled_solve1D_new(double (*f)(double,double,double,int,int,double,double,double),double a,double tol,double mi,double gamma0,int n,int k,double Vs,double Vt,double tp);

double* coupled_gamma0_solve1D_zbr(double* (*f)(double,double,double,double,double,double,double,int,int,double,double,double,double,bool),double a, double b,double fgmh,double fg,double tol,double Vs,double Vt,double nl,double tp,double sxc,double syc,int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet);

double sc_integrate1D_gl_gpu(Under_Integral_Func f,double a, double b, int n, int k,double Tc,double mi,double gamma0,double tp);

__device__ double device_function (double x, double a, double b, double c, double d);
__device__ double another_device_function (double x, double a, double b, double c, double d);

