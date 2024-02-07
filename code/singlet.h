#include <cuda_runtime.h>

double evth(double x);

double g(double x,double gamma0);

double ek(double x,double tp);

__device__ double sc_uifunc(double x,double Tc,double mi, double gamma0,double tp);
__device__ double singlet_uifunt(double x,double Tc,double mi,double gamma0,double tp);

double singlet_gap(double Tc,double mi,double gamma0,int n,int k,double Vs,double tp);

double sc_occ(double mi,double Tc,double gamma0,int n,int k,double nl,double tp);

double* singlet_get_res(double Vs,double nl,double tp,double gamma0,double xt,double yt,double xc,double yc,int n,int k);

double* sc_dos(double o,double g0,double tp);

double fs(double gamma0,double mi,double lambda);



