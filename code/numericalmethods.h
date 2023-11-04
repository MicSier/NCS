#include <cmath>
#include <string>
#include <fstream>

using namespace std;


double occ_solve1D_zbr(double (*f)(double,double,double,double*,double*,int,int,double,double),double a, double b,double tol,double mi,double gamma0,double xk[],double ak[],int n,int k,double Vs,double tp);

double* sc_solve1D_zbr(double (*f)(double,double,double,double*,double*,int,int,double,double),double a, double b,double tol,double gamma0,double xk[],double ak[],int n,int k,double Vs,double sxc,double syc,double nl,double tp);

double coupled_solve1D_zbr(double (*f)(double,double,double,double*,double*,int,int,double,double,double),double a, double b,double tol,double gamma0,double xk[],double ak[],int n,int k,double Vs,double Vt,double sxc,double syc,double nl,double tp);

double sc_solve1D_new(double (*f)(double,double,double,double*,double*,int,int,double),double a,double tol,double mi,double gamma0,double xk[],double ak[],int n,int k,double Vs,double tp);

double coupled_solve1D_new(double (*f)(double,double,double,double*,double*,int,int,double,double,double),double a,double tol,double mi,double gamma0,double xk[],double ak[],int n,int k,double Vs,double Vt,double tp);

double* coupled_gamma0_solve1D_zbr(double* (*f)(double,double,double,double,double,double,double,double*,double*,int,int,double,double,double,double,bool),double a, double b,double fgmh,double fg,double tol,double Vs,double Vt,double nl,double tp,double sxc,double syc,double xk[],double ak[],int n,int k,double pts,double ptt,double pchs,double pcht,bool isfromsinglet);

double sc_integrate1D_gl_parallel(double (*f)(double,double,double,double,double),double a, double b, int n, int k,double* xk,double* ak,double Tc,double mi,double gamma0,double tp);
